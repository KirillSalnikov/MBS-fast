#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>
#include <map>
#include <string>

namespace FullAutoScheduling {
// A startup work estimate, never an accuracy certificate. Incoherent PO
// retains aperture diffraction and therefore is not scale-invariant like GO.
// Bound the pilot so a large size does not spend the final sample budget
// before the measured eight-replica allocator has any evidence.
inline double IncoherentInitial(double diameter,double wavelength,double imaginaryIndex) {
    const double pi=3.14159265358979323846;
    const double size=pi*diameter/wavelength;
    const double optical=4*pi*std::fabs(imaginaryIndex)*diameter/wavelength;
    return std::min(8192.,std::max(512.,std::max(size,512*(1+std::sqrt(std::max(0.,optical))))));
}
inline int CorrectionProbeCount(int current,int proposed,int maximum) {
    const int doubled=current>maximum/2?maximum:2*current;
    const int cap=current>maximum/4?maximum:4*current;
    return std::max(doubled,std::min(cap,proposed));
}
// Schedule a preparation recipe's angular packets consecutively. Large
// groups still start first, but eight seeds no longer evict one another's
// angle-independent beams between subdivisions of the same request.
struct NativePacket {std::string preparation;double cost;};
inline std::vector<size_t> NativePacketOrder(const std::vector<NativePacket>& packets) {
    struct Group {std::vector<size_t> ids;double cost=0;};
    std::vector<Group> groups;std::map<std::string,size_t> index;
    for(size_t i=0;i<packets.size();++i){
        const auto inserted=index.emplace(packets[i].preparation,groups.size());
        if(inserted.second)groups.emplace_back();
        auto& group=groups[inserted.first->second];group.ids.push_back(i);group.cost+=packets[i].cost;
    }
    std::stable_sort(groups.begin(),groups.end(),[](const Group& a,const Group& b){return a.cost>b.cost;});
    std::vector<size_t> order;
    for(auto& group:groups){
        std::stable_sort(group.ids.begin(),group.ids.end(),[&](size_t a,size_t b){return packets[a].cost>packets[b].cost;});
        order.insert(order.end(),group.ids.begin(),group.ids.end());
    }
    return order;
}
// An unbounded ratio interval is a sampling failure, never a work estimate.
inline bool Resolved(double error) { return std::isfinite(error)&&error>=0&&error<1e50; }
inline bool PreferNextDepth(double mean,double interval,double budget,int attempt,bool terminal) {
    // Scheduling heuristic only. After two paired pilots, a mean above the
    // budget cannot be rescued by shrinking its interval. If the interval is
    // over twice the remaining gap, the usual N^-1/2 forecast needs over 4x
    // more samples. Try a deeper candidate first. NEVER use this to accept a
    // depth; the terminal candidate retains its complete sampling search.
    if(terminal||attempt<1||!(budget>0))return false;
    return !Resolved(mean)||!Resolved(interval)||mean>=budget||interval>2*(budget-mean);
}
inline int PilotLimit(int count,int maximum,bool fixed) {
    if(fixed)return count;
    return std::max(count,std::min(maximum,count>4096/8?4096:count*8));
}
inline int NextPilot(int count,int limit) { return count>=limit?count:std::min(limit,count>limit/2?limit:2*count); }
inline std::vector<std::vector<double>> AngleBatches(const std::vector<double>& angles,size_t cap) {
    if(angles.empty())return {};
    if(!cap||angles.size()<=cap)return {angles};
    cap=std::max(size_t(2),cap);
    // Balance the last two packets: never emit a singleton or exceed the cap.
    const size_t packets=(angles.size()+cap-1)/cap;
    const size_t width=angles.size()/packets,extra=angles.size()%packets;
    std::vector<std::vector<double>> result;size_t begin=0;
    for(size_t i=0;i<packets;++i){const size_t end=begin+width+(i<extra);
        result.emplace_back(angles.begin()+begin,angles.begin()+end);begin=end;}
    return result;
}
inline int ShortPacketCandidate(size_t rows,bool fp32) {
    // Very narrow angular packets cannot hide long producer bursts. This is
    // only a candidate: the controller still checks two independent pilots
    // for all16/ratio agreement AND a measured reduction in processing time.
    return fp32||rows<=8?512:1024;
}
inline bool ShortPacket(size_t rows,int phi,int count,int chunk,bool explicitChunk,int candidateChunk=1024) {
    return !explicitChunk&&chunk!=candidateChunk&&count>=2048&&rows>0&&rows<=256&&phi>0&&rows<=size_t(1024/phi);
}
inline bool RetryPhi(int oldPeak,int peak,int attempts) {
    return oldPeak>0&&attempts<3&&peak/oldPeak>=4;
}
}
