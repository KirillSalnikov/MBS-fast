#pragma once
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

// Scheduling only. All acceptance tests use independently sampled matrices.
namespace FullAutoAuditPlan {
using Budgets=std::vector<std::vector<int>>;
struct StencilGroup {std::set<size_t> rows;std::vector<size_t> intervals;};
using StencilGroups=std::map<std::vector<int>,StencilGroup>;
inline void CheckIds(const std::vector<size_t>& ids,size_t size){
    if(ids.empty()||!std::is_sorted(ids.begin(),ids.end())||std::adjacent_find(ids.begin(),ids.end())!=ids.end()||ids.back()>=size)
        throw std::invalid_argument("invalid independent audit angular indices");
}
inline void Extend(Budgets& counts,const std::vector<size_t>& production,const std::vector<size_t>& audit){
    if(counts.empty())throw std::invalid_argument("empty independent audit budgets");
    CheckIds(production,counts[0].size());CheckIds(audit,counts[0].size());
    for(const auto& level:counts)if(level.size()!=counts[0].size())throw std::invalid_argument("inconsistent audit budget dimensions");
    // Read only the production neighbours. New quarter guards must not spread
    // a large budget recursively across a sequence of newly inserted nodes.
    const Budgets before=counts;
    for(size_t id:audit)if(!std::binary_search(production.begin(),production.end(),id)){
        const auto upper=std::upper_bound(production.begin(),production.end(),id);
        const size_t right=upper==production.end()?production.back():*upper;
        const size_t left=upper==production.begin()?production.front():*(upper-1);
        for(size_t l=0;l<counts.size();++l)counts[l][id]=std::max(before[l][id],std::max(before[l][left],before[l][right]));
    }
}
inline Budgets Scale(const Budgets& production,const std::vector<size_t>& ids,double fraction,
                     int initial,int minimum,bool local){
    if(production.empty()||!(fraction>0&&fraction<=1)||initial<=0||minimum<=0)
        throw std::invalid_argument("invalid independent audit allocation");
    CheckIds(ids,production[0].size());Budgets result=production;
    for(size_t l=0;l<production.size();++l){
        if(production[l].size()!=production[0].size())throw std::invalid_argument("inconsistent audit budget dimensions");
        int peak=0;for(int n:production[l]){
            if(n<=0||(n&(n-1)))throw std::invalid_argument("audit budget must be a positive power of two");
            peak=std::max(peak,n);
        }
        for(size_t id:ids){
            const int n=local?production[l][id]:peak;
            const int floor=std::min(n,l?minimum:initial);
            int scaled=1;while(scaled<n*fraction&&scaled<n)scaled*=2;
            result[l][id]=std::min(n,std::max(floor,scaled));
        }
    }return result;
}
inline StencilGroups Stencils(const std::vector<size_t>& knots,const Budgets& counts,const std::vector<size_t>& ids){
    if(counts.empty()||knots.size()<2)throw std::invalid_argument("invalid common-prefix audit stencil");
    CheckIds(knots,counts[0].size());CheckIds(ids,counts[0].size());
    for(size_t k:knots)if(!std::binary_search(ids.begin(),ids.end(),k))throw std::invalid_argument("missing audit support node");
    StencilGroups groups;const size_t width=std::min(size_t(4),knots.size());
    for(size_t interval=0;interval+1<knots.size();++interval){
        const size_t left=std::min(interval?interval-1:0,knots.size()-width),right=std::min(interval,knots.size()-width);
        std::set<size_t> nodes;for(size_t k=left;k<right+width;++k)nodes.insert(knots[k]);
        const auto begin=std::lower_bound(ids.begin(),ids.end(),knots[interval]);
        const auto end=std::upper_bound(begin,ids.end(),knots[interval+1]);nodes.insert(begin,end);
        std::vector<int> budget(counts.size(),0);
        for(size_t l=0;l<counts.size();++l)for(size_t id:nodes)budget[l]=std::max(budget[l],counts[l][id]);
        auto& group=groups[budget];group.rows.insert(nodes.begin(),nodes.end());group.intervals.push_back(interval);
    }return groups;
}
}
