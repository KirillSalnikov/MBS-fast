#pragma once

#include <vector>
#include <cstddef>
#include <functional>
#include <string>

// Experimental resident convex tracer. Geometry and optics always use FP64,
// independently of the diffraction storage precision. Limits cause a complete
// CPU retrace of the affected orientation, never a truncated beam set.
namespace GpuConvexQueue {
constexpr int MaxVertices = 24;
// Droxtals have 20 faces. Keep a bounded input while covering them and other
// small convex polyhedra; larger shapes still take the complete CPU path.
constexpr int MaxFacets = 32;
constexpr int LegacyMaxDepth = 18;
constexpr int MaxDepth = 22;
constexpr int TrackCapacity = MaxDepth + 2;
// Depth-first traversal adds at most faceCount-1 pending children per level.
// Cover the deepest opt-in reference plus the initial face seeds.
constexpr int StackCapacity = MaxFacets * (MaxDepth + 2);
constexpr int LegacyStackCapacity = MaxFacets * (LegacyMaxDepth + 2);
struct Vec { double x, y, z; };
struct Cx { double re, im; };
struct Beam {
    Vec vertices[MaxVertices], direction, basis;
    Cx j[4];
    double path, front;
    int count, facet, depth;
    unsigned char track[TrackCapacity];
};
struct Face {
    Vec vertices[MaxVertices], inward, outward;
    double plane;
    int count;
};
struct Input {
    Face faces[MaxFacets];
    Beam seeds[MaxFacets];
    double ri, lengthTolerance, areaTolerance, depthTolerance;
    double reference[3], cutoff[3];
    int faceCount, seedCount, maxDepth;
};
struct Output { Beam beam; int orientation, order; };
struct Status {
    int failure, count;
    unsigned long long evaluated, rejected, rejectedJones, rejectedArea, rejectedImportance;
};
struct Metrics { double kernelSeconds = 0, wallSeconds = 0, downloadSeconds = 0, sortSeconds = 0; size_t bytes = 0; };

// Capacity errors must not hide an invalid polygon or a nonfinite beam: only
// pure capacity failures may be retried with larger temporary buffers.
#ifdef __CUDACC__
__host__ __device__
#endif
inline int CombineFailure(int a,int b){
    if(a==4||b==4)return 4;
    if(a==1||b==1)return 1;
    return a>b?a:b;
}

// Device/allocation failures return false, with no accepted partial results.
void SetCacheLimit(size_t bytes);
// Boundary-only cleanup; false preserves failed-free ownership and blocks reuse.
bool ReleaseWorkspace();
struct ChunkLimit {int count=1;size_t budgetBytes=0,reservedBytes=0;bool queried=false,fits=false;};
// The actual trace rechecks available memory; this is a conservative scheduling
// bound for one resident packet, including both permitted capacity retries.
ChunkLimit ResidentChunkLimit(int requested,int faces,int maximumDepth,bool allowWide=false);
bool Trace(const std::vector<Input>& inputs, std::vector<Output>& output,
           std::vector<size_t>& order, std::vector<Status>& status, Metrics& metrics, int threads,
           const std::vector<Output>* external = nullptr,
           const std::function<bool(const Input*,const Output*,size_t)>& consumeDevice = {},
           const std::string& capacityIdentity = std::string());
}
