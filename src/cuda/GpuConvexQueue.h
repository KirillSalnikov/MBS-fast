#pragma once

#include <vector>
#include <cstddef>
#include <functional>

// Experimental resident convex tracer. Geometry and optics always use FP64,
// independently of the diffraction storage precision. Limits cause a complete
// CPU retrace of the affected orientation, never a truncated beam set.
namespace GpuConvexQueue {
constexpr int MaxVertices = 24;
constexpr int MaxFacets = 16;
constexpr int StackCapacity = 256;
struct Vec { double x, y, z; };
struct Cx { double re, im; };
struct Beam {
    Vec vertices[MaxVertices], direction, basis;
    Cx j[4];
    double path, front;
    int count, facet, depth;
    unsigned char track[20];
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

// Device/allocation failures return false, with no accepted partial results.
bool Trace(const std::vector<Input>& inputs, std::vector<Output>& output,
           std::vector<size_t>& order, std::vector<Status>& status, Metrics& metrics, int threads,
           const std::vector<Output>* external = nullptr,
           const std::function<bool(const Input*,const Output*,size_t)>& consumeDevice = {});
}
