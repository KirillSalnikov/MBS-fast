#pragma once
#include <array>
#include <cstddef>
#include <utility>
#include <vector>
class ArgPP;
struct RunConfig;

namespace FullAuto {
typedef std::array<double,16> Mueller;
typedef std::vector<Mueller> Matrix;
typedef std::vector<Matrix> Samples;
Matrix CubicInterpolate(const std::vector<double>& x, const Matrix& y,
                        const std::vector<double>& query);
Matrix LocalCubicInterpolate(const std::vector<double>& x,const Matrix& y,
                            const std::vector<double>& query);
// Sparse weights of the very same local Hermite interpolant (up to 5 nodes).
std::vector<std::pair<size_t,double>> LocalCubicWeights(const std::vector<double>& x,double query);
struct Statistics {
    Matrix mean, interval, ratioInterval;
    double maxInterval = 0, maxRatioInterval = 0;
};
Statistics Estimate(const Samples& samples);
// Numerical equivalence of geometry-verified rotations. This audits paired
// blocks, not convergence of the orientation average; no Fieller gate.
double RotationEquivalenceBound(const Samples& candidate,const Samples& reference,bool ratios=true);
double Difference(const Matrix& a, const Matrix& b, const Matrix& reference);
std::vector<int> Allocate(const std::vector<Samples>& levels,
                          const std::vector<double>& costs,
                          const std::vector<int>& counts,
                          const Matrix& reference, double tolerance);
std::vector<int> PrioritizeAllocation(const std::vector<Samples>& levels,
                          const std::vector<double>& costs,
                          const std::vector<int>& counts,
                          const std::vector<int>& proposed,
                          double tolerance, bool ratios=true,bool covariance=true);
std::vector<int> AllocateCovariant(const std::vector<Samples>& levels,
                          const std::vector<double>& costs,const std::vector<int>& counts,
                          const Matrix& reference,double tolerance,bool ratios,int maximum,double decisionTolerance=0);
int Run(const ArgPP& args, const RunConfig& config, int argc, const char* argv[],
        double maximumEdge=0,double maximumDiameter=0,double particleVolume=0);
}
