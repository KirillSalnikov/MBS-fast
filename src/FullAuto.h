#pragma once
#include <array>
#include <vector>
class ArgPP;
struct RunConfig;

namespace FullAuto {
typedef std::array<double,16> Mueller;
typedef std::vector<Mueller> Matrix;
typedef std::vector<Matrix> Samples;
Matrix CubicInterpolate(const std::vector<double>& x, const Matrix& y,
                        const std::vector<double>& query);
struct Statistics {
    Matrix mean, interval;
    double maxInterval = 0;
};
Statistics Estimate(const Samples& samples);
double Difference(const Matrix& a, const Matrix& b, const Matrix& reference);
std::vector<int> Allocate(const std::vector<Samples>& levels,
                          const std::vector<double>& costs,
                          const std::vector<int>& counts,
                          const Matrix& reference, double tolerance);
int Run(const ArgPP& args, const RunConfig& config, int argc, const char* argv[]);
}
