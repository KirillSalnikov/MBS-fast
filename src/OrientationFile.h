#pragma once

#include <string>
#include <vector>

struct EulerOrientationRadians
{
    double beta;
    double gamma;
};

struct WeightedEulerOrientationRadians
{
    double beta, gamma, alpha, weight;
};
// Explicit Haar cubature: BETA_DEG GAMMA_DEG ALPHA_DEG WEIGHT, sum(weight)=1.
std::vector<WeightedEulerOrientationRadians> ReadSO3RuleDegrees(const std::string& path);

// Orientation files are a user-facing format and therefore store angles in
// degrees. The returned values are ready for Particle::Rotate(), which expects
// radians.
std::vector<EulerOrientationRadians> ReadOrientationFileDegrees(
    const std::string &path);
