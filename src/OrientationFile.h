#pragma once

#include <string>
#include <vector>

struct EulerOrientationRadians
{
    double beta;
    double gamma;
};

// Orientation files are a user-facing format and therefore store angles in
// degrees. The returned values are ready for Particle::Rotate(), which expects
// radians.
std::vector<EulerOrientationRadians> ReadOrientationFileDegrees(
    const std::string &path);
