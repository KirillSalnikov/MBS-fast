#pragma once

#include "matrix.hpp"

#include <cfloat>
#include <cmath>

namespace PoleMueller
{
inline double AngularTolerance()
{
    return 256.0 * DBL_EPSILON;
}

inline bool IsForward(double theta)
{
    return std::fabs(theta) <= AngularTolerance();
}

inline bool IsBackward(double theta)
{
    return std::fabs(theta - M_PI) <= AngularTolerance();
}

inline void ApplyForward(matrix &m)
{
    const double M00 = m[0][0];
    const double M11 = 0.5 * (m[1][1] + m[2][2]);
    const double M33 = m[3][3];
    const double M03 = m[0][3];

    m.Fill(0.0);
    m[0][0] = M00;
    m[1][1] = M11;
    m[2][2] = M11;
    m[3][3] = M33;
    m[0][3] = M03;
    m[3][0] = M03;
}

inline void ApplyBackward(matrix &m)
{
    const double M00 = m[0][0];
    const double M11 = 0.5 * (m[1][1] - m[2][2]);
    const double M33 = m[3][3];
    const double M03 = m[0][3];

    m.Fill(0.0);
    m[0][0] = M00;
    m[1][1] = M11;
    m[2][2] = -M11;
    m[3][3] = M33;
    m[0][3] = M03;
    m[3][0] = M03;
}
}
