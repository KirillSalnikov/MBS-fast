#pragma once

#include "geometry_lib.h"

#include <algorithm>
#include <cfloat>
#include <cmath>

// All topology predicates and generated vertices use the same double storage.
// Tolerances are relative to the physical scale of the corresponding length,
// area, or scalar product instead of being absolute model-space constants.
// Keep roughly 40 significant binary bits for topology decisions.  This is
// well below physical/model resolution, but above accumulated roundoff from
// repeated projection and clipping in double precision.
constexpr double GEOMETRY_RELATIVE_EPSILON = 0x1p-40;
constexpr double EPS_PROJECTION = GEOMETRY_RELATIVE_EPSILON;

inline double geometry_length_tolerance(double scale)
{
    return GEOMETRY_RELATIVE_EPSILON
        * std::max(std::fabs(scale), DBL_MIN);
}

inline double geometry_area_tolerance(double scale)
{
    const double length = std::max(std::fabs(scale), DBL_MIN);
    return GEOMETRY_RELATIVE_EPSILON * length * length;
}

inline double geometry_depth_tolerance(double scale)
{
    // A depth classification combines a rotated plane parameter, a dot
    // product and a projection.  Allow their accumulated rounding while still
    // retaining about 37 significant binary bits of geometric separation.
    return 8.0 * geometry_length_tolerance(scale);
}

inline double geometry_parallel_tolerance(double productScale)
{
    return GEOMETRY_RELATIVE_EPSILON
        * std::max(std::fabs(productScale), DBL_MIN);
}

inline double geometry_length3(double x, double y, double z)
{
    return std::sqrt(x*x + y*y + z*z);
}

inline bool geometry_points_distinct(const Point3f &a, const Point3f &b,
                                     double geometryScale)
{
    const double dx = a.cx - b.cx;
    const double dy = a.cy - b.cy;
    const double dz = a.cz - b.cz;
    return geometry_length3(dx, dy, dz)
        > geometry_length_tolerance(geometryScale);
}

// Clip a convex polygon by an affine scalar field sampled at its vertices.
// The scalar is interpolated along every edge, so the generated crossing lies
// on the exact zero level set.  Values within tolerance are snapped to zero;
// this keeps topology decisions scale invariant without displacing the plane.
inline int clip_polygon_by_nonnegative_scalar(const Point3f *input,
                                               const double *scalar,
                                               int inputSize,
                                               double tolerance,
                                               Point3f *output)
{
    if (inputSize <= 0)
        return 0;
    if (inputSize > MAX_VERTEX_NUM)
        throw std::runtime_error("invalid polygon size in half-space clipping");

    const auto snapped = [tolerance](double value) {
        return std::fabs(value) <= tolerance ? 0.0 : value;
    };
    const auto append = [](Point3f *vertices, int &size,
                           const Point3f &point) {
        if (size >= MAX_VERTEX_NUM)
            throw std::runtime_error(
                "half-space clipping exceeds the 64-vertex geometry limit");
        vertices[size++] = point;
    };

    int outputSize = 0;
    Point3f start = input[inputSize - 1];
    double startValue = snapped(scalar[inputSize - 1]);
    bool startInside = startValue >= 0.0;

    for (int index = 0; index < inputSize; ++index)
    {
        const Point3f end = input[index];
        const double endValue = snapped(scalar[index]);
        const bool endInside = endValue >= 0.0;

        if (startInside != endInside)
        {
            const double denominator = startValue - endValue;
            if (std::fabs(denominator) > DBL_MIN)
            {
                const double fraction = std::max(0.0, std::min(1.0,
                    startValue / denominator));
                append(output, outputSize,
                       start + (end - start)*fraction);
            }
        }
        if (endInside)
            append(output, outputSize, end);

        start = end;
        startValue = endValue;
        startInside = endInside;
    }
    return outputSize;
}

// ProjectToFacetPlane maps p to p - t*direction, where
// t = plane(p)/(direction.planeNormal).  In the beam-clipping callers the
// projection direction points from a target facet back to the source aperture,
// so only t <= 0 belongs to the forward ray half-space.  Clip the polygon
// before projection: accepting the complete polygon when only one vertex is
// forward lets its behind-plane part shadow the aperture.
inline int clip_polygon_by_nonpositive_projection_parameter(
    const Point3f *input, int inputSize,
    const Point3f &planeNormal, const Point3f &projectionDirection,
    double depthTolerance, Point3f *output)
{
    if (inputSize <= 0)
        return 0;
    if (inputSize > MAX_VERTEX_NUM)
        throw std::runtime_error("invalid polygon size in projection-depth clipping");

    const double denominator = DotProduct(projectionDirection, planeNormal);
    const double denominatorScale = geometry_length3(
        projectionDirection.cx, projectionDirection.cy,
        projectionDirection.cz) * geometry_length3(
            planeNormal.cx, planeNormal.cy, planeNormal.cz);
    if (denominatorScale <= DBL_MIN
        || std::fabs(denominator)
            <= geometry_parallel_tolerance(denominatorScale))
    {
        return 0;
    }

    const double denominatorSign = denominator >= 0.0 ? 1.0 : -1.0;
    double forwardScalar[MAX_VERTEX_NUM];
    for (int i = 0; i < inputSize; ++i)
    {
        const double planeValue = DotProduct(input[i], planeNormal)
            + planeNormal.d_param;
        // -t >= 0, without dividing by a possibly small denominator.
        forwardScalar[i] = -planeValue * denominatorSign;
    }
    return clip_polygon_by_nonnegative_scalar(
        input, forwardScalar, inputSize,
        depthTolerance * std::fabs(denominator), output);
}

inline bool is_inside_i(const Point3f &x, const Point3f &p1,
                        const Point3f &p2, const Point3f &normal)
{
    const double ex = p2.cx - p1.cx;
    const double ey = p2.cy - p1.cy;
    const double ez = p2.cz - p1.cz;
    const double qx = x.cx - p1.cx;
    const double qy = x.cy - p1.cy;
    const double qz = x.cz - p1.cz;

    const double crossX = ey*qz - ez*qy;
    const double crossY = ez*qx - ex*qz;
    const double crossZ = ex*qy - ey*qx;
    const double side = crossX*normal.cx
        + crossY*normal.cy + crossZ*normal.cz;
    const double edgeLength = geometry_length3(ex, ey, ez);
    const double offsetLength = geometry_length3(qx, qy, qz);
    const double normalLength = Length(normal);
    const double areaScale = edgeLength
        * std::max(edgeLength, offsetLength) * normalLength;
    const double tolerance = GEOMETRY_RELATIVE_EPSILON * areaScale;
    return side >= -tolerance;
}

inline bool is_layOnLine_i(const Point3f &x, const Point3f &a,
                           const Point3f &b)
{
    const Point3f ab = b - a;
    const Point3f ax = x - a;
    const double length2 = Norm(ab);
    if (length2 <= DBL_MIN)
        return false;

    const double t = DotProduct(ax, ab) / length2;
    return t >= -EPS_PROJECTION && t <= 1.0 + EPS_PROJECTION;
}

inline Point3f intersect_geometry_lines(const Point3f &a1,
                                        const Point3f &b1,
                                        const Point3f &va,
                                        const Point3f &vb,
                                        const Point3f &normal,
                                        bool &ok)
{
    const Point3f transverse = CrossProduct(vb, normal);
    const double denominator = DotProduct(va, transverse);
    const double denominatorScale = Length(va)*Length(transverse);
    if (denominatorScale <= DBL_MIN
        || std::fabs(denominator)
            <= geometry_parallel_tolerance(denominatorScale))
    {
        ok = false;
        return Point3f();
    }

    const double t = DotProduct(a1 - b1, transverse) / denominator;
    const Point3f result = a1 - va*t;
    if (!std::isfinite(result.cx) || !std::isfinite(result.cy)
        || !std::isfinite(result.cz))
    {
        ok = false;
        return Point3f();
    }

    ok = true;
    return result;
}

inline Point3f intersect_i(const Point3f &a1, const Point3f &a2,
                           const Point3f &b1, const Point3f &b2,
                           const Point3f &normal, bool &ok)
{
    return intersect_geometry_lines(a1, b1, a2 - a1, b2 - b1, normal, ok);
}

inline Point3f intersect_iv(const Point3f &a1, const Point3f &b1,
                            const Point3f &va, const Point3f &vb,
                            const Point3f &normal, bool &ok)
{
    return intersect_geometry_lines(a1, b1, va, vb, normal, ok);
}

void computeIntersection(const Point3f &s, const Point3f &e,
                         const Point3f &p1, const Point3f &p2,
                         const Point3f &normal, Point3f &x);
