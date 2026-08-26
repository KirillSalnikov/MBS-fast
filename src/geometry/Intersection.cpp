#include "Intersection.h"

void computeIntersection(const Point3f &s, const Point3f &e,
                         const Point3f &p1, const Point3f &p2,
                         const Point3f &normal, Point3f &x)
{
    const Point3f segment = e - s;
    const Point3f cuttingPlaneNormal = CrossProduct(p2 - p1, normal);
    const double denominator = DotProduct(segment, cuttingPlaneNormal);
    const double denominatorScale = std::sqrt(
        DotProduct(segment, segment)
        * DotProduct(cuttingPlaneNormal, cuttingPlaneNormal));
    if (denominatorScale <= DBL_MIN
        || std::fabs(denominator)
            <= geometry_parallel_tolerance(denominatorScale))
    {
        x = s;
        return;
    }

    const double t = DotProduct(p1 - s, cuttingPlaneNormal) / denominator;
    x = s + segment*t;
}
