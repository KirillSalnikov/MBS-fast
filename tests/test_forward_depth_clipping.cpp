#include "geometry/Intersection.h"
#include "Splitting.h"

#include <cassert>
#include <cmath>

namespace
{
void CheckScale(double scale)
{
    const Point3f polygon[] = {
        Point3f(-2*scale, -scale, 0),
        Point3f( 2*scale, -scale, 0),
        Point3f( 2*scale,  scale, 0),
        Point3f(-2*scale,  scale, 0)
    };
    // An affine depth field whose zero line crosses two edges without passing
    // through an existing vertex.  The retained rectangle has half the area.
    const double depth[] = {-2*scale, 2*scale, 2*scale, -2*scale};
    Point3f output[MAX_VERTEX_NUM];
    const int size = clip_polygon_by_nonnegative_scalar(
        polygon, depth, 4, geometry_length_tolerance(4*scale), output);
    assert(size == 4);

    double twiceArea = 0.0;
    for (int i = 0; i < size; ++i)
    {
        const Point3f &first = output[i];
        const Point3f &second = output[(i + 1) % size];
        twiceArea += first.cx*second.cy - second.cx*first.cy;
    }
    const double expectedArea = 4*scale*scale;
    assert(std::fabs(0.5*std::fabs(twiceArea) - expectedArea)
           <= expectedArea*1e-12);

    const double allForward[] = {scale, scale, scale, scale};
    assert(clip_polygon_by_nonnegative_scalar(
               polygon, allForward, 4,
               geometry_length_tolerance(4*scale), output) == 4);
    const double allBehind[] = {-scale, -scale, -scale, -scale};
    assert(clip_polygon_by_nonnegative_scalar(
               polygon, allBehind, 4,
               geometry_length_tolerance(4*scale), output) == 0);

    const double nearPlane[] = {
        -0.5*geometry_depth_tolerance(4*scale), scale, scale, scale
    };
    assert(clip_polygon_by_nonnegative_scalar(
               polygon, nearPlane, 4,
               geometry_depth_tolerance(4*scale), output) == 4);
}

void CheckPhasePlaneDistance()
{
    const Point3f direction(0.6, 0.0, -0.8);
    const Point3f phaseReference(2.0, -7.0, 3.0);
    const double front = -(direction.cx*phaseReference.cx
                         + direction.cy*phaseReference.cy
                         + direction.cz*phaseReference.cz);
    const Point3f target = phaseReference + direction*4.25;
    assert(std::fabs(SignedPhasePlaneDistance(direction, front, target) - 4.25)
           <= 1e-14);
    const Point3f behind = phaseReference - direction*1.75;
    assert(std::fabs(SignedPhasePlaneDistance(direction, front, behind) + 1.75)
           <= 1e-14);
}
}

int main()
{
    CheckScale(1.0);
    CheckScale(1.0e-6);
    CheckScale(1.0e6);
    CheckPhasePlaneDistance();
    return 0;
}
