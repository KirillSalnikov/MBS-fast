#include "particle/Bullet.h"
#include "particle/ConcaveHexagonal.h"
#include "particle/Hexagonal.h"
#include "particle/HexagonalAggregate.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

namespace
{
void WriteText(const char *path, const char *text)
{
    std::ofstream output(path);
    output << text;
    assert(output.good());
}

bool LoadFails(const char *path, const char *expected)
{
    try
    {
        std::unique_ptr<Particle> particle(new Particle());
        particle->SetFromFile(path);
    }
    catch (const std::runtime_error &error)
    {
        return std::string(error.what()).find(expected) != std::string::npos;
    }
    return false;
}
}

int main()
{
    const complex refractiveIndex(1.3116, 0.0);

    Polygon edited;
    edited.AddVertex(Point3f(0, 0, 0));
    edited.AddVertex(Point3f(1, 0, 0));
    edited.AddVertex(Point3f(2, 0, 0));
    edited.AddVertex(Point3f(3, 0, 0));
    edited.DeleteVertex(1);
    assert(edited.nVertices == 3);
    assert(edited.arr[0].cx == 0 && edited.arr[1].cx == 2
           && edited.arr[2].cx == 3);
    edited.InsertVertex(1, Point3f(1, 0, 0));
    assert(edited.nVertices == 4);
    for (int i = 0; i < edited.nVertices; ++i)
        assert(edited.arr[i].cx == i);

    PolygonArray emptyPolygons;
    bool emptyPopRejected = false;
    try
    {
        emptyPolygons.Pop();
    }
    catch (const std::out_of_range &)
    {
        emptyPopRejected = true;
    }
    assert(emptyPopRejected);

    std::unique_ptr<Hexagonal> column(
        new Hexagonal(refractiveIndex, 10.0, 7.0));
    assert(std::fabs(column->MaximumEdgeLength() - 7.0) < 1e-12);
    assert(column->MaximalDimention() > column->MaximumEdgeLength());
    column->Resize(24.0);
    const double expectedEdge = 7.0 * 24.0 / std::sqrt(149.0);
    assert(std::fabs(column->MaximumEdgeLength() - expectedEdge) < 1e-5);

    std::unique_ptr<ConcaveHexagonal> concave(
        new ConcaveHexagonal(refractiveIndex, 35.0, 50.0, 55.0079798014));

    for (int i = 0; i < 6; ++i)
    {
        assert(concave->facets[i].isVisibleIn);
        assert(!concave->facets[i].isVisibleOut);
    }

    for (int i = 6; i < 12; ++i)
    {
        assert(!concave->facets[i].isVisibleIn);
        assert(concave->facets[i].isVisibleOut);
    }

    for (int i = 12; i < 18; ++i)
    {
        assert(concave->facets[i].isVisibleIn);
        assert(!concave->facets[i].isVisibleOut);
    }

    const char *roundTripPath = "/tmp/mbs_fast_concave_visibility.dat";
    concave->Output(roundTripPath);
    std::unique_ptr<Particle> loaded(new Particle());
    loaded->SetFromFile(roundTripPath);
    assert(loaded->nFacets == concave->nFacets);
    for (int i = 0; i < concave->nFacets; ++i)
    {
        assert(loaded->facets[i].isVisibleIn == concave->facets[i].isVisibleIn);
        assert(loaded->facets[i].isVisibleOut == concave->facets[i].isVisibleOut);
    }
    std::remove(roundTripPath);
    std::remove("/tmp/mbs_fast_concave_visibility.dat.visibility");

    const char *validTetrahedron =
        "1\n0\n90 60\n"
        "0 0 0\n0 1 0\n1 0 0\n\n"
        "0 0 0\n1 0 0\n0 0 1\n\n"
        "0 0 0\n0 0 1\n0 1 0\n\n"
        "1 0 0\n0 1 0\n0 0 1\n";
    const char *validPath = "/tmp/mbs_fast_valid_tetrahedron.dat";
    WriteText(validPath, validTetrahedron);
    std::unique_ptr<Particle> valid(new Particle());
    valid->SetFromFile(validPath);
    assert(valid->nFacets == 4);
    std::remove(validPath);

    const char *openPath = "/tmp/mbs_fast_open_surface.dat";
    WriteText(openPath,
        "1\n0\n90 60\n"
        "0 0 0\n0 1 0\n1 0 0\n\n"
        "0 0 0\n1 0 0\n0 0 1\n\n"
        "0 0 0\n0 0 1\n0 1 0\n");
    assert(LoadFails(openPath, "open or non-manifold"));
    std::remove(openPath);

    // A coordinate mismatch smaller than one float ulp must not be collapsed
    // by the manifold vertex map now that geometry is stored in double.
    const char *subFloatGapPath = "/tmp/mbs_fast_sub_float_gap.dat";
    WriteText(subFloatGapPath,
        "1\n0\n90 60\n"
        "0 0 0\n0 1 0\n1 0 0\n\n"
        "0 0 0\n1.000000001 0 0\n0 0 1\n\n"
        "0 0 0\n0 0 1\n0 1 0\n\n"
        "1 0 0\n0 1 0\n0 0 1\n");
    assert(LoadFails(subFloatGapPath, "open or non-manifold"));
    std::remove(subFloatGapPath);

    const char *badEdgePath = "/tmp/mbs_fast_bad_edge_winding.dat";
    WriteText(badEdgePath,
        "1\n0\n90 60\n"
        "0 0 0\n0 1 0\n1 0 0\n\n"
        "0 0 0\n1 0 0\n0 0 1\n\n"
        "0 0 0\n0 1 0\n0 0 1\n\n"
        "1 0 0\n0 1 0\n0 0 1\n");
    assert(LoadFails(badEdgePath, "same direction"));
    std::remove(badEdgePath);

    const char *inwardPath = "/tmp/mbs_fast_inward_surface.dat";
    WriteText(inwardPath,
        "1\n0\n90 60\n"
        "1 0 0\n0 1 0\n0 0 0\n\n"
        "0 0 1\n1 0 0\n0 0 0\n\n"
        "0 1 0\n0 0 1\n0 0 0\n\n"
        "0 0 1\n0 1 0\n1 0 0\n");
    assert(LoadFails(inwardPath, "wound inward"));
    std::remove(inwardPath);

    const char *intersectingComponentsPath =
        "/tmp/mbs_fast_intersecting_components.dat";
    WriteText(intersectingComponentsPath,
        "1\n1 4\n90 60\n"
        "0 0 0\n0 1 0\n1 0 0\n\n"
        "0 0 0\n1 0 0\n0 0 1\n\n"
        "0 0 0\n0 0 1\n0 1 0\n\n"
        "1 0 0\n0 1 0\n0 0 1\n\n"
        "0.2 0.2 0.2\n0.2 1.2 0.2\n1.2 0.2 0.2\n\n"
        "0.2 0.2 0.2\n1.2 0.2 0.2\n0.2 0.2 1.2\n\n"
        "0.2 0.2 0.2\n0.2 0.2 1.2\n0.2 1.2 0.2\n\n"
        "1.2 0.2 0.2\n0.2 1.2 0.2\n0.2 0.2 1.2\n");
    assert(LoadFails(intersectingComponentsPath,
                     "non-adjacent surface facets intersect"));
    std::remove(intersectingComponentsPath);

    // Two closed pyramids touch only through crossing coplanar base patches.
    // No vertex of either base lies inside the other, so this specifically
    // exercises the positive-area 2-D overlap test rather than edge piercing.
    const char *coplanarOverlapPath =
        "/tmp/mbs_fast_coplanar_overlap.dat";
    WriteText(coplanarOverlapPath,
        "1\n1 5\n90 60\n"
        "-2 -0.5 0\n-2 0.5 0\n2 0.5 0\n2 -0.5 0\n\n"
        "-2 0.5 0\n-2 -0.5 0\n0 0 1\n\n"
        "2 0.5 0\n-2 0.5 0\n0 0 1\n\n"
        "2 -0.5 0\n2 0.5 0\n0 0 1\n\n"
        "-2 -0.5 0\n2 -0.5 0\n0 0 1\n\n"
        "-0.5 -2 0\n0.5 -2 0\n0.5 2 0\n-0.5 2 0\n\n"
        "0.5 -2 0\n-0.5 -2 0\n0 0 -1\n\n"
        "0.5 2 0\n0.5 -2 0\n0 0 -1\n\n"
        "-0.5 2 0\n0.5 2 0\n0 0 -1\n\n"
        "-0.5 -2 0\n-0.5 2 0\n0 0 -1\n");
    assert(LoadFails(coplanarOverlapPath, "coplanar surface facets overlap"));
    std::remove(coplanarOverlapPath);

    const char *selfIntersectingPath =
        "/tmp/mbs_fast_self_intersecting_facet.dat";
    WriteText(selfIntersectingPath,
        "1\n0\n90 60\n"
        "0 1 0\n0.587785 -0.809017 0\n-0.951057 0.309017 0\n"
        "0.951057 0.309017 0\n-0.587785 -0.809017 0\n");
    assert(LoadFails(selfIntersectingPath, "self-intersect"));
    std::remove(selfIntersectingPath);

    std::unique_ptr<Bullet> bullet(
        new Bullet(refractiveIndex, 35.0, 50.0, 20.0));
    for (int i = 0; i < bullet->nFacets; ++i)
    {
        assert(!bullet->facets[i].isVisibleIn);
        assert(!bullet->facets[i].isVisibleOut);
    }

    std::unique_ptr<HexagonalAggregate> aggregate(
        new HexagonalAggregate(refractiveIndex, 35.0, 50.0, 2));
    const int hiddenOut[] = {1, 2, 3, 7, 10, 11, 12, 15};
    for (int i = 0; i < aggregate->nFacets; ++i)
    {
        bool expectedOut = true;
        for (unsigned j = 0; j < sizeof(hiddenOut) / sizeof(hiddenOut[0]); ++j)
            if (hiddenOut[j] == i)
                expectedOut = false;
        assert(!aggregate->facets[i].isVisibleIn);
        assert(aggregate->facets[i].isVisibleOut == expectedOut);
    }

    return 0;
}
