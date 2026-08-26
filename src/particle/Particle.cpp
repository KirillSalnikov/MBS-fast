#include "Particle.h"
#include <algorithm>
#include <array>
#include <cfloat>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>
#include "global.h"
#include "Intersection.h"

using namespace std;
using ::complex;

namespace
{

struct ParticleFileLine
{
    int number;
    std::string text;
};

[[noreturn]] void ParticleFileError(const std::string &path, int line,
                                    const std::string &problem,
                                    const std::string &fix)
{
    std::ostringstream message;
    message << "invalid particle file '" << path << "'";
    if (line > 0)
        message << " at line " << line;
    message << ": " << problem << "\n  Fix: " << fix;
    throw std::runtime_error(message.str());
}

std::string Trim(const std::string &value)
{
    const std::string whitespace = " \t\r\n";
    const size_t first = value.find_first_not_of(whitespace);
    if (first == std::string::npos)
        return std::string();
    const size_t last = value.find_last_not_of(whitespace);
    return value.substr(first, last - first + 1);
}

std::vector<std::string> DataTokens(const std::string &text)
{
    std::vector<std::string> tokens;
    std::istringstream input(text);
    std::string token;
    while (input >> token)
    {
        if (!token.empty() && token[0] == '#')
            break;
        tokens.push_back(token);
    }
    return tokens;
}

long ParseIntegerToken(const std::string &token, const std::string &path,
                       int line, const std::string &name)
{
    char *end = nullptr;
    errno = 0;
    const long value = std::strtol(token.c_str(), &end, 10);
    if (errno == ERANGE || end == token.c_str() || *end != '\0')
    {
        ParticleFileError(path, line,
                          name + " must be an integer; got '" + token + "'.",
                          "replace it with an integer in the documented particle-file header.");
    }
    return value;
}

double ParseDoubleToken(const std::string &token, const std::string &path,
                        int line, const std::string &name)
{
    char *end = nullptr;
    errno = 0;
    const double value = std::strtod(token.c_str(), &end);
    if (errno == ERANGE || end == token.c_str() || *end != '\0'
        || !std::isfinite(value))
    {
        ParticleFileError(path, line,
                          name + " must be a finite number; got '" + token + "'.",
                          "replace it with a finite decimal value.");
    }
    return value;
}

ParticleFileLine NextHeaderLine(std::ifstream &input, const std::string &path,
                                int &lineNumber, const std::string &name)
{
    std::string text;
    while (std::getline(input, text))
    {
        ++lineNumber;
        const std::string trimmed = Trim(text);
        if (trimmed.empty() || trimmed[0] == '#')
            continue;
        return ParticleFileLine{lineNumber, trimmed};
    }
    ParticleFileError(path, lineNumber,
                      "missing " + name + " header.",
                      "provide the three required header records before the facet vertices.");
}

void ValidateFacetPolygon(const std::vector<Point3f> &vertices,
                          const std::string &path, int line)
{
    Facet facet;
    for (const Point3f &point : vertices)
        facet.AddVertex(point);

    const Point3f normal = facet.Normal();
    double diameter = 0.0;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        for (size_t j = i + 1; j < vertices.size(); ++j)
            diameter = std::max(diameter, Length(vertices[j] - vertices[i]));
    }
    const double lengthTolerance = geometry_length_tolerance(diameter);
    const double areaTolerance = geometry_area_tolerance(diameter);
    if (Length(normal) <= EPS_PROJECTION)
    {
        ParticleFileError(
            path, line,
            "a facet has no stable normal, usually because its boundary self-intersects.",
            "write one simple convex boundary or split the facet into convex polygons.");
    }

    for (size_t i = 0; i < vertices.size(); ++i)
    {
        for (size_t j = i + 1; j < vertices.size(); ++j)
        {
            if (Length(vertices[j] - vertices[i]) <= lengthTolerance)
            {
                ParticleFileError(
                    path, line,
                    "a facet contains duplicate or numerically indistinguishable vertices.",
                    "remove repeated vertices and keep every edge longer than the geometry tolerance.");
            }
        }

        const double planeDistance = std::fabs(
            DotProduct(vertices[i] - vertices[0], normal));
        if (planeDistance > lengthTolerance)
        {
            ParticleFileError(
                path, line,
                "a facet is not planar within the geometry tolerance.",
                "triangulate the facet or project all of its vertices onto one plane.");
        }

        const Point3f edgeA = vertices[(i + 1) % vertices.size()] - vertices[i];
        const Point3f edgeB = vertices[(i + 2) % vertices.size()]
                            - vertices[(i + 1) % vertices.size()];
        const double turn = DotProduct(CrossProduct(edgeA, edgeB), normal);
        if (turn < -areaTolerance)
        {
            ParticleFileError(
                path, line,
                "a facet is non-convex or its vertex order self-intersects.",
                "split it into convex planar facets with one consistent winding.");
        }
    }

    const double ax = std::fabs(static_cast<double>(normal.cx));
    const double ay = std::fabs(static_cast<double>(normal.cy));
    const double az = std::fabs(static_cast<double>(normal.cz));
    const int drop = ax > ay && ax > az ? 0 : (ay > az ? 1 : 2);
    const auto project = [drop](const Point3f &point) {
        if (drop == 0)
            return std::make_pair(static_cast<double>(point.cy),
                                  static_cast<double>(point.cz));
        if (drop == 1)
            return std::make_pair(static_cast<double>(point.cx),
                                  static_cast<double>(point.cz));
        return std::make_pair(static_cast<double>(point.cx),
                              static_cast<double>(point.cy));
    };
    const auto orient2d = [](const std::pair<double, double> &a,
                             const std::pair<double, double> &b,
                             const std::pair<double, double> &c) {
        return (b.first - a.first)*(c.second - a.second)
             - (b.second - a.second)*(c.first - a.first);
    };
    const auto rangesOverlap = [lengthTolerance](double a0, double a1,
                                                  double b0, double b1) {
        if (a0 > a1) std::swap(a0, a1);
        if (b0 > b1) std::swap(b0, b1);
        return a1 + lengthTolerance >= b0
            && b1 + lengthTolerance >= a0;
    };
    const auto side = [areaTolerance](double value) {
        return value > areaTolerance ? 1 : (value < -areaTolerance ? -1 : 0);
    };

    for (size_t first = 0; first < vertices.size(); ++first)
    {
        const size_t firstEnd = (first + 1) % vertices.size();
        const std::pair<double, double> a = project(vertices[first]);
        const std::pair<double, double> b = project(vertices[firstEnd]);
        for (size_t second = first + 1; second < vertices.size(); ++second)
        {
            const size_t secondEnd = (second + 1) % vertices.size();
            if (second == firstEnd || secondEnd == first)
                continue;

            const std::pair<double, double> c = project(vertices[second]);
            const std::pair<double, double> d = project(vertices[secondEnd]);
            if (!rangesOverlap(a.first, b.first, c.first, d.first)
                || !rangesOverlap(a.second, b.second, c.second, d.second))
                continue;

            const int s1 = side(orient2d(a, b, c));
            const int s2 = side(orient2d(a, b, d));
            const int s3 = side(orient2d(c, d, a));
            const int s4 = side(orient2d(c, d, b));
            if (s1*s2 <= 0 && s3*s4 <= 0)
            {
                ParticleFileError(
                    path, line,
                    "a facet boundary self-intersects or contains overlapping non-adjacent edges.",
                    "write a simple convex boundary and split complex faces into separate facets.");
            }
        }
    }
}

struct ProjectedPoint
{
    double x;
    double y;
};

double Cross2d(const ProjectedPoint &a, const ProjectedPoint &b,
               const ProjectedPoint &c)
{
    return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

double SignedArea2d(const std::vector<ProjectedPoint> &polygon)
{
    double twiceArea = 0.0;
    for (size_t i = 0; i < polygon.size(); ++i)
    {
        const ProjectedPoint &a = polygon[i];
        const ProjectedPoint &b = polygon[(i + 1) % polygon.size()];
        twiceArea += a.x*b.y - a.y*b.x;
    }
    return 0.5*twiceArea;
}

std::vector<ProjectedPoint> ProjectFacet2d(
    const std::vector<Point3f> &facet, int drop)
{
    std::vector<ProjectedPoint> projected;
    projected.reserve(facet.size());
    for (const Point3f &point : facet)
    {
        if (drop == 0)
            projected.push_back({point.cy, point.cz});
        else if (drop == 1)
            projected.push_back({point.cx, point.cz});
        else
            projected.push_back({point.cx, point.cy});
    }
    return projected;
}

std::vector<ProjectedPoint> IntersectConvexPolygons2d(
    const std::vector<ProjectedPoint> &subject,
    const std::vector<ProjectedPoint> &clip,
    double areaTolerance)
{
    std::vector<ProjectedPoint> result = subject;
    const double orientation = SignedArea2d(clip) >= 0.0 ? 1.0 : -1.0;
    for (size_t edge = 0; edge < clip.size() && !result.empty(); ++edge)
    {
        const ProjectedPoint &a = clip[edge];
        const ProjectedPoint &b = clip[(edge + 1) % clip.size()];
        std::vector<ProjectedPoint> input;
        input.swap(result);
        ProjectedPoint previous = input.back();
        double previousSide = orientation*Cross2d(a, b, previous);
        bool previousInside = previousSide >= -areaTolerance;
        for (const ProjectedPoint &current : input)
        {
            const double currentSide = orientation*Cross2d(a, b, current);
            const bool currentInside = currentSide >= -areaTolerance;
            if (currentInside != previousInside)
            {
                const double denominator = previousSide - currentSide;
                if (std::fabs(denominator) > DBL_MIN)
                {
                    const double t = previousSide/denominator;
                    result.push_back({
                        previous.x + t*(current.x - previous.x),
                        previous.y + t*(current.y - previous.y)});
                }
            }
            if (currentInside)
                result.push_back(current);
            previous = current;
            previousSide = currentSide;
            previousInside = currentInside;
        }
    }
    return result;
}

bool HasCoplanarAreaOverlap(const std::vector<Point3f> &first,
                            const Point3f &firstNormal,
                            const std::vector<Point3f> &second,
                            const Point3f &secondNormal,
                            double lengthTolerance,
                            double areaTolerance)
{
    if (Length(CrossProduct(firstNormal, secondNormal))
        > GEOMETRY_RELATIVE_EPSILON)
        return false;
    if (std::fabs(DotProduct(second[0] - first[0], firstNormal))
        > lengthTolerance)
        return false;

    const double ax = std::fabs(firstNormal.cx);
    const double ay = std::fabs(firstNormal.cy);
    const double az = std::fabs(firstNormal.cz);
    const int drop = ax >= ay && ax >= az ? 0 : (ay >= az ? 1 : 2);
    const std::vector<ProjectedPoint> intersection = IntersectConvexPolygons2d(
        ProjectFacet2d(first, drop), ProjectFacet2d(second, drop), areaTolerance);
    return intersection.size() >= 3
        && std::fabs(SignedArea2d(intersection)) > areaTolerance;
}

void ValidateClosedOrientedSurface(
    const std::vector<std::vector<Point3f> > &facets,
    const std::vector<int> &facetLines,
    const std::string &path)
{
    typedef std::tuple<double, double, double> VertexKey;
    struct EdgeUse
    {
        int count = 0;
        int direction = 0;
        int firstFacet = -1;
        int secondFacet = -1;
    };

    std::map<VertexKey, int> vertexIds;
    std::map<std::pair<int, int>, EdgeUse> edges;
    std::vector<std::vector<int> > facetVertexIds;
    facetVertexIds.reserve(facets.size());
    int nextVertexId = 0;
    for (size_t facetIndex = 0; facetIndex < facets.size(); ++facetIndex)
    {
        const std::vector<Point3f> &facet = facets[facetIndex];
        std::vector<int> ids;
        ids.reserve(facet.size());
        for (const Point3f &point : facet)
        {
            const VertexKey key(point.cx, point.cy, point.cz);
            std::map<VertexKey, int>::iterator found = vertexIds.find(key);
            if (found == vertexIds.end())
                found = vertexIds.insert(std::make_pair(key, nextVertexId++)).first;
            ids.push_back(found->second);
        }

        for (size_t i = 0; i < ids.size(); ++i)
        {
            const int from = ids[i];
            const int to = ids[(i + 1) % ids.size()];
            const std::pair<int, int> key = from < to
                ? std::make_pair(from, to) : std::make_pair(to, from);
            EdgeUse &use = edges[key];
            if (use.count == 0)
                use.firstFacet = static_cast<int>(facetIndex);
            else if (use.count == 1)
                use.secondFacet = static_cast<int>(facetIndex);
            ++use.count;
            use.direction += from < to ? 1 : -1;
        }
        facetVertexIds.push_back(ids);
    }

    std::vector<std::vector<int> > adjacency(facets.size());
    for (std::map<std::pair<int, int>, EdgeUse>::const_iterator it = edges.begin();
         it != edges.end(); ++it)
    {
        const EdgeUse &use = it->second;
        const int facet = std::max(0, use.firstFacet);
        if (use.count != 2)
        {
            ParticleFileError(
                path, facetLines[facet],
                "the surface is open or non-manifold: an edge is used "
                    + std::to_string(use.count) + " times instead of two.",
                "close every component and make each undirected edge belong to exactly two facets.");
        }
        if (use.direction != 0)
        {
            ParticleFileError(
                path, facetLines[facet],
                "adjacent facets traverse a shared edge in the same direction.",
                "reverse one facet so all outward normals use a consistent winding.");
        }
        adjacency[use.firstFacet].push_back(use.secondFacet);
        adjacency[use.secondFacet].push_back(use.firstFacet);
    }

    // Edge incidence alone does not reject two otherwise valid closed
    // components that pass through each other. Catch every non-coplanar
    // edge/facet penetration between facets that do not share a vertex.
    double surfaceMin[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double surfaceMax[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    std::vector<Point3f> facetNormals;
    facetNormals.reserve(facets.size());
    for (const std::vector<Point3f> &facet : facets)
    {
        Facet polygon;
        for (const Point3f &point : facet)
        {
            polygon.AddVertex(point);
            for (int coordinate = 0; coordinate < 3; ++coordinate)
            {
                surfaceMin[coordinate] = std::min(
                    surfaceMin[coordinate], point.coordinates[coordinate]);
                surfaceMax[coordinate] = std::max(
                    surfaceMax[coordinate], point.coordinates[coordinate]);
            }
        }
        facetNormals.push_back(polygon.Normal());
    }
    const double surfaceDx = surfaceMax[0] - surfaceMin[0];
    const double surfaceDy = surfaceMax[1] - surfaceMin[1];
    const double surfaceDz = surfaceMax[2] - surfaceMin[2];
    const double surfaceScale = std::sqrt(
        surfaceDx*surfaceDx + surfaceDy*surfaceDy + surfaceDz*surfaceDz);
    const double lengthTolerance = geometry_length_tolerance(surfaceScale);
    const double areaTolerance = geometry_area_tolerance(surfaceScale);

    const auto pointInsideFacet = [areaTolerance](
        const Point3f &point, const std::vector<Point3f> &facet,
        const Point3f &normal) {
        for (size_t edge = 0; edge < facet.size(); ++edge)
        {
            const Point3f &a = facet[edge];
            const Point3f &b = facet[(edge + 1) % facet.size()];
            if (DotProduct(CrossProduct(b - a, point - a), normal)
                < -areaTolerance)
                return false;
        }
        return true;
    };
    const auto edgePiercesFacet = [lengthTolerance, &pointInsideFacet](
        const Point3f &begin, const Point3f &end,
        const std::vector<Point3f> &facet, const Point3f &normal) {
        const double beginDistance = DotProduct(begin - facet[0], normal);
        const double endDistance = DotProduct(end - facet[0], normal);
        if ((beginDistance > lengthTolerance && endDistance > lengthTolerance)
            || (beginDistance < -lengthTolerance
                && endDistance < -lengthTolerance))
            return false;

        // Coplanar overlap requires a 2-D test and is diagnosed separately by
        // the facet/edge validators when boundaries coincide. Do not divide by
        // a near-zero plane-distance difference here.
        const double denominator = beginDistance - endDistance;
        if (std::fabs(denominator) <= lengthTolerance)
            return false;
        const double position = beginDistance / denominator;
        if (position < -GEOMETRY_RELATIVE_EPSILON
            || position > 1.0 + GEOMETRY_RELATIVE_EPSILON)
            return false;
        const Point3f point = begin + (end - begin)*position;
        return pointInsideFacet(point, facet, normal);
    };

    for (size_t first = 0; first < facets.size(); ++first)
    {
        for (size_t second = first + 1; second < facets.size(); ++second)
        {
            if (HasCoplanarAreaOverlap(
                    facets[first], facetNormals[first],
                    facets[second], facetNormals[second],
                    lengthTolerance, areaTolerance))
            {
                ParticleFileError(
                    path, facetLines[first],
                    "non-adjacent coplanar surface facets overlap with positive area.",
                    "remove duplicate or overlapping coplanar patches before tracing.");
            }

            bool sharesVertex = false;
            for (int firstVertex : facetVertexIds[first])
            {
                if (std::find(facetVertexIds[second].begin(),
                              facetVertexIds[second].end(), firstVertex)
                    != facetVertexIds[second].end())
                {
                    sharesVertex = true;
                    break;
                }
            }
            if (sharesVertex)
                continue;

            const auto pierces = [&facets, &facetNormals, &edgePiercesFacet](
                size_t edgeFacet, size_t targetFacet) {
                const std::vector<Point3f> &edges = facets[edgeFacet];
                for (size_t edge = 0; edge < edges.size(); ++edge)
                {
                    if (edgePiercesFacet(
                            edges[edge], edges[(edge + 1) % edges.size()],
                            facets[targetFacet], facetNormals[targetFacet]))
                        return true;
                }
                return false;
            };
            if (pierces(first, second) || pierces(second, first))
            {
                ParticleFileError(
                    path, facetLines[first],
                    "non-adjacent surface facets intersect.",
                    "remove overlapping components and facet penetrations before tracing.");
            }
        }
    }

    std::vector<unsigned char> visited(facets.size(), 0);
    for (size_t seed = 0; seed < facets.size(); ++seed)
    {
        if (visited[seed])
            continue;

        std::vector<int> component(1, static_cast<int>(seed));
        visited[seed] = 1;
        for (size_t cursor = 0; cursor < component.size(); ++cursor)
        {
            const int facet = component[cursor];
            for (int neighbor : adjacency[facet])
            {
                if (!visited[neighbor])
                {
                    visited[neighbor] = 1;
                    component.push_back(neighbor);
                }
            }
        }

        const Point3f origin = facets[component[0]][0];
        double minCoord[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
        double maxCoord[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
        double signedVolume = 0.0;
        for (int facetIndex : component)
        {
            const std::vector<Point3f> &facet = facets[facetIndex];
            for (const Point3f &point : facet)
            {
                for (int coordinate = 0; coordinate < 3; ++coordinate)
                {
                    minCoord[coordinate] = std::min(
                        minCoord[coordinate],
                        static_cast<double>(point.coordinates[coordinate]));
                    maxCoord[coordinate] = std::max(
                        maxCoord[coordinate],
                        static_cast<double>(point.coordinates[coordinate]));
                }
            }
            const Point3f a = facet[0] - origin;
            for (size_t i = 1; i + 1 < facet.size(); ++i)
            {
                const Point3f b = facet[i] - origin;
                const Point3f c = facet[i + 1] - origin;
                signedVolume += DotProduct(a, CrossProduct(b, c)) / 6.0;
            }
        }

        const double dx = maxCoord[0] - minCoord[0];
        const double dy = maxCoord[1] - minCoord[1];
        const double dz = maxCoord[2] - minCoord[2];
        const double scale = std::sqrt(dx*dx + dy*dy + dz*dz);
        const double volumeTolerance = GEOMETRY_RELATIVE_EPSILON
                                     * scale * scale * scale;
        if (std::fabs(signedVolume) <= volumeTolerance)
        {
            ParticleFileError(
                path, facetLines[component[0]],
                "a closed surface component has zero or numerically unstable enclosed volume.",
                "remove self-intersections and give the component a finite three-dimensional thickness.");
        }
        if (signedVolume < 0.0)
        {
            ParticleFileError(
                path, facetLines[component[0]],
                "a closed surface component is wound inward.",
                "reverse every facet of this component so its normals point outward.");
        }
    }
}

double MaximumPairDistance(std::vector<Point3f> &vertices)
{
    if (vertices.size() < 2)
        return 0.0;

    std::sort(vertices.begin(), vertices.end(),
        [](const Point3f &a, const Point3f &b) {
            if (a.cx != b.cx)
                return a.cx < b.cx;
            if (a.cy != b.cy)
                return a.cy < b.cy;
            return a.cz < b.cz;
        });
    vertices.erase(std::unique(vertices.begin(), vertices.end(),
        [](const Point3f &a, const Point3f &b) {
            return a.cx == b.cx && a.cy == b.cy && a.cz == b.cz;
        }), vertices.end());

    double maximumSquared = 0.0;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        for (size_t j = i + 1; j < vertices.size(); ++j)
        {
            const Point3f delta = vertices[j] - vertices[i];
            maximumSquared = std::max(maximumSquared, Norm(delta));
        }
    }
    return std::sqrt(maximumSquared);
}

} // namespace

Particle::Particle()
{
    isConcave = false;
}

void Particle::SetFromFile(const std::string &filename)
{
    std::ifstream pfile(filename, std::ios::in);

    if (!pfile.is_open())
    {
        ParticleFileError(filename, 0, "the file cannot be opened.",
                          "check that the path exists and is readable.");
    }

    int lineNumber = 0;
    const ParticleFileLine concavityLine = NextHeaderLine(
        pfile, filename, lineNumber, "concavity");
    const std::vector<std::string> concavity = DataTokens(concavityLine.text);
    if (concavity.size() != 1)
        ParticleFileError(filename, concavityLine.number,
                          "the concavity header expects one value (0 or 1).",
                          "write '0' for convex geometry or '1' for nonconvex geometry.");
    const long concavityValue = ParseIntegerToken(
        concavity[0], filename, concavityLine.number, "concavity flag");
    if (concavityValue != 0 && concavityValue != 1)
        ParticleFileError(filename, concavityLine.number,
                          "the concavity flag must be 0 or 1.",
                          "write '0' for convex geometry or '1' for nonconvex geometry.");

    const ParticleFileLine aggregateLine = NextHeaderLine(
        pfile, filename, lineNumber, "aggregate");
    const std::vector<std::string> aggregate = DataTokens(aggregateLine.text);
    if (aggregate.empty() || aggregate.size() > 2)
        ParticleFileError(filename, aggregateLine.number,
                          "the aggregate header expects '0' or '1 FACETS_PER_PART'.",
                          "write '0' for one particle, or for example '1 8' for an aggregate.");
    const long aggregateValue = ParseIntegerToken(
        aggregate[0], filename, aggregateLine.number, "aggregate flag");
    if (aggregateValue != 0 && aggregateValue != 1)
        ParticleFileError(filename, aggregateLine.number,
                          "the aggregate flag must be 0 or 1.",
                          "write '0' for one particle or '1 FACETS_PER_PART' for an aggregate.");
    if (aggregateValue == 0 && aggregate.size() != 1)
        ParticleFileError(filename, aggregateLine.number,
                          "FACETS_PER_PART is present while the aggregate flag is 0.",
                          "remove the second value or change the header to '1 FACETS_PER_PART'.");
    if (aggregateValue == 1 && aggregate.size() != 2)
        ParticleFileError(filename, aggregateLine.number,
                          "an aggregate requires FACETS_PER_PART.",
                          "write the number of facets in each component, for example '1 8'.");

    int parsedFacetsPerPart = 0;
    if (aggregateValue == 1)
    {
        const long value = ParseIntegerToken(
            aggregate[1], filename, aggregateLine.number, "FACETS_PER_PART");
        if (value < 1 || value > MAX_FACET_NUM)
            ParticleFileError(filename, aggregateLine.number,
                              "FACETS_PER_PART is outside the supported range.",
                              "use an integer from 1 to " + std::to_string(MAX_FACET_NUM) + ".");
        parsedFacetsPerPart = static_cast<int>(value);
    }

    const ParticleFileLine symmetryLine = NextHeaderLine(
        pfile, filename, lineNumber, "symmetry");
    const std::vector<std::string> symmetry = DataTokens(symmetryLine.text);
    if (symmetry.size() != 2)
        ParticleFileError(filename, symmetryLine.number,
                          "the symmetry header expects BETA_DEG GAMMA_DEG.",
                          "provide two positive angular-domain widths, for example '90 60'.");
    const double betaDeg = ParseDoubleToken(
        symmetry[0], filename, symmetryLine.number, "beta symmetry");
    const double gammaDeg = ParseDoubleToken(
        symmetry[1], filename, symmetryLine.number, "gamma symmetry");
    if (!(betaDeg > 0.0 && betaDeg <= 180.0))
        ParticleFileError(filename, symmetryLine.number,
                          "beta symmetry must be in (0, 180] degrees.",
                          "use the particle's positive beta fundamental-domain width.");
    if (!(gammaDeg > 0.0 && gammaDeg <= 360.0))
        ParticleFileError(filename, symmetryLine.number,
                          "gamma symmetry must be in (0, 360] degrees.",
                          "use the particle's positive gamma fundamental-domain width.");

    typedef std::array<double, 3> ParsedCoordinate;
    std::vector<std::vector<ParsedCoordinate> > parsedCoordinates;
    std::vector<std::vector<Point3f> > parsedFacets;
    std::vector<int> parsedFacetLines;
    std::vector<ParsedCoordinate> vertices;
    int facetStartLine = 0;
    const auto finishFacet = [&](int separatorLine) {
        if (vertices.empty())
            return;
        if (vertices.size() < 3)
            ParticleFileError(filename, facetStartLine,
                              "a facet has fewer than three vertices.",
                              "add vertices or remove the incomplete facet before line "
                                  + std::to_string(separatorLine) + ".");
        if (parsedCoordinates.size() >= MAX_FACET_NUM)
            ParticleFileError(filename, facetStartLine,
                              "the particle exceeds the facet limit.",
                              "reduce the geometry to at most "
                                  + std::to_string(MAX_FACET_NUM) + " facets.");

        parsedCoordinates.push_back(vertices);
        parsedFacetLines.push_back(facetStartLine);
        vertices.clear();
        facetStartLine = 0;
    };

    std::string text;
    while (std::getline(pfile, text))
    {
        ++lineNumber;
        const std::string trimmed = Trim(text);
        if (trimmed.empty())
        {
            finishFacet(lineNumber);
            continue;
        }
        if (trimmed[0] == '#')
            continue;

        const std::vector<std::string> coordinates = DataTokens(trimmed);
        if (coordinates.size() != 3)
            ParticleFileError(filename, lineNumber,
                              "a vertex expects exactly three coordinates X Y Z.",
                              "write one finite three-dimensional vertex per line.");
        if (vertices.size() >= MAX_VERTEX_NUM)
            ParticleFileError(filename, lineNumber,
                              "a facet exceeds the vertex limit.",
                              "split or simplify it to at most "
                                  + std::to_string(MAX_VERTEX_NUM) + " vertices.");
        if (vertices.empty())
            facetStartLine = lineNumber;
        const ParsedCoordinate coordinate = {{
            ParseDoubleToken(coordinates[0], filename, lineNumber, "X coordinate"),
            ParseDoubleToken(coordinates[1], filename, lineNumber, "Y coordinate"),
            ParseDoubleToken(coordinates[2], filename, lineNumber, "Z coordinate")}};
        vertices.push_back(coordinate);
    }
    if (pfile.bad())
        ParticleFileError(filename, lineNumber,
                          "an I/O error interrupted the particle body.",
                          "verify the file storage and permissions, then read the complete file again.");
    finishFacet(lineNumber + 1);

    if (parsedCoordinates.empty())
        ParticleFileError(filename, lineNumber,
                          "the file contains no facets.",
                          "add facet polygons after the three header records; separate facets with blank lines.");

    // Preserve the source coordinates in double until an order-independent
    // local origin is known. Choosing the first file vertex as origin changes
    // edge rounding when
    // facets or cyclic vertex starts are serialized in another order.
    double minCoordinate[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double maxCoordinate[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    for (const std::vector<ParsedCoordinate> &facet : parsedCoordinates)
    {
        for (const ParsedCoordinate &point : facet)
        {
            for (int coordinate = 0; coordinate < 3; ++coordinate)
            {
                minCoordinate[coordinate] = std::min(
                    minCoordinate[coordinate], point[coordinate]);
                maxCoordinate[coordinate] = std::max(
                    maxCoordinate[coordinate], point[coordinate]);
            }
        }
    }
    const double canonicalOrigin[3] = {
        0.5*(minCoordinate[0] + maxCoordinate[0]),
        0.5*(minCoordinate[1] + maxCoordinate[1]),
        0.5*(minCoordinate[2] + maxCoordinate[2])};
    const double spanX = maxCoordinate[0] - minCoordinate[0];
    const double spanY = maxCoordinate[1] - minCoordinate[1];
    const double spanZ = maxCoordinate[2] - minCoordinate[2];
    const double canonicalScale = std::sqrt(
        spanX*spanX + spanY*spanY + spanZ*spanZ);
    if (!(canonicalScale > 0.0) || !std::isfinite(canonicalScale))
        ParticleFileError(filename, 0, "the particle extent is zero or non-finite.",
                          "provide a finite three-dimensional closed surface.");

    // Text coordinates translated far from the origin lose a few low bits
    // before they reach this parser. Snap only those nonphysical serialization
    // bits to a power-of-two grid, so equivalent translated/rewritten meshes
    // produce identical local geometry and deterministic facet ordering.
    const int scaleExponent = std::ilogb(canonicalScale);
    const double canonicalQuantum = std::ldexp(1.0, scaleExponent - 40);
    const auto canonicalCoordinate = [canonicalQuantum](double value) {
        return std::nearbyint(value/canonicalQuantum)*canonicalQuantum;
    };
    parsedFacets.reserve(parsedCoordinates.size());
    for (size_t facetIndex = 0; facetIndex < parsedCoordinates.size(); ++facetIndex)
    {
        std::vector<Point3f> facet;
        facet.reserve(parsedCoordinates[facetIndex].size());
        for (const ParsedCoordinate &point : parsedCoordinates[facetIndex])
        {
            facet.push_back(Point3f(
                canonicalCoordinate(point[0] - canonicalOrigin[0]),
                canonicalCoordinate(point[1] - canonicalOrigin[1]),
                canonicalCoordinate(point[2] - canonicalOrigin[2])));
        }
        Facet candidate;
        for (const Point3f &point : facet)
            candidate.AddVertex(point);
        const double area = candidate.Area();
        if (!(area > 0.0) || !std::isfinite(area))
            ParticleFileError(filename, parsedFacetLines[facetIndex],
                              "a facet has zero or non-finite area.",
                              "remove duplicate/collinear vertices and keep a nondegenerate polygon.");
        ValidateFacetPolygon(facet, filename, parsedFacetLines[facetIndex]);
        parsedFacets.push_back(facet);
    }

    if (aggregateValue == 1
        && parsedFacets.size() % static_cast<size_t>(parsedFacetsPerPart) != 0)
    {
        ParticleFileError(filename, aggregateLine.number,
                          "the total facet count is not divisible by FACETS_PER_PART.",
                          "correct FACETS_PER_PART or complete every aggregate component.");
    }
    ValidateClosedOrientedSurface(parsedFacets, parsedFacetLines, filename);

    // Facet ids are not part of the particle-file format, yet near-coplanar
    // beam subtraction uses them as a deterministic final tie-break. Assign
    // ids from geometry so rewriting the same surface in another facet order
    // cannot change the beam tree. Aggregate files retain their component
    // blocks because FACETS_PER_PART makes that ordering semantic.
    std::vector<size_t> facetOrder(parsedFacets.size());
    for (size_t i = 0; i < facetOrder.size(); ++i)
        facetOrder[i] = i;
    const auto pointLess = [](const Point3f &a, const Point3f &b) {
        if (a.cz != b.cz) return a.cz > b.cz;
        if (a.cy != b.cy) return a.cy < b.cy;
        return a.cx < b.cx;
    };
    for (std::vector<Point3f> &facet : parsedFacets)
    {
        const std::vector<Point3f>::iterator first = std::min_element(
            facet.begin(), facet.end(), pointLess);
        std::rotate(facet.begin(), first, facet.end());
    }
    if (aggregateValue == 0)
    {
        const auto center = [](const std::vector<Point3f> &facet) {
            double x = 0.0, y = 0.0, z = 0.0;
            for (const Point3f &point : facet)
            {
                x += point.cx;
                y += point.cy;
                z += point.cz;
            }
            const double inverse = 1.0/facet.size();
            return Point3f(x*inverse, y*inverse, z*inverse);
        };
        std::sort(facetOrder.begin(), facetOrder.end(),
            [&parsedFacets, &pointLess, &center](size_t leftIndex,
                                                size_t rightIndex) {
                const std::vector<Point3f> &left = parsedFacets[leftIndex];
                const std::vector<Point3f> &right = parsedFacets[rightIndex];
                const Point3f leftCenter = center(left);
                const Point3f rightCenter = center(right);
                if (pointLess(leftCenter, rightCenter)) return true;
                if (pointLess(rightCenter, leftCenter)) return false;
                if (left.size() != right.size()) return left.size() < right.size();

                size_t leftFirst = 0;
                size_t rightFirst = 0;
                for (size_t i = 1; i < left.size(); ++i)
                    if (pointLess(left[i], left[leftFirst])) leftFirst = i;
                for (size_t i = 1; i < right.size(); ++i)
                    if (pointLess(right[i], right[rightFirst])) rightFirst = i;
                for (size_t i = 0; i < left.size(); ++i)
                {
                    const Point3f &leftPoint =
                        left[(leftFirst + i) % left.size()];
                    const Point3f &rightPoint =
                        right[(rightFirst + i) % right.size()];
                    if (pointLess(leftPoint, rightPoint)) return true;
                    if (pointLess(rightPoint, leftPoint)) return false;
                }
                return leftIndex < rightIndex;
            });

        std::vector<std::vector<Point3f> > orderedFacets;
        std::vector<int> orderedLines;
        orderedFacets.reserve(parsedFacets.size());
        orderedLines.reserve(parsedFacetLines.size());
        for (size_t originalIndex : facetOrder)
        {
            orderedFacets.push_back(parsedFacets[originalIndex]);
            orderedLines.push_back(parsedFacetLines[originalIndex]);
        }
        parsedFacets.swap(orderedFacets);
        parsedFacetLines.swap(orderedLines);
    }

    nFacets = static_cast<int>(parsedFacets.size());
    isConcave = concavityValue != 0;
    isAggregated = aggregateValue != 0;
    nFacetsInPart = parsedFacetsPerPart;
    m_symmetry.beta = DegToRad(betaDeg);
    m_symmetry.gamma = DegToRad(gammaDeg);
    m_symmetry.alpha = 0.0;

    for (int i = 0; i < nFacets; ++i)
    {
        defaultFacets[i].Clear();
        defaultFacets[i].isVisibleIn = true;
        defaultFacets[i].isVisibleOut = true;
        facets[i].Clear();
        facets[i].isVisibleIn = true;
        facets[i].isVisibleOut = true;
        for (const Point3f &point : parsedFacets[i])
            defaultFacets[i].AddVertex(point);
    }

    if (isConcave)
    {
//        RemoveWalls();
    }

    SetDefaultNormals();
    Reset();
    SetDParams();
    SetDefaultCenters();

    if (isConcave || isAggregated)
    {
        for (int i = 0; i < nFacets; ++i)
        {
            defaultFacets[i].isVisibleIn = false;
            defaultFacets[i].isVisibleOut = false;
            facets[i].isVisibleIn = false;
            facets[i].isVisibleOut = false;
        }
    }

    const std::string visibilityPath = filename + ".visibility";
    std::ifstream visibilityFile(visibilityPath.c_str(), std::ios::in);
    if (visibilityFile.is_open())
    {
        std::vector<std::pair<bool, bool> > visibilityRecords;
        visibilityRecords.reserve(nFacets);
        int storedFacetCount = 0;
        if (!(visibilityFile >> storedFacetCount) || storedFacetCount != nFacets)
        {
            ParticleFileError(
                visibilityPath, 1,
                "the visibility facet count does not match the geometry.",
                "regenerate the .visibility sidecar together with its particle file.");
        }
        for (int i = 0; i < nFacets; ++i)
        {
            int visibleIn = -1;
            int visibleOut = -1;
            if (!(visibilityFile >> visibleIn >> visibleOut)
                || (visibleIn != 0 && visibleIn != 1)
                || (visibleOut != 0 && visibleOut != 1))
            {
                ParticleFileError(
                    visibilityPath, i + 2,
                    "a visibility record must contain two values (0 or 1).",
                    "write one 'VISIBLE_IN VISIBLE_OUT' record per facet.");
            }
            visibilityRecords.push_back(
                std::make_pair(visibleIn != 0, visibleOut != 0));
        }
        std::string extra;
        if (visibilityFile >> extra)
        {
            ParticleFileError(
                visibilityPath, nFacets + 2,
                "the visibility sidecar contains extra records.",
                "keep exactly one visibility record per facet.");
        }
        for (int i = 0; i < nFacets; ++i)
        {
            const size_t originalIndex = facetOrder[i];
            defaultFacets[i].isVisibleIn = visibilityRecords[originalIndex].first;
            defaultFacets[i].isVisibleOut = visibilityRecords[originalIndex].second;
            facets[i].isVisibleIn = visibilityRecords[originalIndex].first;
            facets[i].isVisibleOut = visibilityRecords[originalIndex].second;
        }
    }
}

void Particle::Init(int facetCount, const ::complex &refrIndex)
{
    nFacets = facetCount;
    m_refractiveIndex = refrIndex;
}

void Particle::RotateCenters()
{
    for (int i = 0; i < nFacets; ++i)
    {
        RotatePoint(defaultFacets[i].center, facets[i].center);
    }
}

void Particle::Rotate(double beta, double gamma, double alpha)
{
    rotAngle = Angle{alpha, beta, gamma};
    SetRotateMatrix(beta, gamma, alpha);
    ApplyRotateMatrix();
}

void Particle::RotateQuaternion(double qx, double qy, double qz, double qw)
{
    rotAngle = Angle{0.0, 0.0, 0.0};
    SetRotateMatrixFromQuaternion(qx, qy, qz, qw);
    ApplyRotateMatrix();
}

void Particle::ApplyRotateMatrix()
{
    // REF: слить всё в один цикл
    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < facets[i].nVertices; ++j)
        {
            RotatePoint(defaultFacets[i].arr[j], facets[i].arr[j]);
        }
    }

    RotateNormals();
    RotateCenters();
}

void Particle::Fix()
{
    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < facets[i].nVertices; ++j)
        {
            defaultFacets[i].arr[j] = facets[i].arr[j];
        }
    }
}

void Particle::Scale(double ratio)
{
    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < defaultFacets[i].nVertices; ++j)
        {
            defaultFacets[i].arr[j].cx *= ratio;
            defaultFacets[i].arr[j].cy *= ratio;
            defaultFacets[i].arr[j].cz *= ratio;
        }
        defaultFacets[i].center.cx *= ratio;
        defaultFacets[i].center.cy *= ratio;
        defaultFacets[i].center.cz *= ratio;
        defaultFacets[i].in_normal.d_param *= ratio;
        defaultFacets[i].ex_normal.d_param *= ratio;
    }
}

// vector<unsigned> Particle::FindFlippedFacets() const
// {
//     vector<unsigned> facets;

//     if (isAggregated)
//     {
//         vector<Particle> parts;
//         DisconnectParts(parts);

//         unsigned totalNFacets = 0;

//         for (unsigned i = 0; i < parts.size(); ++i)
//         {
//             vector<unsigned> partFacets;
//             parts[i].FindFlippedFacets(partFacets);

//             for (unsigned j = 0; j < partFacets.size(); ++j)
//             {
//                 partFacets[j] += totalNFacets;
//                 facets.push_back(partFacets[j]);
//             }

//             totalNFacets += parts[i].nFacets;
//         }
//     }
//     else
//     {
//         FindFlippedFacets(facets);
//     }

//     return facets;
// }

// void Particle::DisconnectParts(std::vector<Particle> &parts) const
// {
//     for (unsigned i = 0; i < partIndices.size(); ++i)
//     {
//         Particle part;

//         for (unsigned j = partIndices[i].from; j < partIndices[i].to; ++j)
//         {
//             part.AddFacet(facets[j]);
//         }

//         part.UpdateFacets();
//         parts.push_back(part);
//     }
// }

void Particle::Resize(double size)
{
    double dmax = MaximalDimention();
    double ratio = size / dmax;
    Scale(ratio);
    double actual = MaximalDimention();
    if (actual > DBL_EPSILON)
        Scale(size / actual);
    Reset();
}

Facet Particle::CreateBase(int nVertices, double radius, double zCoordinate)
{
    Facet base;
    double step = M_2PI/nVertices;

    int begin;
    int end;
    int inc;

    if (zCoordinate >= 0)
    {
        begin = 0;
        end = nVertices;
        inc = 1;
        base.normal[0] = Point3f(0, 0, -1);
        base.normal[1] = Point3f(0, 0, 1);
    }
    else  // reverse vertex order
    {
        begin = nVertices-1;
        end = -1;
        inc = -1;
        base.normal[0] = Point3f(0, 0, 1);
        base.normal[1] = Point3f(0, 0, -1);
    }

    for (int i = begin; i != end; i += inc)
    {
        double a = step*i;
        base.AddVertex(Point3f(radius*cos(a), radius*sin(a), zCoordinate));
    }

    return base;
}

void Particle::CreateSideFacets(const Facet &top, const Facet &bottom,
                                int &iFacet)
{
    std::vector<Facet> sideFacets;
    sideFacets = ConnectBases(top, bottom);

    for (size_t i = 0; i < sideFacets.size(); ++i)
    {
        defaultFacets[++iFacet] = sideFacets[i];
    }
}

std::vector<Facet> Particle::ConnectBases(const Facet &top,
                                          const Facet &bottom)
{
    std::vector<Facet> faces;

    if (top.nVertices == bottom.nVertices)
    {
        bool isInverse = DotProduct(bottom.normal[1], top.normal[1]) < 0;
        int end = top.nVertices-1;

        int i1 = end;
        int i2 = 0;

        if (top.normal[1].coordinates[2] >= 0)
        {
            for (int i = 0; i < top.nVertices; ++i)
            {
                Facet f;

                f.AddVertex(top.arr[i2]);
                f.AddVertex(top.arr[i1]);

                if (isInverse)
                {
                    f.AddVertex(bottom.arr[end-i1]);
                    f.AddVertex(bottom.arr[end-i2]);
                }
                else
                {
                    f.AddVertex(bottom.arr[i1]);
                    f.AddVertex(bottom.arr[i2]);
                }

                faces.push_back(f);

                i1 = i2;
                ++i2;
            }
        }
        else
        {
            for (int i = 0; i < top.nVertices; ++i)
            {
                Facet f;

                f.AddVertex(top.arr[i1]);
                f.AddVertex(top.arr[i2]);

                if (isInverse)
                {
                    f.AddVertex(bottom.arr[end-i2]);
                    f.AddVertex(bottom.arr[end-i1]);
                }
                else
                {
                    f.AddVertex(bottom.arr[i2]);
                    f.AddVertex(bottom.arr[i1]);
                }

                faces.push_back(f);

                i1 = i2;
                ++i2;
            }
        }
    }
    else
    {
        throw "Number of base vertices are not equal";
    }

    return faces;
}

void Particle::FindFlippedFacets(std::vector<unsigned> &facetIndices) const
{
    Point3f globalCenter = Center();

    for (int i = 0; i < nFacets; ++i)
    {
        if (!facets[i].IsConormal(facets[i].center - globalCenter))
        {
            facetIndices.push_back(i);
        }
    }
}

void Particle::Concate(const std::vector<Particle> &parts)
{
    int i = 0;
    nFacets = 0;

    for (const Particle &part : parts)
    {
        nFacets += part.nFacets;

        for (int j = 0; j < part.nFacets; ++j)
        {
            defaultFacets[i++] = part.facets[j];
        }
    }

    isAggregated = !parts.empty();
    nFacetsInPart = parts.empty() ? 0 : parts.front().nFacets;
    for (const Particle &part : parts)
    {
        if (part.nFacets != nFacetsInPart)
        {
            nFacetsInPart = 0;
            break;
        }
    }
}

double Particle::LongRadius() const
{
    Point3f p0(0, 0, 0);

    double radius = 0;

    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < facets[i].nVertices; ++j)
        {
            Point3f v_len = facets[i].arr[j] - p0;
            double len = Length(v_len);

            if (len > radius)
            {
                radius = len;
            }
        }
    }

    return radius;
}

double Particle::MaximalDimention() const
{
    std::vector<Point3f> vertices;
    size_t vertexCount = 0;
    for (int i = 0; i < nFacets; ++i)
        vertexCount += defaultFacets[i].nVertices;
    vertices.reserve(vertexCount);

    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < defaultFacets[i].nVertices; ++j)
        {
            vertices.push_back(defaultFacets[i].arr[j]);
        }
    }
    return MaximumPairDistance(vertices);
}

double Particle::MaximalDimentionPart() const
{
    std::vector<Point3f> vertices;
    size_t vertexCount = 0;
    for (int i = 0; i < nFacetsInPart; ++i)
        vertexCount += defaultFacets[i].nVertices;
    vertices.reserve(vertexCount);

    for (int i = 0; i < nFacetsInPart; ++i)
    {
        for (int j = 0; j < defaultFacets[i].nVertices; ++j)
            vertices.push_back(defaultFacets[i].arr[j]);
    }
    return MaximumPairDistance(vertices);
}

double Particle::MaximumEdgeLength() const
{
    double maxLength = 0.0;

    for (int i = 0; i < nFacets; ++i)
    {
        const Facet &facet = defaultFacets[i];
        if (facet.nVertices < 2)
            continue;

        for (int j = 0; j < facet.nVertices; ++j)
        {
            const int next = (j + 1) % facet.nVertices;
            maxLength = std::max(
                maxLength, Length(facet.arr[next] - facet.arr[j]));
        }
    }

    return maxLength;
}

double Particle::Area()
{
    double area = 0;

    for (int i = 0; i < nFacets; ++i)
    {
        area += facets[i].Area();
    }

    return area;
}

Point3f Particle::Center() const
{
    Point3f center = Point3f(0, 0, 0);
    int nVertices = 0;

    for (int i = 0; i < nFacets; ++i)
    {
        auto &facet = defaultFacets[i];
        nVertices += facet.nVertices;

        for (int j = 0; j < facet.nVertices; ++j)
        {
            center = center + facet.arr[j];
        }
    }

    center = center/nVertices;
    return center;
}

double Particle::Volume()
{
    double volume = 0;

    for (int i = 0; i < nFacets; ++i)
    {
        const Facet &facet = defaultFacets[i];
        if (facet.nVertices < 3)
            continue;

        const Point3f &p0 = facet.arr[0];
        for (int j = 1; j + 1 < facet.nVertices; ++j)
        {
            Point3f cross = CrossProduct(facet.arr[j], facet.arr[j + 1]);
            volume += DotProduct(p0, cross) / 6.0;
        }
    }

    return fabs(volume);
}

const ::complex &Particle::GetRefractiveIndex() const
{
    return m_refractiveIndex;
}

const Symmetry &Particle::GetSymmetry() const
{
    return m_symmetry;
}

void Particle::Move(double dx, double dy, double dz)
{
    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < defaultFacets[i].nVertices; ++j)
        {
            facets[i].arr[j] = defaultFacets[i].arr[j] + Point3f(dx, dy, dz);
        }
    }
}

bool Particle::IsConcave() const
{
    return isConcave;
}

void Particle::Output(const std::string &name) const
{
    if (isAggregated
        && (nFacetsInPart <= 0 || nFacets % nFacetsInPart != 0))
    {
        throw std::runtime_error(
            "cannot save aggregate geometry to '" + name
            + "': FACETS_PER_PART is missing or inconsistent.\n"
              "  Fix: define an equal positive facet count for every aggregate component.");
    }

    std::ofstream file(name, std::ios::out);
    if (!file.is_open())
    {
        throw std::runtime_error(
            "cannot save particle geometry to '" + name + "'.\n"
            "  Fix: create its parent directory and verify write permissions and free disk space.");
    }

    file << std::setprecision(17);

    file << (int)isConcave << std::endl;
    file << (int)isAggregated;

    if (isAggregated)
    {
        file << ' ' << nFacetsInPart;
    }

    file << std::endl;

    file << RadToDeg(m_symmetry.beta) << ' '
           << RadToDeg(m_symmetry.gamma)
           << std::endl << std::endl;

    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < facets[i].nVertices; ++j)
        {
            const Point3f &p = facets[i].arr[j];
            file << p.coordinates[0] << ' '
                            << p.coordinates[1] << ' '
                            << p.coordinates[2] << ' '
                            /*<< i */;
            file << std::endl;
        }

        if ((i + 1) != nFacets)
        {
            file << std::endl;
        }
    }

    file.close();
    if (!file)
    {
        throw std::runtime_error(
            "failed while writing particle geometry to '" + name + "'.\n"
              "  Fix: verify free disk space and that the destination remains writable.");
    }

    const std::string visibilityPath = name + ".visibility";
    std::ofstream visibilityFile(visibilityPath.c_str(), std::ios::out);
    if (!visibilityFile.is_open())
    {
        throw std::runtime_error(
            "cannot save facet visibility metadata to '" + visibilityPath
            + "'.\n  Fix: verify write permissions for the geometry directory.");
    }
    visibilityFile << nFacets << '\n';
    for (int i = 0; i < nFacets; ++i)
    {
        visibilityFile << static_cast<int>(facets[i].isVisibleIn) << ' '
                       << static_cast<int>(facets[i].isVisibleOut) << '\n';
    }
    visibilityFile.close();
    if (!visibilityFile)
    {
        throw std::runtime_error(
            "failed while writing facet visibility metadata to '"
            + visibilityPath
            + "'.\n  Fix: verify free disk space and destination permissions.");
    }
}

void Particle::SetRefractiveIndex(const ::complex &value)
{
    m_refractiveIndex = value;
}

void Particle::SetDefaultNormals()
{
    for (int i = 0; i < nFacets; ++i)
    {
        defaultFacets[i].SetNormal();
    }
}

void Particle::Reset()
{
    for (int i = 0; i < nFacets; ++i)
    {
        facets[i] = defaultFacets[i];
    }
}

void Particle::SetDefaultCenters()
{
    for (int i = 0; i < nFacets; ++i)
    {
        defaultFacets[i].SetCenter();
    }
}

void Particle::SetRotateMatrix(double beta, double gamma, double alpha)
{
    double cosA, cosB, cosG,
            sinA, sinB, sinG;

    sincos(alpha, &sinA, &cosA);
    sincos(beta,  &sinB, &cosB);
    sincos(gamma, &sinG, &cosG);

    double cosAcosB = cosA*cosB;
    double sinAcosG = sinA*cosG;
    double sinAsinG = sinA*sinG;

    m_rotMatrix[0][0] = cosAcosB*cosG - sinAsinG;
    m_rotMatrix[1][0] = sinAcosG*cosB + cosA*sinG;
    m_rotMatrix[2][0] = -sinB*cosG;

    m_rotMatrix[0][1] = -(cosAcosB*sinG + sinAcosG);
    m_rotMatrix[1][1] = cosA*cosG - sinAsinG*cosB;
    m_rotMatrix[2][1] = sinB*sinG;

    m_rotMatrix[0][2] = cosA*sinB;
    m_rotMatrix[1][2] = sinA*sinB;
    m_rotMatrix[2][2] = cosB;
}

void Particle::SetRotateMatrixFromQuaternion(double qx, double qy, double qz, double qw)
{
    double norm = sqrt(qx*qx + qy*qy + qz*qz + qw*qw);
    if (norm <= DBL_EPSILON)
    {
        qx = qy = qz = 0.0;
        qw = 1.0;
    }
    else
    {
        qx /= norm;
        qy /= norm;
        qz /= norm;
        qw /= norm;
    }

    const double xx = qx*qx, yy = qy*qy, zz = qz*qz;
    const double xy = qx*qy, xz = qx*qz, yz = qy*qz;
    const double wx = qw*qx, wy = qw*qy, wz = qw*qz;

    m_rotMatrix[0][0] = 1.0 - 2.0*(yy + zz);
    m_rotMatrix[0][1] = 2.0*(xy - wz);
    m_rotMatrix[0][2] = 2.0*(xz + wy);

    m_rotMatrix[1][0] = 2.0*(xy + wz);
    m_rotMatrix[1][1] = 1.0 - 2.0*(xx + zz);
    m_rotMatrix[1][2] = 2.0*(yz - wx);

    m_rotMatrix[2][0] = 2.0*(xz - wy);
    m_rotMatrix[2][1] = 2.0*(yz + wx);
    m_rotMatrix[2][2] = 1.0 - 2.0*(xx + yy);
}

void Particle::RemoveWalls()
{
    for (int i = 0; i < nFacets; ++i)
    {
        auto &f1 = defaultFacets[i];

        for (int j = i+1; j < nFacets; ++j)
        {
            auto &f2 = defaultFacets[j];

            bool isFoundF = true;

            if (f1.nVertices == f2.nVertices)
            {
                for (int m = 0; m < f1.nVertices; ++m)
                {
                    bool isFoundV = false;

                    for (int k = 0; k < f2.nVertices; ++k)
                    {
                        if (f1.arr[m].IsEqualTo(f2.arr[k], 0.005))
                        {
                            isFoundV = true;
                            break;
                        }
                    }

                    if (!isFoundV)
                    {
                        isFoundF = false;
                        break;
                    }
                }

                if (isFoundF)
                {
                    defaultFacets[i].nVertices = 0;
                    defaultFacets[j].nVertices = 0;
                    break;
                }
            }
        }
    }

    for (int i = 0; i < nFacets; ++i)
    {
        if (defaultFacets[i].nVertices == 0)
        {
            for (int j = i+1; j < nFacets; ++j)
            {
                defaultFacets[j-1] = defaultFacets[j];
            }

            --nFacets;
            --i;
        }
    }
}

void Particle::RotateNormals()
{
    for (int i = 0; i < nFacets; ++i)
    {
        RotatePoint(defaultFacets[i].in_normal, facets[i].in_normal);
    }

    SetDParams();

    for (int i = 0; i < nFacets; ++i)
    {
        facets[i].ex_normal = -facets[i].in_normal;
        facets[i].ex_normal.d_param = -facets[i].in_normal.d_param;
    }
}

void Particle::SetDParams()
{
    for (int i = 0; i < nFacets; ++i)
    {
        double d = DotProduct(facets[i].arr[0], facets[i].in_normal);
        facets[i].in_normal.d_param = -d;
    }
}

void Particle::RotatePoint(const Point3f &point, Point3f &result)
{
    result.cx = point.cx*m_rotMatrix[0][0] + point.cy*m_rotMatrix[0][1] + point.cz*m_rotMatrix[0][2];
    result.cy = point.cx*m_rotMatrix[1][0] + point.cy*m_rotMatrix[1][1] + point.cz*m_rotMatrix[1][2];
    result.cz = point.cx*m_rotMatrix[2][0] + point.cy*m_rotMatrix[2][1] + point.cz*m_rotMatrix[2][2];
}

void Particle::SetSymmetry(double beta, double gamma, double alpha)
{
    m_symmetry.beta = beta;
    m_symmetry.gamma = gamma;
    m_symmetry.alpha = alpha;
}
