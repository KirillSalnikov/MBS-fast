#include "Particle.h"
#include <cfloat>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <vector>
#include "global.h"

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

    std::vector<std::vector<Point3f> > parsedFacets;
    std::vector<Point3f> vertices;
    int facetStartLine = 0;
    const auto finishFacet = [&](int separatorLine) {
        if (vertices.empty())
            return;
        if (vertices.size() < 3)
            ParticleFileError(filename, facetStartLine,
                              "a facet has fewer than three vertices.",
                              "add vertices or remove the incomplete facet before line "
                                  + std::to_string(separatorLine) + ".");
        if (parsedFacets.size() >= MAX_FACET_NUM)
            ParticleFileError(filename, facetStartLine,
                              "the particle exceeds the facet limit.",
                              "reduce the geometry to at most "
                                  + std::to_string(MAX_FACET_NUM) + " facets.");

        Facet candidate;
        for (const Point3f &point : vertices)
            candidate.AddVertex(point);
        const double area = candidate.Area();
        if (!(area > 0.0) || !std::isfinite(area))
            ParticleFileError(filename, facetStartLine,
                              "a facet has zero or non-finite area.",
                              "remove duplicate/collinear vertices and keep a nondegenerate polygon.");
        parsedFacets.push_back(vertices);
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
        if (vertices.size() >= MAX_VERTEX_NUM - 1)
            ParticleFileError(filename, lineNumber,
                              "a facet exceeds the vertex limit.",
                              "split or simplify it to at most "
                                  + std::to_string(MAX_VERTEX_NUM - 1) + " vertices.");
        if (vertices.empty())
            facetStartLine = lineNumber;
        vertices.push_back(Point3f(
            ParseDoubleToken(coordinates[0], filename, lineNumber, "X coordinate"),
            ParseDoubleToken(coordinates[1], filename, lineNumber, "Y coordinate"),
            ParseDoubleToken(coordinates[2], filename, lineNumber, "Z coordinate")));
    }
    if (pfile.bad())
        ParticleFileError(filename, lineNumber,
                          "an I/O error interrupted the particle body.",
                          "verify the file storage and permissions, then read the complete file again.");
    finishFacet(lineNumber + 1);

    if (parsedFacets.empty())
        ParticleFileError(filename, lineNumber,
                          "the file contains no facets.",
                          "add facet polygons after the three header records; separate facets with blank lines.");
    if (aggregateValue == 1
        && parsedFacets.size() % static_cast<size_t>(parsedFacetsPerPart) != 0)
    {
        ParticleFileError(filename, aggregateLine.number,
                          "the total facet count is not divisible by FACETS_PER_PART.",
                          "correct FACETS_PER_PART or complete every aggregate component.");
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

    for (unsigned i = 0; i < nFacets; ++i)
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

    isAggregated = true;
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
    double Dmax = 0;
    double newDmax;

    std::vector<Point3f> vertices;

    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < defaultFacets[i].nVertices; ++j)
        {
            vertices.push_back(defaultFacets[i].arr[j]);
        }
    }

    for (int i = 0; i < vertices.size(); ++i)
    {
        for (int j = 0; j < vertices.size(); ++j)
        {
            newDmax = Length(vertices[j] - vertices[i]);

            if (newDmax > Dmax)
            {
                Dmax = newDmax;
            }
        }
    }

    return Dmax;
}

double Particle::MaximalDimentionPart() const
{
    double Dmax = 0;
    double newDmax;

    Polygon512 pol;

    for (int i = 0; i < nFacetsInPart; ++i)
    {
        pol.Concat(defaultFacets[i]);
    }

    for (int i = 0; i < pol.nVertices; ++i)
    {
        for (int j = 0; j < pol.nVertices; ++j)
        {
            newDmax = Length(pol.arr[j] - pol.arr[i]);

            if (newDmax > Dmax)
            {
                Dmax = newDmax;
            }
        }
    }

    return Dmax;
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

void Particle::Move(float dx, float dy, float dz)
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

void Particle::Output(std::string name)
{
    std::ofstream file(name, std::ios::out);

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
            Point3f p = facets[i].arr[j];
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
