#include "Polygon.h"
#include <math.h>
#include <algorithm>
#include <cfloat>

#define EPS_NORMAL 0.1

Polygon::Polygon(const Polygon &other)
{
	if (other.nVertices < 0 || other.nVertices > MAX_VERTEX_NUM)
		throw std::runtime_error("invalid polygon vertex count during copy");
	nVertices = other.nVertices;
	std::copy_n(other.arr, nVertices, arr);
}

Polygon::Polygon(Polygon &&other)
{
	if (other.nVertices < 0 || other.nVertices > MAX_VERTEX_NUM)
		throw std::runtime_error("invalid polygon vertex count during move");
	nVertices = other.nVertices;
	std::copy_n(other.arr, nVertices, arr);

	other.nVertices = 0;
}

void Polygon::AddVertex(const Point3f &v)
{
	if (nVertices >= MAX_VERTEX_NUM)
	{
		throw std::runtime_error(
			"polygon exceeds the 64-vertex geometry limit");
	}
	arr[nVertices++] = v;
}

Polygon &Polygon::operator =(const Polygon &other)
{
	if (this != &other)
	{
		if (other.nVertices < 0 || other.nVertices > MAX_VERTEX_NUM)
			throw std::runtime_error("invalid polygon vertex count during assignment");
		nVertices = other.nVertices;
		std::copy_n(other.arr, nVertices, arr);
	}

	return *this;
}

Polygon &Polygon::operator = (Polygon &&other)
{
	if (this != &other)
	{
		if (other.nVertices < 0 || other.nVertices > MAX_VERTEX_NUM)
			throw std::runtime_error("invalid polygon vertex count during move assignment");
		nVertices = other.nVertices;
		std::copy_n(other.arr, nVertices, arr);

		other.nVertices = 0;
	}

	return *this;
}

std::ostream &operator <<(std::ostream &os, const Polygon &beam)
{
	using namespace std;

	os << "polygon: {" << endl;

	for (int i = 0; i < beam.nVertices; ++i)
	{
		os << "\t" << i << ": "
		   << beam.arr[i].cx << ", "
		   << beam.arr[i].cy << ", "
		   << beam.arr[i].cz << ", "
		   << beam.arr[i].d_param << endl;
	}

	os << "}" << endl << endl;
	return os;
}

double Polygon::Area() const
{
    if (nVertices < MIN_VERTEX_NUM)
        return 0.0;

    double square = 0.0;
    const Point3f &basePoint = arr[0];
    double p1x = static_cast<double>(arr[1].cx) - basePoint.cx;
    double p1y = static_cast<double>(arr[1].cy) - basePoint.cy;
    double p1z = static_cast<double>(arr[1].cz) - basePoint.cz;

    for (int i = 2; i < nVertices; ++i)
    {
        const double p2x = static_cast<double>(arr[i].cx) - basePoint.cx;
        const double p2y = static_cast<double>(arr[i].cy) - basePoint.cy;
        const double p2z = static_cast<double>(arr[i].cz) - basePoint.cz;
        const double crossX = p1y*p2z - p1z*p2y;
        const double crossY = p1z*p2x - p1x*p2z;
        const double crossZ = p1x*p2y - p1y*p2x;
        square += std::sqrt(crossX*crossX + crossY*crossY
                           + crossZ*crossZ);
        p1x = p2x;
        p1y = p2y;
        p1z = p2z;
	}

	return square/2.0;
}

Point3f Polygon::Center() const
{
    if (nVertices <= 0)
        return Point3f(0, 0, 0);

    double x = 0.0, y = 0.0, z = 0.0;

    for (int i = 0; i < nVertices; ++i)
    {
        x += arr[i].cx;
        y += arr[i].cy;
        z += arr[i].cz;
    }

    const double inv = 1.0/nVertices;
    return Point3f(x*inv, y*inv, z*inv);
}

Point3f Polygon::Normal() const
{
    if (nVertices < MIN_VERTEX_NUM)
        return Point3f(0, 0, 0);

    // Newell's formula uses every edge and is much less sensitive than
    // normalizing the first non-collinear triangle.
    double nx = 0.0, ny = 0.0, nz = 0.0;
    for (int i = 0; i < nVertices; ++i)
    {
        const Point3f &a = arr[i];
        const Point3f &b = arr[(i + 1) % nVertices];
        nx += (static_cast<double>(a.cy) - b.cy)*(a.cz + b.cz);
        ny += (static_cast<double>(a.cz) - b.cz)*(a.cx + b.cx);
        nz += (static_cast<double>(a.cx) - b.cx)*(a.cy + b.cy);
    }
    const double length = std::sqrt(nx*nx + ny*ny + nz*nz);
    if (!(length > DBL_MIN) || !std::isfinite(length))
        return Point3f(0, 0, 0);
    return Point3f(nx/length, ny/length, nz/length);
}
