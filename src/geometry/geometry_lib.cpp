#include <math.h>
#include "Particle.h"
#include "Intersection.h"

double DotProductD(const Vector3d &v1, const Vector3d &v2)
{
	return	  v1.x * v2.x
			+ v1.y * v2.y
			+ v1.z * v2.z;
}

double Norm(const Vector3f &p)
{
	return	  (double)p.cx * (double)p.cx
			+ (double)p.cy * (double)p.cy
			+ (double)p.cz * (double)p.cz;
}

double NormD(const Vector3d &p)
{
	return	  p.x * p.x
			+ p.y * p.y
			+ p.z * p.z;
}

void CrossProduct(const Vector3f &v1, const Vector3f &v2, Vector3f &res)
{
	res.cx = v1.cy*v2.cz - v1.cz*v2.cy;
	res.cy = v1.cz*v2.cx - v1.cx*v2.cz;
	res.cz = v1.cx*v2.cy - v1.cy*v2.cx;
	res.d_param = 0.0;
}

Point3f IntersectVectors(const Point3f &c1, const Point3f &c2,
						 const Point3f &v1, const Point3f &v2,
						 const Point3f &normalToFacet, bool &isOk)
{
	return intersect_iv(c1, c2, v1, v2, normalToFacet, isOk);
}

// REF: try to move to Point3f
Point3f CrossProduct(const Point3f &v1, const Point3f &v2)
{
	return Point3f(v1.cy*v2.cz - v1.cz*v2.cy,
				   v1.cz*v2.cx - v1.cx*v2.cz,
				   v1.cx*v2.cy - v1.cy*v2.cx);
}

std::ostream &operator <<(std::ostream &os, const Point3f &p)
{
	os << p.coordinates[0]
			<< " " << p.coordinates[1]
			<< " " << p.coordinates[2];
	return os;
}

// OPT:
Point3d CrossProductD(const Point3d &v1, const Point3d &v2)
{
	Point3d res;
	res.x = v1.y*v2.z - v1.z*v2.y;
	res.y = v1.z*v2.x - v1.x*v2.z;
	res.z = v1.x*v2.y - v1.y*v2.x;
	return res;
}

double Length(const Vector3f &v)
{
	return sqrt(Norm(v));
}

double LengthD(const Vector3d &v)
{
	return sqrt(NormD(v));
}

void Normalize(Vector3f &v)
{
	const double length = sqrt(Norm(v));
	if (!(length > DBL_MIN) || !std::isfinite(length))
	{
		v = Vector3f(0, 0, 0);
		return;
	}
	v.cx /= length;
	v.cy /= length;
	v.cz /= length;
}

Vector3d NormalizeD(const Vector3d &v)
{
	Vector3d res;
	const double length = sqrt(NormD(v));
	if (!(length > DBL_MIN) || !std::isfinite(length))
		return Vector3d(0, 0, 0);
	res.x = v.x / length;
	res.y = v.y / length;
	res.z = v.z / length;
	return res;
}

bool ProjectPointToPlane(const Point3f &point, const Vector3f &direction,
						 const Vector3f &planeNormal, Point3f &projection)
{
	double tmp = DotProduct(point, planeNormal);
	double dp  = DotProduct(direction, planeNormal);
	const double scale = Length(direction)*Length(planeNormal);
	if (scale <= DBL_MIN
		|| std::fabs(dp) <= geometry_parallel_tolerance(scale))
		return false;
	tmp += planeNormal.d_param;
	tmp /= dp;
	const double x = static_cast<double>(point.cx) - direction.cx*tmp;
	const double y = static_cast<double>(point.cy) - direction.cy*tmp;
	const double z = static_cast<double>(point.cz) - direction.cz*tmp;
	if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
		return false;
	projection = Point3f(x, y, z);
	return true;
}
