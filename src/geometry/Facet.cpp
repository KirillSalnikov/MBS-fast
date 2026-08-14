#include "Facet.h"
#include "Intersection.h"

void Facet::SetNormal()
{
	ex_normal = Normal();
	in_normal = -ex_normal;
}

void Facet::SetCenter()
{
    center = Center();
}

bool Facet::IsConormal(Point3f normal) const
{
    return DotProduct(this->normal[Out], normal) > EPS_PROJECTION;
}

Facet &Facet::operator =(const Facet &other)
{
	if (this != &other)
	{
		Polygon::operator =(other);
		in_normal = other.in_normal;
		ex_normal = other.ex_normal;
		center = other.center;
		isVisibleIn = other.isVisibleIn;
		isVisibleOut = other.isVisibleOut;
	}

	return *this;
}
