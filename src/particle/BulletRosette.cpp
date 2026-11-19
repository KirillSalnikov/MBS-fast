#include "BulletRosette.h"
#include "Bullet.h"
#include "global.h"

BulletRosette::BulletRosette()
{
	SetSymmetry(M_PI/2, M_PI);
}

BulletRosette::BulletRosette(const complex &refrIndex, double diameter,
							 double height, double peakHeight, bool joinedTips)
	: BulletRosette()
{
	isConcave = true;
	isAggregated = true;

	Init(8, refrIndex);

	std::vector<Particle> bullets;
	const double centerShift = height / 2 + peakHeight + (joinedTips ? 0.0 : 1.0);

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -centerShift);
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -centerShift);
		b.Fix();
		b.Rotate(DegToRad(90), 0, 0);
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -centerShift);
		b.Fix();
		b.Rotate(DegToRad(180), 0, 0);
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -centerShift);
		b.Fix();
		b.Rotate(DegToRad(270), 0, 0);
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -centerShift);
		b.Fix();
		b.Rotate(DegToRad(90), 0, 0);
		b.Fix();
		b.Rotate(0, 0, DegToRad(90));
		b.Fix();
		bullets.push_back(b);
	}

	{
		Bullet b(refrIndex, diameter, height, peakHeight);
		b.Move(0, 0, -centerShift);
		b.Fix();
		b.Rotate(DegToRad(90), 0, 0);
		b.Fix();
		b.Rotate(0, 0, DegToRad(270));
		b.Fix();
		bullets.push_back(b);
	}

	Concate(bullets);

	SetDefaultNormals();
	Reset();
	SetDParams();
	SetDefaultCenters();

	for (int i = 0; i < nFacets; ++i)
	{
		defaultFacets[i].isVisibleIn = false;
		defaultFacets[i].isVisibleOut = false;
		facets[i].isVisibleIn = false;
		facets[i].isVisibleOut = false;
	}
}

void BulletRosette::GetParticalFacetIdRangeByFacetId(int id, int &begin, int &end) const
{
	int patN = id/13;
	begin = patN*13;
	end = begin+13;
}
