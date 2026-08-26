#pragma once

#include "Beam.h"

#include <cfloat>

#define EPS_COS_00		0.99999999998254670756866631966593		//1 - cos(89.99999999)
// Reject only numerically tangent incidence. Directions and normals are
// double precision, so the former float-sized angular dead zone is too wide.
constexpr double EPS_GRAZING_INCIDENCE = 0x1p-40;

inline double SignedPhasePlaneDistance(const Point3f &direction,
                                       double front,
                                       const Point3f &point)
{
    return static_cast<double>(direction.cx)*point.cx
         + static_cast<double>(direction.cy)*point.cy
         + static_cast<double>(direction.cz)*point.cz + front;
}

class Splitting
{
public:
	Splitting(bool isOpticalPath);
	void ComputeRiParams(const complex &ri);

	void ComputeCosA(const Point3f &normal, const Point3f &incidentDir);
	void ComputeSplittingParams(const Point3f &dir, const Point3f &normal);

	bool IsCompleteReflection();
	bool IsNormalIncidence();
	bool IsIncident();

	double ComputeEffectiveReRi() const;

	double ComputeIncidentOpticalPath(const Point3f &direction,
									  const Point3f &facetPoint);
	double ComputeOutgoingOpticalPath(const Beam &beam);
	double ComputeSegmentOpticalPath(const Beam &beam,
									 const Point3f &facetPoint) const;

	void ComputeCRBeamParams(const Point3f &normal, const Beam &incidentBeam,
							 Beam &inBeam);

	void ComputeNormalBeamParams(const Beam &incidentBeam,
								 Beam &inBeam, Beam &outBeam);
	void ComputeNormalBeamParamsExternal(const Light &incidentLight,
										 Beam &inBeam, Beam &outBeam);

	void ComputeRegularBeamsParams(const Point3f &normal,
								   const Beam &incidentBeam,
								   Beam &inBeam, Beam &outBeam);
	void ComputeRegularBeamParamsExternal(const Point3f &facetNormal,
										  Beam &incidentBeam,
										  Beam &inBeam, Beam &outBeam);

	Point3f ChangeBeamDirection(const Vector3f &oldDir, const Vector3f &normal,
								Location oldLoc, Location loc);
private:
	Point3f r{}; // tangential component of the incident direction
	double reRiEff;
	double s;
//	double cosA;
	bool m_isOpticalPath;

	complex m_ri;	//  refractive index
	double m_cRiRe;
	double m_cRiRe2;
	double m_cRiIm;

public:
	double cosA;

	complex GetRi() const;

private:
	void ComputeCRJonesParams(complex &cv, complex &ch);

	void ComputeRegularJonesParams(const Point3f &normal,
								   const Beam &incidentBeam,
								   Beam &inBeam, Beam &outBeam);
	void ComputeInternalRefractiveDirection(const Vector3f &tangential,
											const Vector3f &normal, Vector3f &dir);

};
