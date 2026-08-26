#pragma once

#include "Tracks.h"
#include "Beam.h"
#include "Particle.h"
#include "Intersection.h"
#include "Splitting.h"
#include "CutoffStatistics.h"

#include <float.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

class TraceLimitExceeded : public std::runtime_error
{
public:
    explicit TraceLimitExceeded(const std::string &message)
        : std::runtime_error(message) {}
};

class Scattering
{
public:
    // REF: убрать в протектед потом
    Particle *m_particle;			///< scattering particle (crystal)
    double m_wave;
    double restriction = 100;

protected:
    Facet *m_facets;
    double m_geometryScale;
    double m_geometryLengthTolerance;
    double m_geometryAreaTolerance;
    double m_geometryDepthTolerance;
    bool m_unitFacetNormals;
    Splitting splitting;
    Light *m_incidentLight;

    Vector3f m_incidentDir;
    Vector3f m_polarBasis;
    int m_nActs;

    std::vector<Beam> m_beamTree;	///< tree of beams (works like stack), grows lazily
    int m_treeSize;
    double m_incidentEnergy;
    double m_traceRefJNorm = 0;
    double m_traceRefArea = 0;
    double m_traceRefImportance = 0;
    bool m_trackIdsRequired = false;
    std::shared_ptr<TraceCutoffStatistics> m_traceCutoffStatistics;

public:
    Scattering(Particle *particle, Light *incidentLight, bool isOpticalPath,
               int nActs);
    virtual ~Scattering();

    /// Create a clone that uses a different particle copy.
    /// Caller owns the returned pointer.
    virtual Scattering* CloneFor(Particle *newParticle, Light *light) {
        Scattering *copy = new Scattering(newParticle, light, true, m_nActs);
        copy->CopyRuntimeOptionsFrom(*this);
        return copy;
    }
    virtual void PrepareForParallelTrace() {}
    void SetMaxReflections(int n) { m_nActs = n; }
    int GetMaxReflections() const { return m_nActs; }
    double m_traceCutoffJRel = 0;
    double m_traceCutoffAreaRel = 0;
    double m_traceCutoffImportanceRel = 0;
    int m_traceMaxBeams = 0;
    int m_traceLimitRetries = 0;
    double m_traceRetryFactor = 10.0;
    std::string m_cutoffProfileName = "legacy-default";
    bool m_gpuTracePrefilter = false;
    bool m_traceCpuProjectedPrefilter = true;
    double m_traceCpuProjectedPrefilterMargin = -1.0;
    bool m_tracePrefilterStats = false;
    void SetTrackIdsRequired(bool value) { m_trackIdsRequired = value; }
    void CopyRuntimeOptionsFrom(const Scattering &source);
    std::string TraceCutoffReport() const;

    virtual bool ScatterLight(double /*beta*/, double /*gamma*/, std::vector<Beam> &/*scaterredBeams*/)
    {
        return false;
    }

    virtual bool ScatterLight(double beta, double gamma, const std::vector<std::vector<int>> &tracks,
                                     std::vector<Beam> &scaterredBeams);

    bool ScatterLightWithLimitRetry(double beta, double gamma,
                                    std::vector<Beam> &scatteredBeams,
                                    bool formShadow);

    virtual void FormShadowBeam(std::vector<Beam> &scaterredBeams);

    double GetIncedentEnergy() const;

    double ComputeInternalOpticalPath(const Beam &beam, const Point3f sourcePoint,
                                      const std::vector<int> &track);
//	double CrossSection(const Point3f &beamDir) const;

protected:
    void SetIncidentBeamOpticalParams(unsigned facetId, Beam &inBeam, Beam &outBeam);

    void Difference(const Polygon &subject, const Vector3f &subjNormal,
                    const Polygon &clip, const Vector3f &clipNormal,
                    const Vector3f &clipDir, PolygonArray &difference,
                    Polygon *overlap = nullptr) const;

    void Intersect(int facetId, const Beam& beam, Polygon &intersection) const;

    void SetPolygonByFacet(int facetId, Polygon &polygon) const;

    bool IsTerminalAct(const Beam &beam);
    bool IsTracePruned(const Beam &beam) const;
    bool EnsureBeamTree();
    void ResetTraceReference();
    void UpdateTraceReference(const Beam &beam);

    void SplitLightToBeams(int facetId, Beam &inBeam, Beam &outBeam);

    void ComputePolarisationParams(const Vector3f &dir,
                                   const Vector3f &facetNormal, Beam &beam);

    void AddProjectedIncidentEnergy(int facetId, const Polygon &lightedPolygon);


    bool PushBeamToTree(Beam &beam, int facetId, int level, Location location);


    IdType RecomputeTrackId(const IdType &oldId, int facetId);

    void OrderVertices2f(std::vector<Point2f> &vertices,
                         Polygon &orderedPolygon);

    void ProjectParticleToXY(std::vector<Point2f> &projected);

    void RemoveDublicatedVertices2f(const std::vector<Point2f> &projected,
                                  std::vector<Point2f> &cleared);

    void SetOutputPolygon(const Point3f *outputPoints, int outputSize,
                          Polygon &polygon) const;

    void SetConvexOutputPolygon(const Point3f *outputPoints, int outputSize,
                                const Point3f &normal,
                                Polygon &polygon) const;

    bool ProjectToFacetPlane(const Polygon &polygon, const Vector3f &dir,
                             const Point3f &normal, Point3f *projection) const;

};
