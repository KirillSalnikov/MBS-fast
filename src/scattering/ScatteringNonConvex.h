#pragma once

#include "Scattering.h"

#ifdef USE_CUDA
struct GpuTraceExactHit;
#endif

/** NOTE: пучки выходят со случайно ориентированным порядком вершин */
class ScatteringNonConvex : public Scattering
{
public:
    ScatteringNonConvex(Particle *particle, Light *incidentLight,
                        bool isOpticalPath, int nActs);

    Scattering* CloneFor(Particle *p, Light *l) override {
        ScatteringNonConvex *copy = new ScatteringNonConvex(p, l, true, m_nActs);
        copy->CopyRuntimeOptionsFrom(*this);
        return copy;
    }
    void PrepareForParallelTrace() override;

    bool ScatterLight(double beta, double gamma, std::vector<Beam> &scaterredBeams) override;
    bool ScatterLight(double beta, double gamma, const std::vector<std::vector<int>> &tracks,
                             std::vector<Beam> &scaterredBeams) override;
private:
    void SortFacets_faster(const Point3f &beamDir, IntArray &facetIDs);
    bool CutBeamByFacet(const Polygon &intersection, int facetId, Beam &beam,
                        Polygon &reachedIntersection, PolygonArray &result);

    void CutExternalBeam(const Beam &beam, std::vector<Beam> &scaterredBeams,
                         const IntArray *readyFacetIds = nullptr);

    void FindVisibleFacets(const Beam &beam, IntArray &facetIds);
    void FindVisibleFacetsForLight(IntArray &facetIDs);
    void BuildFacetVisibilityCache();
    void BuildFacetNormalCache();
    void BuildFacetProjectionCache();
    void CopyVisibilityCacheFrom(const ScatteringNonConvex &source);

    void SelectVisibleFacets(const Beam &beam, IntArray &facetIDs);
    void SortVisibleFacets(const Beam &beam, IntArray &facetIDs);
    void SelectVisibleFacetsForLight(IntArray &facetIDs);

    bool SetOpticalBeamParams(const Facet &facet, const Beam &incidentBeam,
                              Beam &inBeam, Beam &outBeam);

    void IntersectWithFacet(const IntArray &facetIds, int prevFacetNum,
                            PolygonArray &resFacets);

    bool SplitLightToBeams();

    bool IsOutgoingBeam(Beam &incidentBeam);

    void TraceFirstBeamFixedFacet(int facetID, bool &isIncident);

    bool PushBeamsToTree(int facetId, const PolygonArray &polygons,
                         Beam &inBeam, Beam &outBeam);

    bool IsVisibleFacet(int facetID, const Beam &beam);
    bool MayBeamIntersectFacetProjected(const Beam &beam, int facetId) const;

    bool SplitByFacet(const IntArray &facetIDs, int facetIndex);

    bool SplitBeamByFacet(const Polygon &intersection, int facetId,
                          Beam &beam, bool &ok
#ifdef USE_CUDA
                          , const GpuTraceExactHit *preparedHit = nullptr
#endif
                          );

    void CutPolygonByFacets(const Polygon &pol,
                            const IntArray &facetIds, size_t size,
                            const Vector3f &polNormal, const Vector3f &clipNormal,
                            const Vector3f &dir,
                            PolygonArray &pols);
    void CutPolygonByAggregateParts(const Polygon &pol, int facetId,
                                    const Vector3f &polNormal,
                                    const Vector3f &dir,
                                    PolygonArray &pols);

    bool PushBeamPartsToTree(const Beam &beam,
                             const PolygonArray &parts);
    template<class T>
    bool PushBeamToTree(Beam &beam, const Beam &oldBeam,
                        const T &newId, int facetId,
                        Location loc);

protected:
    bool SplitBeams(std::vector<Beam> &scaterredBeams);

private:
    PolygonArray m_polygonBuffer;
    PolygonArray m_polygonResultBuffer;
    int m_visibleFacetCache[2][MAX_FACET_NUM][MAX_FACET_NUM];
    size_t m_visibleFacetCacheSize[2][MAX_FACET_NUM];
    double m_facetNormalComponents[2][MAX_FACET_NUM][3];
    double m_facetProjectionBounds[3][MAX_FACET_NUM][4];
    int m_facetProjectionDrop[2][MAX_FACET_NUM];
    size_t m_traceSortCalls = 0;
    size_t m_traceSortCandidates = 0;
    size_t m_traceSortOverlapCalls = 0;
    const ScatteringNonConvex *m_visibilityCacheOwner = nullptr;
    bool m_visibilityCacheBuilt = false;
};
