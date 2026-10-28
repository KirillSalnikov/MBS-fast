#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, Light *incidentLight,
                                   bool isOpticalPath, int nActs)
    : Scattering(particle, incidentLight, isOpticalPath, nActs)
{
}

bool ScatteringConvex::ScatterLight(double /*beta*/, double /*gamma*/,
                                    std::vector<Beam> &outBeams)
{
    return PrepareExternalBeams(outBeams) && TraceInternalBeams(outBeams);
}

bool ScatteringConvex::PrepareExternalBeams(std::vector<Beam>& outBeams)
{
    // m_particle->Rotate(beta, gamma, 0);

    m_incidentEnergy = 0;
    m_treeSize = 0;
    m_beamTree.clear();
    ResetTraceReference();

    /// first extermal beam
    for (int facetID = 0; facetID < m_particle->nFacets; ++facetID)
    {
        const Point3f &inNormal = m_facets[facetID].in_normal;
        splitting.ComputeCosA(m_incidentDir, inNormal);

        if (!splitting.IsIncident()) /// beam is not incident to this facet
        {
            continue;
        }

        Beam inBeam, outBeam;
        SplitLightToBeams(facetID, inBeam, outBeam);

        auto newId = RecomputeTrackId(0, facetID);

        outBeam.id = newId;
        outBeam.lastFacetId = facetID;
        outBeam.nActs = 0;
        UpdateTraceReference(outBeam);
        UpdateTraceReference(inBeam);
        outBeams.push_back(outBeam);

        inBeam.id = newId;
        if (!PushBeamToTree(inBeam, facetID, 0, Location::In))
            return false;

        AddProjectedIncidentEnergy(facetID, outBeam);
    }

    return true;
}

#ifdef USE_CUDA
namespace {
GpuConvexQueue::Vec PackVec(const Point3f& p) { return {p.cx,p.cy,p.cz}; }
Point3f UnpackVec(const GpuConvexQueue::Vec& p) { return Point3f(p.x,p.y,p.z); }
GpuConvexQueue::Beam PackBeam(const Beam& b) {
    GpuConvexQueue::Beam p={};p.count=b.nVertices;p.depth=b.lastFacetId==__INT_MAX__?0:b.nActs;p.facet=b.lastFacetId;
    p.track[0]=b.lastFacetId;
    for(int v=0;v<p.count;++v)p.vertices[v]=PackVec(b.arr[v]);
    p.direction=PackVec(b.direction);p.basis=PackVec(b.polarizationBasis);p.path=b.opticalPath;p.front=b.lastFacetId==__INT_MAX__?0:b.front;
    const complex j[]={b.J.m11,b.J.m12,b.J.m21,b.J.m22};
    for(int i=0;i<4;++i)p.j[i]={real(j[i]),imag(j[i])};return p;
}
}
bool ScatteringConvex::EncodeGpuBeam(const Beam& beam,GpuConvexQueue::Beam& output) {
    if(beam.nVertices<3||beam.nVertices>GpuConvexQueue::MaxVertices)return false;
    // Shadow beams do not initialize the internal ray's front/depth fields.
    output=PackBeam(beam);return true;
}
bool ScatteringConvex::PrepareGpuInput(GpuConvexQueue::Input& input,std::vector<Beam>& external) {
    // Track/absorption reconstruction and nonconvex clipping are deliberately
    // excluded from this first prototype. Unsupported cases keep the CPU path.
    const complex ri=m_particle->GetRefractiveIndex();
    if(m_trackIdsRequired||imag(ri)!=0||!(real(ri)>0)||m_traceMaxBeams!=0
       ||m_nActs<0||m_nActs>18||m_particle->nFacets>GpuConvexQueue::MaxFacets)return false;
    for(int f=0;f<m_particle->nFacets;++f)
        if(m_facets[f].nVertices<3||m_facets[f].nVertices>GpuConvexQueue::MaxVertices)return false;
    if(!PrepareExternalBeams(external))return false;
    for(const auto& b:m_beamTree)if(b.id!=0)return false;
    input.faceCount=m_particle->nFacets;input.seedCount=m_beamTree.size();input.maxDepth=m_nActs;input.ri=real(ri);
    input.lengthTolerance=m_geometryLengthTolerance;input.areaTolerance=m_geometryAreaTolerance;input.depthTolerance=m_geometryDepthTolerance;
    input.reference[0]=m_traceRefJNorm;input.reference[1]=m_traceRefArea;input.reference[2]=m_traceRefImportance;
    input.cutoff[0]=m_traceCutoffJRel;input.cutoff[1]=m_traceCutoffAreaRel;input.cutoff[2]=m_traceCutoffImportanceRel;
    for(int f=0;f<input.faceCount;++f){auto& out=input.faces[f];const Facet& in=m_facets[f];out.count=in.nVertices;
        out.inward=PackVec(in.in_normal);out.outward=PackVec(in.ex_normal);out.plane=in.in_normal.d_param;
        for(int v=0;v<in.nVertices;++v)out.vertices[v]=PackVec(in.arr[v]);}
    for(int i=0;i<input.seedCount;++i)input.seeds[i]=PackBeam(m_beamTree[i]);
    return true;
}
Beam ScatteringConvex::DecodeGpuBeam(const GpuConvexQueue::Beam& p,int diagnosticFacets) {
    Beam b;b.nVertices=p.count;for(int v=0;v<p.count;++v)b.arr[v]=UnpackVec(p.vertices[v]);
    b.direction=UnpackVec(p.direction);b.polarizationBasis=UnpackVec(p.basis);b.opticalPath=p.path;b.front=p.front;
    b.J.m11=complex(p.j[0].re,p.j[0].im);b.J.m12=complex(p.j[1].re,p.j[1].im);
    b.J.m21=complex(p.j[2].re,p.j[2].im);b.J.m22=complex(p.j[3].re,p.j[3].im);
    b.id=0;if(diagnosticFacets)for(int i=0;i<=p.depth;++i)
        b.id.multiplyAddNonnegative(diagnosticFacets+1,(p.track[i]+1)*(diagnosticFacets+1));
    b.locations=0;b.lastFacetId=p.facet;b.nActs=p.depth;b.location=Location::Out;return b;
}
void ScatteringConvex::AddGpuStatistics(const GpuConvexQueue::Status& s) {
    m_traceCutoffStatistics->evaluated.fetch_add(s.evaluated,std::memory_order_relaxed);
    m_traceCutoffStatistics->rejected.fetch_add(s.rejected,std::memory_order_relaxed);
    m_traceCutoffStatistics->rejectedJones.fetch_add(s.rejectedJones,std::memory_order_relaxed);
    m_traceCutoffStatistics->rejectedArea.fetch_add(s.rejectedArea,std::memory_order_relaxed);
    m_traceCutoffStatistics->rejectedImportance.fetch_add(s.rejectedImportance,std::memory_order_relaxed);
}
#endif

bool ScatteringConvex::ScatterLight(double, double, const std::vector<std::vector<int>> &/*tracks*/, std::vector<Beam> &)
{
    return true;
}

bool ScatteringConvex::TraceInternalBeams(std::vector<Beam> &outBeams)
{
    while (m_treeSize != 0)
    {
        Beam beam = m_beamTree.back();
        m_beamTree.pop_back();
        m_treeSize = (int)m_beamTree.size();

        if (IsTerminalAct(beam))
        {
            continue;
        }

        // A rejected facet does not produce a child. Reuse the scratch child
        // instead of constructing its Jones matrix and track id for every
        // candidate facet. SplitSecondaryBeams overwrites it on every hit.
        Beam inBeam;
        for (int id = 0; id < m_particle->nFacets; ++id)
        {
            if (id == beam.lastFacetId)
            {
                continue;
            }

            bool isIncident = SplitSecondaryBeams(beam, id, inBeam, outBeams);

            if (!isIncident)
            {
                continue;
            }

            inBeam.id = RecomputeTrackId(beam.id, id);
            inBeam.locations = beam.locations;
            if (!IsTracePruned(inBeam)
                && !PushBeamToTree(inBeam, id, beam.nActs+1, Location::In))
                return false;
        }
    }

    return true;
}

bool ScatteringConvex::SplitSecondaryBeams(Beam &incidentBeam, int facetID,
                                           Beam &inBeam, std::vector<Beam> &outBeams)
{
    const Point3f &incidentDir = incidentBeam.direction;

    // ext. normal uses in this calculating
    const Point3f &normal = m_facets[facetID].ex_normal;
    splitting.ComputeCosA(normal, incidentDir);

    if (!splitting.IsIncident())
    {
        return false;
    }

    // Most candidate facets miss the aperture. Clipping only needs a polygon;
    // defer the optical state (Jones matrix, track id, etc.) until a hit.
    Polygon aperture;
    Intersect(facetID, incidentBeam, aperture);

    if (aperture.nVertices < MIN_VERTEX_NUM)
    {
        return false;
    }

    Beam outBeam;
    outBeam = aperture;

#ifdef _DEBUG // DEB
    auto id = RecomputeTrackId(incidentBeam.id, facetID);
    if (id == 423)
        int fff = 0;
#endif
    inBeam = outBeam;
    splitting.ComputeSplittingParams(incidentBeam.direction, normal);
//    ComputePolarisationParams(incidentBeam.direction, normal, incidentBeam);

    if (!splitting.IsNormalIncidence())
    {	// regular incidence
        Beam incBeam = incidentBeam;

        ComputePolarisationParams(incidentBeam.direction, normal, incBeam);

        if (!splitting.IsCompleteReflection())
        {
            outBeam.id = RecomputeTrackId(incidentBeam.id, facetID);

            splitting.ComputeRegularBeamsParams(normal, incBeam, inBeam, outBeam);
            outBeam.nActs = incidentBeam.nActs + 1;
            outBeam.opticalPath += splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
            outBeam.lastFacetId = facetID;
            outBeams.push_back(outBeam);
        }
        else // complete internal reflection incidence
        {
            splitting.ComputeCRBeamParams(normal, incBeam, inBeam);
        }
    }
    else
    {	// normal incidence
        splitting.ComputeNormalBeamParams(incidentBeam, inBeam, outBeam);

        outBeam.nActs = incidentBeam.nActs + 1;
        outBeam.id = RecomputeTrackId(incidentBeam.id, facetID);
        double path = splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
        outBeam.opticalPath += path;
        outBeam.lastFacetId = facetID;
        outBeams.push_back(outBeam);
    }

    return true;
}
