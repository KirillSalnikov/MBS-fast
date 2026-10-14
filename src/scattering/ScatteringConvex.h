#pragma once

#include "Scattering.h"
#ifdef USE_CUDA
#include "cuda/GpuConvexQueue.h"
#endif

class ScatteringConvex : public Scattering
{
public:
    ScatteringConvex(Particle *particle, Light *incidentLight,
                     bool isOptionalPath, int nActs);

    Scattering* CloneFor(Particle *p, Light *l) override {
        ScatteringConvex *copy = new ScatteringConvex(p, l, true, m_nActs);
        copy->CopyRuntimeOptionsFrom(*this);
        return copy;
    }

    bool ScatterLight(double beta, double gamma, std::vector<Beam> &outBeams) override;
    bool ScatterLight(double, double, const std::vector<std::vector<int>> &,
                      std::vector<Beam> &) override; ///> for predefined trajectories

#ifdef USE_CUDA
    bool PrepareGpuInput(GpuConvexQueue::Input& input, std::vector<Beam>& external);
    void AddGpuStatistics(const GpuConvexQueue::Status& status);
    static Beam DecodeGpuBeam(const GpuConvexQueue::Beam& beam, int diagnosticFacets = 0);
#endif

protected:
    bool PrepareExternalBeams(std::vector<Beam>& outBeams);
    bool TraceInternalBeams(std::vector<Beam> &outBeams);

    bool SplitSecondaryBeams(Beam &incidentBeam, int facetID,
                             Beam &inBeam, std::vector<Beam> &outBeams);
};
