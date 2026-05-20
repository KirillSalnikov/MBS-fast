#pragma once

#include "Scattering.h"

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

protected:
    void TraceInternalBeams(std::vector<Beam> &outBeams);

    bool SplitSecondaryBeams(Beam &incidentBeam, int facetID,
                             Beam &inBeam, std::vector<Beam> &outBeams);
};
