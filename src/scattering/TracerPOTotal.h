#pragma once

#include "TracerPO.h"
#include "BeamCache.h"
#include <vector>
#include <string>

class TracerPOTotal : public TracerPO
{
public:
	TracerPOTotal(Particle *particle, int nActs, const std::string &resultFileName);
	void TraceRandom(const AngleRange &betaRange,
					 const AngleRange &gammaRange) override;
    void TraceMonteCarlo(const AngleRange &betaRange,
                         const AngleRange &gammaRange, int nOrientations);
    void TraceFromFile(const std::string &orientFile);

    /// Multi-size computation: trace once, compute diffraction for multiple sizes
    void TraceFromFileMultiSize(const std::string &orientFile,
                                const std::vector<double> &x_sizes,
                                double x_ref);

    /// Sobol quasi-random orientation averaging with particle symmetry
    void TraceFromSobol(int nOrient, double betaSym, double gammaSym);

    /// Adaptive convergence mode
    void TraceAdaptive(double eps, double betaSym, double gammaSym, int maxOrientOverride = 0);
};
