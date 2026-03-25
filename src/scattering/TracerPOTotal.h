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

    /// Multi-size with Sobol: trace reference size, cache, recompute for all x_sizes
    void TraceSobolMultiSize(int nOrient, double betaSym, double gammaSym,
                              const std::vector<double> &x_sizes, double x_ref);

    /// Sobol quasi-random orientation averaging with particle symmetry
    void TraceFromSobol(int nOrient, double betaSym, double gammaSym);

    /// Adaptive convergence mode
    void TraceAdaptive(double eps, double betaSym, double gammaSym, int maxOrientOverride = 0);

    /// Full 3D sequential optimization: n → N_phi → N_orient
    void TraceAutoFull(double eps, double betaSym, double gammaSym, int maxOrientOverride,
                       Particle *particle, double wave, ScatteringRange &conus,
                       class HandlerPOTotal *handler);

    /// Coherent across orientations (legacy mode, physically incorrect for random)
    bool m_cohOrient = false;

    /// MPI rank and size (default: single process)
    int m_mpiRank = 0;
    int m_mpiSize = 1;
    void SetMPI(int rank, int size) { m_mpiRank = rank; m_mpiSize = size; }
};
