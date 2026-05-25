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
    void TraceEulerQuadratureMultiSize(int nBeta, int nGamma,
                                       double betaSym, double gammaSym,
                                       const std::vector<double> &x_sizes,
                                       double x_ref);
    void TraceLatticeMultiSize(int nOrient, double betaSym, double gammaSym,
                               const std::vector<double> &x_sizes,
                               double x_ref, int generator = 0);

    /// Multi-size with regular beta/gamma grid: trace once at reference size,
    /// cache scale-invariant beam geometry, recompute diffraction for all sizes.
    void TraceRandomMultiSize(const AngleRange &betaRange,
                              const AngleRange &gammaRange,
                              const std::vector<double> &x_sizes,
                              const std::vector<std::string> &labels);

    /// Sobol quasi-random orientation averaging with particle symmetry
    void TraceFromSobol(int nOrient, double betaSym, double gammaSym);
    void TraceFromSobolSeed(int nOrient, unsigned int seed,
                            double betaSym, double gammaSym);
    void TraceFromSobolRing(int nBeta, int nGamma,
                            double betaSym, double gammaSym);
    void TraceFromHammersley(int nOrient, double betaSym, double gammaSym);
    void TraceFromLattice(int nOrient, double betaSym, double gammaSym);
    void TraceFromLatticeGenerator(int nOrient, int generator,
                                   double betaSym, double gammaSym);
    void TraceFromEulerQuadrature(int nBeta, int nGamma,
                                  double betaSym, double gammaSym);
    double TraceFromSobolVariablePhi(int nOrient, double betaSym, double gammaSym,
                                     const std::vector<int> &rowNphi,
                                     const std::vector<int> &rowNorient,
                                     int outputNphi, double fftEps = -1.0,
                                     const std::vector<unsigned int> &owenSeeds =
                                         std::vector<unsigned int>(),
                                     const std::vector<double> &rowBeamCutoff =
                                         std::vector<double>());

    /// Adaptive theta grid: trace once, build theta grid by recursive bisection,
    /// then full diffraction on the found grid.
    /// gridOnly=true: only build theta grid, don't compute full Mueller
    void TraceAdaptiveTheta(int nOrient, double betaSym, double gammaSym,
                             double eps = 0.05, int maxDepth = 8, bool gridOnly = false);

    /// Adaptive convergence mode
    int TraceAdaptive(double eps, double betaSym, double gammaSym,
                      int maxOrientOverride = 0, bool runFinalDiffraction = true,
                      bool strictConvergence = false);

    /// Full 3D sequential optimization: n → N_phi → N_orient
    void TraceAutoFull(double eps, double betaSym, double gammaSym, int maxOrientOverride,
                       Particle *particle, double wave, ScatteringRange &conus,
                       class HandlerPOTotal *handler);

private:
    struct WeightedOrientation
    {
        double beta;
        double gamma;
        double weight;
    };

    void TraceWeightedOrientations(const std::vector<WeightedOrientation> &orientations,
                                   const std::string &label,
                                   double betaSym, double gammaSym);
    void EvaluateOwenControlSample(int nOrient, unsigned int seed,
                                   double betaSym, double gammaSym,
                                   const ScatteringRange &ctrlSphere,
                                   const std::vector<int> &ctrlIdx,
                                   int nAz,
                                   std::vector<double> &ctrlAvg,
                                   double &energyAvg);
    void EvaluateOwenControlBatch(int orientStart, int orientCount,
                                  int sobolCount, unsigned int seed,
                                  double sampleWeight,
                                  double betaSym, double gammaSym,
                                  const ScatteringRange &ctrlSphere,
                                  const std::vector<int> &ctrlIdx,
                                  int nAz,
                                  std::vector<double> &ctrlSum,
                                  double &energySum);

public:
    /// Coherent across orientations (legacy mode, physically incorrect for random)
    bool m_cohOrient = false;

    /// Save intermediate Mueller per beta (--save_betas)
    bool m_saveBetas = false;

    /// Enable checkpoint save/resume for long --orientfile runs.
    bool m_enableCheckpoint = false;

    /// Fast pole handling: compute one gamma at beta poles with gamma weight.
    /// Default false uses honest all-gamma averaging; --pole enables this shortcut.
    bool m_fastPoleGamma = false;

    /// Mirror gamma fundamental domain. The omitted mirror half is restored as
    /// M(phi) -> P M(-phi) P, P=diag(1,1,-1,-1), for full Mueller output.
    bool m_mirrorGamma = false;

    /// Optional manual cap for Sobol streaming chunks. 0 means auto.
    int m_sobolChunkSize = 0;

    /// Points per diffraction ring used by physics-based orientation estimates.
    int m_ringPoints = 3;

    /// Optional final averaging seeds for --autofull nested Owen Sobol.
    std::vector<unsigned int> m_owenAverageSeeds;

    /// Use --autofull search steps, but finish with regular oldauto-style
    /// beta/gamma quadrature instead of Sobol/Owen orientations.
    bool m_oldAutoFullFinal = false;

    /// Last adaptive run's per-theta orientation requirements. 0 means unknown.
    std::vector<int> m_lastAdaptiveRowOrient;

    /// Last adaptive run's per-theta orientation error diagnostics.
    std::vector<double> m_lastAdaptiveRowOrientError;
    std::vector<int> m_lastAdaptiveRowOrientElement;

    /// True when --autofull reached the oldauto-based orientation ceiling and
    /// accepted it as the physics-based maximum instead of chasing noisy tails.
    bool m_lastAdaptiveAcceptedOldautoCap = false;

    /// MPI rank and size (default: single process)
    int m_mpiRank = 0;
    int m_mpiSize = 1;
    void SetMPI(int rank, int size) { m_mpiRank = rank; m_mpiSize = size; }
};
