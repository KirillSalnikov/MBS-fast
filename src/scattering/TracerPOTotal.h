#pragma once

#include "TracerPO.h"
#include "BeamCache.h"
#include "AdaptiveConfig.h"
#include <vector>
#include <string>

struct PreparedOrientation;

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

    /// Conservative multi-size regular beta/gamma grid: resize and retrace
    /// every size independently.  This is slower than TraceRandomMultiSize
    /// but avoids reusing reference-size beam topology for non-convex shapes.
    void TraceRandomMultiSizeIndependent(const AngleRange &betaRange,
                                         const AngleRange &gammaRange,
                                         const std::vector<double> &x_sizes,
                                         const std::vector<std::string> &labels,
                                         double x_ref);

    /// Sobol quasi-random orientation averaging with particle symmetry
    void TraceFromSobol(int nOrient, double betaSym, double gammaSym);
    void TraceFromSO3Quaternion(int nOrient);
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
    void TraceFromEulerAdaptiveGamma(int nBeta, int nGammaMax,
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
                            double eps = 0.05, int maxDepth = 8,
                            bool gridOnly = false);

    /// Converge the azimuthal scattering grid while keeping the current theta
    /// grid and reflection depth fixed. Returns the selected N_phi.
    int TraceAdaptivePhi(double eps, int nOrient, double betaSym,
                         double gammaSym);

    /// Converge the maximum internal reflection depth. Returns the selected n.
    int TraceAdaptiveReflections(double eps, int nOrient, double betaSym,
                                 double gammaSym, int startReflection = 0);

    /// Independently refine Gauss beta and periodic gamma counts, then verify
    /// their joint refinement. Returns the selected pair through the outputs.
    void TraceAdaptiveEuler(double eps, double betaSym, double gammaSym,
                            int maxOrientOverride, int &nBeta, int &nGamma,
                            bool runFinalDiffraction = true);

    /// Adaptive convergence mode
    int TraceAdaptive(double eps, double betaSym, double gammaSym,
                      int maxOrientOverride = 0, bool runFinalDiffraction = true,
                      bool strictConvergence = false);

    /// Unified automatic controller used by --auto, --autofull and
    /// --diffraction-autofull.
    void TraceAutoConverged(double eps, double betaSym, double gammaSym,
                            int maxOrientOverride, bool tuneReflections,
                            bool regularEulerFinal, ScatteringRange &conus,
                            bool tunePhi, bool tuneTheta);
    void ResetConvergenceReport(const std::string &mode, double eps);

private:
    struct WeightedOrientation
    {
        double beta;
        double gamma;
        double alpha;
        double qx;
        double qy;
        double qz;
        double qw;
        double weight;
        bool useQuaternion;

        WeightedOrientation()
            : beta(0), gamma(0), alpha(0),
              qx(0), qy(0), qz(0), qw(1),
              weight(1), useQuaternion(false) {}
        WeightedOrientation(double beta_, double gamma_, double weight_)
            : beta(beta_), gamma(gamma_), alpha(0),
              qx(0), qy(0), qz(0), qw(1),
              weight(weight_), useQuaternion(false) {}
        WeightedOrientation(double qx_, double qy_, double qz_, double qw_,
                            double weight_)
            : beta(0), gamma(0), alpha(0),
              qx(qx_), qy(qy_), qz(qz_), qw(qw_),
              weight(weight_), useQuaternion(true) {}
    };

    void TraceWeightedOrientations(
        const std::vector<WeightedOrientation> &orientations,
        const std::string &label, double betaSym, double gammaSym);
    std::vector<WeightedOrientation> BuildSobolOrientations(
        int nOrient, double betaSym, double gammaSym, unsigned int seed) const;
    std::vector<WeightedOrientation> BuildEulerOrientations(
        int nBeta, int nGamma, double betaSym, double gammaSym) const;
    void PrepareOrientationProbe(
        const std::vector<WeightedOrientation> &orientations,
        std::vector<PreparedOrientation> &prepared,
        double &incomingEnergy, double &outputEnergy);
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
    void RecordConvergenceStep(const std::string &parameter, int sweep,
                               int primary, int secondary, double error,
                               double thetaDeg, int element,
                               const std::string &status);
    double m_convergenceTarget = 0.0;

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

    /// Hard safety limits for all adaptive searches. CLI max-* options
    /// override these defaults without changing the requested tolerance.
    AdaptiveConvergenceLimits m_adaptiveLimits;

    /// Optional final averaging seeds for --autofull nested Owen Sobol.
    std::vector<unsigned int> m_owenAverageSeeds;

    /// Last adaptive run's per-theta orientation requirements. 0 means unknown.
    std::vector<int> m_lastAdaptiveRowOrient;

    /// Last adaptive run's per-theta orientation error diagnostics.
    std::vector<double> m_lastAdaptiveRowOrientError;
    std::vector<int> m_lastAdaptiveRowOrientElement;

    /// True when --autofull reached the oldauto-based orientation ceiling and
    /// accepted it as the physics-based maximum instead of chasing noisy tails.
    bool m_lastAdaptiveAcceptedOldautoCap = false;

    /// Last --oldautofull final regular-grid parameters.  Used to reuse the
    /// strict grid for shared multikeq/multigrid after tuning on the largest size.
    bool m_lastOldAutoFullGridValid = false;
    int m_lastOldAutoFullN = 0;
    int m_lastOldAutoFullNphi = 0;
    int m_lastOldAutoFullNBeta = 0;
    int m_lastOldAutoFullNGamma = 0;
    int m_lastOldAutoFullDiv = 0;
    double m_lastOldAutoFullBetaSym = 0.0;
    double m_lastOldAutoFullGammaSym = 0.0;

    /// MPI rank and size (default: single process)
    int m_mpiRank = 0;
    int m_mpiSize = 1;
    void SetMPI(int rank, int size) { m_mpiRank = rank; m_mpiSize = size; }
};
