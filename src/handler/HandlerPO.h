#pragma once

#include "Handler.h"
#include "BeamCache.h"
#include <vector>

/// Preprocessed beam data for parallel direction-loop processing.
/// Contains all scalar data extracted from a beam after sequential preprocessing
/// (RotateSpherical, ComputeBeamInfo, PrecomputeEdgeData, etc.).
/// The direction loop can run on this without touching any Handler member state.
struct PreparedBeam
{
    BeamEdgeData edgeData;
    BeamPolData  polData;
    BeamInfo     info;

    // Beam direction (double precision copy)
    double bdx, bdy, bdz;
    // Aperture axes
    double horAx, horAy, horAz;
    double verAx, verAy, verAz;
    double normx, normy, normz;
    // Center
    double cenx, ceny, cenz;
    double beam_area;
    // PolData scalars (redundant with polData, but pre-extracted for hot loop)
    double pNTx, pNTy, pNTz;
    double pNPx, pNPy, pNPz;
    double pnxDTx, pnxDTy, pnxDTz;
    double pnxDPx, pnxDPy, pnxDPz;
    // J_phased elements
    double jp00r, jp00i, jp01r, jp01i;
    double jp10r, jp10i, jp11r, jp11i;
    bool isExternal;

    // Original beam data needed for fallback path (non-valid edgeData)
    Beam   origBeam;

    // Internal optical path samples used to reapply absorption when a
    // prepared reference-size beam is rescaled in multikeq/multigrid.
    std::vector<double> absorptionPaths;
    double outputCrossSection = 0.0;
    double outputMueller00 = 0.0;
};

/// All preprocessed beams from one orientation, ready for parallel processing.
struct PreparedOrientation
{
    std::vector<PreparedBeam> beams;
    double sinZenith;  // weight for this orientation
    double extinctionOt = 0.0;
};

class HandlerPO : public Handler
{
public:
    HandlerPO(Particle *particle, Light *incidentLight, int nTheta,
              double wavelength);

    void HandleBeams(std::vector<Beam> &beams, double sinZenith) override;

    // --- Multi-size beam caching API ---
    /// Cache beams from one orientation into the BeamCache.
    /// D_ref is the reference particle diameter used during this tracing.
    void CacheBeams(std::vector<Beam> &beams, double weight,
                    double D_ref, double incomingEnergy,
                    OrientationBeams &out);

    /// Compute Mueller matrices for multiple size parameters from cached beams.
    /// x_sizes: vector of size parameters (x = pi*D/lambda)
    /// results_M: output Mueller matrices, one Arr2D per size
    void ComputeFromCache(const BeamCache &cache,
                          const std::vector<double> &x_sizes,
                          std::vector<Arr2D> &results_M,
                          std::vector<double> &results_energy);
    /// Preprocess beams from one orientation into PreparedOrientation.
    /// Must be called sequentially (uses m_isBadBeam, modifies beams).
    void PrepareBeams(std::vector<Beam> &beams, double sinZenith,
                      PreparedOrientation &out);
    double ComputeForwardExtinctionOt(const PreparedOrientation &prepared) const;
    double ComputeForwardExtinctionOtScaled(const PreparedOrientation &prepared,
                                            double scale,
                                            double waveIndex,
                                            double absorptionCoefficient) const;

    /// Copy immutable settings needed by PrepareBeams into a worker-local
    /// handler. The worker gets its own Scattering/Particle state.
    void ConfigureForThreadLocalPrepare(const HandlerPO &source,
                                        Scattering *scattering);

    /// Process prepared beams into a LOCAL Mueller accumulator.
    /// Thread-safe: reads only from handler's immutable data (sphere, wave constants)
    /// and writes only to the provided localM.
    void HandleBeamsToLocal(const PreparedOrientation &prepared,
                            Arr2D &localM,
                            std::vector<Arr2DC> &localJ,
                            std::vector<Arr2DC> *localJ_noshadow = nullptr);
    bool HandleBeamsToLocalGpu(const PreparedOrientation &prepared,
                               Arr2D &localM,
                               Arr2D &localM_noshadow);
    bool HandleOrientationsToLocalGpu(const std::vector<PreparedOrientation> &prepared,
                                      Arr2D &localM,
                                      Arr2D &localM_noshadow);
    bool HandleOrientationsToLocalGpu(const std::vector<PreparedOrientation> &prepared,
                                      int start,
                                      int count,
                                      Arr2D &localM,
                                      Arr2D &localM_noshadow,
                                      double scale = 1.0,
                                      double waveIndex = 0.0);
    bool HandleOrientationsToLocalGpuFftPhi(const std::vector<PreparedOrientation> &prepared,
                                            int start,
                                            int count,
                                            Arr2D &localM,
                                            Arr2D &localM_noshadow,
                                            double scale = 1.0,
                                            double waveIndex = 0.0);
    bool HandleOrientationsToLocalGpuMultiK(const std::vector<PreparedOrientation> &prepared,
                                            int start,
                                            int count,
                                            const std::vector<double> &scales,
                                            double waveIndex,
                                            std::vector<Arr2D> &localMs);
    bool HandleOrientationsToLocalGpuFftPhiMultiK(const std::vector<PreparedOrientation> &prepared,
                                                  int start,
                                                  int count,
                                                  const std::vector<double> &scales,
                                                  double waveIndex,
                                                  std::vector<Arr2D> &localMs);
    int SelectGpuOrientationBatchSize(const std::vector<PreparedOrientation> &prepared,
                                      int start,
                                      int maxCount) const;

    /// Fast diffraction for control points only (4 theta indices, phi=0).
    /// Returns M11 at each control angle. ~40000× faster than full grid.
    void DiffractControlPoints(const PreparedOrientation &prepared,
                                const int *thetaIndices, int nPoints,
                                double *m11_out);

    /// Diffract at arbitrary theta values (radians), phi-averaged.
    /// Like DiffractControlPoints but not tied to existing theta grid.
    void DiffractAtThetas(const PreparedOrientation &prepared,
                           const double *theta_rads, int nPoints,
                           double *m11_out);
    bool DiffractThetasGpu(const std::vector<PreparedOrientation> &prepared,
                           const double *theta_rads,
                           int nPoints,
                           std::vector<double> &m11_out);

    /// Convert coherent Jones (localJ) to Mueller and add to localM.
    static void AddToMuellerLocal(const std::vector<Arr2DC> &localJ,
                                  double normIndex, Arr2D &localM,
                                  int nAz, int nZen);

    void WriteMatricesToFile(std::string &destName, double nrg) override;
    void WriteTotalMatricesToFile(const std::string &destName) override;
    void WriteJonesToFile(const std::string &destName);
    // double ComputeTotalScatteringEnergy() override;

    void SetScatteringSphere(const ScatteringRange &grid) override;

    void SetBackScatteringConus(double radAngle);
    void SetGpuEnabled(bool value);
    bool IsGpuEnabled() const;
    void SetFftEnabled(bool value);
    bool IsFftEnabled() const;
    void SetFftPhiFactor(int value);
    int FftPhiFactor() const;
    void AutoSelectFftPhiFactor(double eps);
    static bool HasNumericFftPhiFactorOverride();
    static int SelectAutoFftPhiFactor(int nPhi, double eps);
    void SetFullOnly(bool value);
    bool IsFullOnly() const;
    bool ComputeNoShadow() const;
    bool HasAbsorptionAccounting() const;

    matrix *m_Lp;
    matrix *m_Ln;
    Arr2D M;				// Mueller matrices (all beams including shadow)
    Arr2D M_noshadow;		// Mueller matrices (without shadow/external beam)

protected:
    virtual void AddToMueller();

    static bool IsParticleBeam(const Beam &beam);
    static bool HasInternalOpticalPath(const Beam &beam);

    void ComputeOpticalLengths(const Beam &beam, BeamInfo &info);

    virtual void RotateJones(const Beam &beam, const BeamInfo &info,
                     const Vector3d &vf, const Vector3d &direction,
                     matrixC &matrix) const;
    static void PrecomputePolData(const Beam &beam, const BeamInfo &info,
                                  BeamPolData &polData);
    static void RotateJonesFast(const BeamPolData &polData,
                                const Vector3d &vf, const Vector3d &direction,
                                matrixC &matrix);
public:
    void CleanJ();
protected:
    matrixC ComputeFnJones(const Matrix2x2c &matrix, const BeamInfo &info,
                           const Vector3d &direction);

    matrixC ApplyDiffraction(const Beam &beam, const BeamInfo &info,
                         const Vector3d &direction, const Vector3d &vf,
                         bool useAbsorptionIntegral = true);
    matrixC ApplyDiffractionFast(const Beam &beam, const BeamInfo &info,
                                 const BeamEdgeData &edgeData,
                                 const Point3d &beamDirD,
                                 const Vector3d &direction, const Vector3d &vf);
    matrixC ApplyDiffractionFast2(const Beam &beam, const BeamInfo &info,
                                  const BeamEdgeData &edgeData,
                                  const Point3d &beamDirD,
                                  const matrixC &J_phased,
                                  bool isExternal,
                                  const Vector3d &direction, const Vector3d &vf);
    matrixC ApplyDiffractionFast3(const BeamPolData &polData,
                                  const BeamInfo &info,
                                  const BeamEdgeData &edgeData,
                                  const Point3d &beamDirD,
                                  const matrixC &J_phased,
                                  bool isExternal,
                                  const Vector3d &direction, const Vector3d &vf);

    BeamInfo ComputeBeamInfo(Beam &beam);



protected:
    std::vector<Arr2D> m_groupMatrices;	//
    std::vector<Arr2DC> m_diffractedMatrices;	// Jones matrices
public:
    bool outputJones = false;
    bool m_gpuEnabled = false;
    bool m_fftEnabled = false;
    int m_fftPhiFactor = 0;
    bool m_fullOnly = true;
protected:
    bool isNanOccured = false;
    bool isNan = false;
    bool isBackScatteringConusEnabled = false;
    double backScatteringConus = 180;
public:
    bool useKarczewski = false;
protected:
    void KarczewskiJones(const Beam &beam, const BeamInfo &info,
                         const Vector3d &vf, const Vector3d &direction,
                         matrixC &matrix) const;

    // Handler interface
public:
    void SetTracks(Tracks *tracks) override;
private:
    void WriteGroupMatrices(Arr2D &matrices, const std::string &name);
};
