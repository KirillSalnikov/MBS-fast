#pragma once

#include "Handler.h"
#include "BeamCache.h"
#include "CutoffStatistics.h"
#include "FftStatistics.h"
#include <memory>
#include <utility>
#include <vector>

struct ScatteringBasisSoA
{
    std::vector<double> vfx;
    std::vector<double> vfy;
    std::vector<double> vfz;
    std::vector<double> vtx;
    std::vector<double> vty;
    std::vector<double> vtz;

    explicit ScatteringBasisSoA(size_t size)
        : vfx(size), vfy(size), vfz(size),
          vtx(size), vty(size), vtz(size)
    {
    }
};

/// Edge data stored with the prepared beam. Most traced polygons have at most
/// eight vertices, so keeping 32 entries inline wastes both RAM and cache.
/// Larger valid polygons retain the full representation transparently.
struct PreparedEdgeData
{
    static const int INLINE_EDGES = 8;

    PreparedEdgeData()
    {
        BindInline();
    }

    PreparedEdgeData(const PreparedEdgeData &other)
    {
        CopyFrom(other);
    }

    PreparedEdgeData &operator=(const PreparedEdgeData &other)
    {
        if (this != &other)
            CopyFrom(other);
        return *this;
    }

    PreparedEdgeData(PreparedEdgeData &&other) noexcept
    {
        MoveFrom(std::move(other));
    }

    PreparedEdgeData &operator=(PreparedEdgeData &&other) noexcept
    {
        if (this != &other)
            MoveFrom(std::move(other));
        return *this;
    }

    void Assign(const BeamEdgeData &source)
    {
        nVertices = source.nVertices;
        valid = source.valid;
        overflow.reset();
        if (!valid)
        {
            nVertices = 0;
            BindInline();
            return;
        }
        if (valid && nVertices > INLINE_EDGES)
        {
            overflow = std::make_shared<BeamEdgeData>(source);
            BindOverflow();
            return;
        }

        BindInline();
        for (int i = 0; i < nVertices; ++i)
        {
            x[i] = source.x[i];
            y[i] = source.y[i];
            slope_yx[i] = source.slope_yx[i];
            slope_xy[i] = source.slope_xy[i];
            edge_valid_x[i] = source.edge_valid_x[i];
            edge_valid_y[i] = source.edge_valid_y[i];
        }
    }

    bool HasOverflow() const
    {
        return static_cast<bool>(overflow);
    }

    double *x;
    double *y;
    double *slope_yx;
    double *slope_xy;
    bool *edge_valid_x;
    bool *edge_valid_y;
    int nVertices = 0;
    bool valid = false;

private:
    void BindInline()
    {
        x = inlineX;
        y = inlineY;
        slope_yx = inlineSlopeYx;
        slope_xy = inlineSlopeXy;
        edge_valid_x = inlineValidX;
        edge_valid_y = inlineValidY;
    }

    void BindOverflow()
    {
        x = overflow->x;
        y = overflow->y;
        slope_yx = overflow->slope_yx;
        slope_xy = overflow->slope_xy;
        edge_valid_x = overflow->edge_valid_x;
        edge_valid_y = overflow->edge_valid_y;
    }

    void CopyFrom(const PreparedEdgeData &other)
    {
        nVertices = other.nVertices;
        valid = other.valid;
        if (other.overflow)
        {
            overflow = std::make_shared<BeamEdgeData>(*other.overflow);
            BindOverflow();
            return;
        }

        overflow.reset();
        BindInline();
        for (int i = 0; i < nVertices; ++i)
        {
            x[i] = other.x[i];
            y[i] = other.y[i];
            slope_yx[i] = other.slope_yx[i];
            slope_xy[i] = other.slope_xy[i];
            edge_valid_x[i] = other.edge_valid_x[i];
            edge_valid_y[i] = other.edge_valid_y[i];
        }
    }

    void ResetMovedFrom(PreparedEdgeData &other)
    {
        other.nVertices = 0;
        other.valid = false;
        other.overflow.reset();
        other.BindInline();
    }

    void MoveFrom(PreparedEdgeData &&other)
    {
        nVertices = other.nVertices;
        valid = other.valid;
        if (other.overflow)
        {
            overflow = std::move(other.overflow);
            BindOverflow();
            ResetMovedFrom(other);
            return;
        }

        overflow.reset();
        BindInline();
        for (int i = 0; i < nVertices; ++i)
        {
            x[i] = other.x[i];
            y[i] = other.y[i];
            slope_yx[i] = other.slope_yx[i];
            slope_xy[i] = other.slope_xy[i];
            edge_valid_x[i] = other.edge_valid_x[i];
            edge_valid_y[i] = other.edge_valid_y[i];
        }
        ResetMovedFrom(other);
    }

    double inlineX[INLINE_EDGES];
    double inlineY[INLINE_EDGES];
    double inlineSlopeYx[INLINE_EDGES];
    double inlineSlopeXy[INLINE_EDGES];
    bool inlineValidX[INLINE_EDGES];
    bool inlineValidY[INLINE_EDGES];
    std::shared_ptr<BeamEdgeData> overflow;
};

struct PreparedBeamFallback
{
    PreparedBeamFallback(const Beam &sourceBeam, const BeamInfo &sourceInfo)
        : beam(sourceBeam), info(sourceInfo)
    {
    }

    Beam beam;
    BeamInfo info;
};

/// Preprocessed beam data for parallel direction-loop processing.
/// Contains all scalar data extracted from a beam after sequential preprocessing
/// (RotateSpherical, ComputeBeamInfo, PrecomputeEdgeData, etc.).
/// The direction loop can run on this without touching any Handler member state.
struct PreparedBeam
{
    // The preprocessing path overwrites every field it consumes.  A
    // user-provided constructor prevents value-initialization from clearing
    // the large fixed edge arrays before they are immediately filled.
    PreparedBeam() {}

    PreparedEdgeData edgeData;

    // Beam direction (double precision copy)
    double bdx, bdy, bdz;
    // Aperture axes
    double horAx, horAy, horAz;
    double verAx, verAy, verAz;
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

    // Compact source data used by scaling and CUDA packing.  Keeping a full
    // Beam here costs more than 2 KiB per output beam, mostly unused polygon
    // storage.  Only invalid edge data needs the legacy diffraction fallback.
    Matrix2x2c rawJ;
    double opticalPath;
    double projLength;
    int nActs;
    std::shared_ptr<PreparedBeamFallback> fallback;

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
    ~HandlerPO() override;

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
    bool HandleOrientationsToLocalGpuOnDevice(
                                      const std::vector<PreparedOrientation> &prepared,
                                      int start,
                                      int count,
                                      int device,
                                      Arr2D &localM,
                                      Arr2D &localM_noshadow,
                                      double scale = 1.0,
                                      double waveIndex = 0.0,
                                      unsigned long long cacheToken = 0);
    bool HandleOrientationsToLocalGpuCached(
                                      const std::vector<PreparedOrientation> &prepared,
                                      int start,
                                      int count,
                                      Arr2D &localM,
                                      Arr2D &localM_noshadow,
                                      double scale,
                                      double waveIndex,
                                      unsigned long long cacheToken);
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
    void SetFftTolerance(double value);
    double FftTolerance() const;
    std::string FftReport() const;
    void AutoSelectFftPhiFactor(double eps);
    static bool HasNumericFftPhiFactorOverride();
    static int SelectAutoFftPhiFactor(int nPhi, double eps);
    void SetFullOnly(bool value);
    bool IsFullOnly() const;
    bool ComputeNoShadow() const;
    bool HasAbsorptionAccounting() const;
    std::string BeamCutoffReport() const;

    matrix *m_Lp;
    matrix *m_Ln;
    Arr2D M;				// Mueller matrices (all beams including shadow)
    Arr2D M_noshadow;		// Mueller matrices (without shadow/external beam)

protected:
    virtual void AddToMueller();

    static bool IsParticleBeam(const Beam &beam);
    static bool HasInternalOpticalPath(const Beam &beam);

    void ComputeOpticalLengths(const Beam &beam, BeamInfo &info,
                               const std::vector<int> *track = nullptr);

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

    BeamInfo ComputeBeamInfo(Beam &beam,
                             const std::vector<int> *track = nullptr);



protected:
    std::vector<Arr2D> m_groupMatrices;	//
    std::vector<Arr2DC> m_diffractedMatrices;	// Jones matrices
public:
    bool outputJones = false;
    bool m_gpuEnabled = false;
    bool m_fftEnabled = false;
    int m_fftPhiFactor = 0;
    double m_fftTolerance = 0.0;
    bool m_fullOnly = true;
    std::string m_cutoffProfileName = "legacy-default";
protected:
    std::shared_ptr<BeamCutoffStatistics> m_beamCutoffStatistics;
    std::shared_ptr<FftInterpolationStatistics> m_fftStatistics;
    std::shared_ptr<const std::vector<Point3d>> m_transverseBasis;
    std::shared_ptr<const ScatteringBasisSoA> m_scatteringBasisSoA;
    int m_transverseThetaStride = 0;
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
