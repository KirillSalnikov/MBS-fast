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
};

/// All preprocessed beams from one orientation, ready for parallel processing.
struct PreparedOrientation
{
    std::vector<PreparedBeam> beams;
    double sinZenith;  // weight for this orientation
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

    /// Process prepared beams into a LOCAL Mueller accumulator.
    /// Thread-safe: reads only from handler's immutable data (sphere, wave constants)
    /// and writes only to the provided localM.
    void HandleBeamsToLocal(const PreparedOrientation &prepared,
                            Arr2D &localM,
                            std::vector<Arr2DC> &localJ);

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

    matrix *m_Lp;
    matrix *m_Ln;
    Arr2D M;				// Mueller matrices

protected:
    virtual void AddToMueller();

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
                         const Vector3d &direction, const Vector3d &vf);
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
protected:
    bool isNanOccured = false;
    bool isNan = false;
    bool isBackScatteringConusEnabled = false;
    double backScatteringConus = 180;
public:
    bool useKarczewski = false;
    bool useExactAbsorption = false;  ///< Use DiffractInclineAbs (gradient absorption across aperture)
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

