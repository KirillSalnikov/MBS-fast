#pragma once

#include "Beam.h"
#include "Scattering.h"
#include "PhysMtr.hpp"
#include "MullerMatrix.h"
#include "Tracks.h"

#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

class ScatteringRange
{
public:
    ScatteringRange(double zenStart, double zenEnd, int nAz, int nZen)
        : zenithStart(zenStart), zenithEnd(zenEnd), nAzimuth(nAz), nZenith(nZen),
          isNonUniform(false)
    {
        azinuthStep = M_2PI/nAz;
        zenithStep = (zenEnd - zenStart)/nZen;
    }

    /// Load non-uniform theta grid from file (theta values in degrees, one per line).
    /// Keeps nAzimuth from the --grid option. Overrides nZenith, zenithStart, zenithEnd, zenithStep.
    bool LoadThetaGrid(const std::string &filename)
    {
        std::ifstream f(filename);
        if (!f.is_open()) return false;

        thetaValues.clear();
        double val;
        std::string line;
        while (std::getline(f, line))
        {
            // skip empty lines and comments
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            if (iss >> val)
                thetaValues.push_back(DegToRad(val));
        }
        f.close();

        if (thetaValues.size() < 2) return false;

        // Sort just in case
        std::sort(thetaValues.begin(), thetaValues.end());

        isNonUniform = true;
        nZenith = (int)thetaValues.size() - 1;
        zenithStart = thetaValues.front();
        zenithEnd = thetaValues.back();
        // zenithStep is meaningless for non-uniform, but set to average for safety
        zenithStep = (zenithEnd - zenithStart) / nZenith;

        return true;
    }

    /// Get zenith angle for bin j (works for both uniform and non-uniform grids)
    double GetZenith(int j) const
    {
        if (isNonUniform)
            return thetaValues[j];
        else
            return zenithStart + j * zenithStep;
    }

    /// Compute 2pi*dcos for bin j (works for both uniform and non-uniform grids)
    double Compute2PiDcos(int j) const
    {
        double radZen = GetZenith(j);

        if (!isNonUniform)
        {
            // Original uniform code
            double _2Pi_dcos = (radZen > M_PI-__FLT_EPSILON__ || radZen < __FLT_EPSILON__) ?
                        1.0-cos(0.5*zenithStep) :
                        cos((j-0.5)*zenithStep)-cos((j+0.5)*zenithStep);
            return _2Pi_dcos * M_2PI;
        }

        // Non-uniform grid: compute bin boundaries at midpoints between neighbors
        double theta_lo, theta_hi;

        if (j == 0)
        {
            theta_lo = thetaValues[0];
            theta_hi = 0.5 * (thetaValues[0] + thetaValues[1]);
        }
        else if (j == nZenith)
        {
            theta_lo = 0.5 * (thetaValues[nZenith - 1] + thetaValues[nZenith]);
            theta_hi = thetaValues[nZenith];
        }
        else
        {
            theta_lo = 0.5 * (thetaValues[j - 1] + thetaValues[j]);
            theta_hi = 0.5 * (thetaValues[j] + thetaValues[j + 1]);
        }

        double _2Pi_dcos = cos(theta_lo) - cos(theta_hi);
        return _2Pi_dcos * M_2PI;
    }

    void ComputeSphereDirections(const Light &incidentLight)
    {
        double sinAz;
        double cosAz;

        double sinZen;
        double cosZen;

        directions.clear();
        vf.clear();
        for (int i = 0; i <= nAzimuth; ++i)
        {
            double az = i * azinuthStep;
            sincos(az, &sinAz, &cosAz);

            directions.push_back(std::vector<Point3d>());
            vf.push_back(std::vector<Point3d>());

            for (int j = 0; j <= nZenith; ++j)
            {
                double zen = GetZenith(j);
                sincos(zen, &sinZen, &cosZen);
#ifdef _DEBUG // DEB
//                if (i == 1 && j == )
//                    int fff = 0;
#endif
                Point3d dir(sinZen*cosAz, sinZen*sinAz, -cosZen);
                directions[i].push_back(dir);

                if (dir.z >= 1-DBL_EPSILON)
                {
                    vf[i].push_back(-incidentLight.polarizationBasis);
                }
                else if (dir.z <= DBL_EPSILON-1)
                {
                    vf[i].push_back(incidentLight.polarizationBasis);
                }
                else
                {
                    vf[i].push_back(Point3d(-sinAz ,cosAz ,0));
                }
            }
        }
    }

public:
    double zenithStart;
    double zenithEnd;
    int nAzimuth;
    int nZenith;
    double azinuthStep;
    double zenithStep;

    bool isNonUniform;
    std::vector<double> thetaValues;  ///< Non-uniform theta grid (radians), sorted

    std::vector<std::vector<Point3d>> vf;
    std::vector<std::vector<Point3d>> directions;
};

class PointContribution
{
public:
    PointContribution(size_t nGroups, double normIndex)
        : m_nGroups(nGroups),
          m_normIndex(normIndex)
    {
        groupJones.resize(nGroups);
        groupMuellers.resize(nGroups);
        ResetJones();
    }

    void AddToMueller(const Matrix2x2c &jones)
    {
        MuellerMatrix m(jones);
        m *= m_normIndex;
//#ifdef _DEBUG // DEB
//        double ddd = m(0,0);
//        if (isnan(ddd))
//            int fff = 0;
//#endif
        rest += m;
    }

    void AddToGroup(const Matrix2x2c &jones, size_t groupId)
    {
        groupJones[groupId] += jones;
    }

    void SumGroupTotal()
    {
        for (size_t gr = 0; gr < m_nGroups; ++gr)
        {
            MuellerMatrix m(groupJones[gr]);
            m *= m_normIndex;
            groupMuellers[gr] += m;
            groupTotal += m;
        }

        ResetJones();
    }

    void SumTotal()
    {
        total += groupTotal;
        total += rest;
    }

    const MuellerMatrix &GetGroupTotal() const
    {
        return groupTotal;
    }

    const MuellerMatrix &GetTotal() const
    {
        return total;
    }

    const MuellerMatrix &GetRest() const
    {
        return rest;
    }

    const MuellerMatrix &GetGroupMueller(size_t groupID)
    {
        return groupMuellers.at(groupID);
    }

    void Reset()
    {
        groupTotal.Reset();
        rest.Reset();
        total.Reset();

        for (MuellerMatrix &m : groupMuellers)
        {
            m.Reset();
        }
    }

private:
    std::vector<Matrix2x2c> groupJones;
    std::vector<MuellerMatrix> groupMuellers;

    MuellerMatrix groupTotal;
    MuellerMatrix rest;
    MuellerMatrix total;

    size_t m_nGroups;
    double m_normIndex;

    void ResetJones()
    {
        for (Matrix2x2c &j : groupJones)
        {
            j.Fill(0.f);
        }
    }
};


class ContributionGO
{
public:
    ContributionGO()
        : muellers(0, 0, 0, 0),
          back(4, 4),
          forward(4, 4)
    {
        back.Fill(0);
        forward.Fill(0);
    }

    void SetStep(int nTh, double radius)
    {
        nTheta = nTh;
        thetaStep = RadToDeg(radius)/nTheta;
        muellers = Arr2D(1, nTheta + 1, 4, 4);
        muellers.ClearArr();
    }

    void AddMueller(float fAngle, const matrix &m)
    {
        double degA = RadToDeg(fAngle);

        int iCell = lround(degA/thetaStep);

        if (degA >= 180-FLT_EPSILON)
        {
            forward += m;
        }
        else if (degA <= FLT_EPSILON)
        {
            back += m;
        }
        else if (iCell < nTheta + 1)
        {
            muellers.insert(0, iCell, m);
#ifdef _DEBUG
            double ddd = muellers(0, iCell, 0, 1);
//            if (isnan(ddd))
           int fff = 0;
#endif
        }
    }

    int nTheta;
    float thetaStep;
    Arr2D muellers;		///< Scattering matrices
    matrix back;		///< Mueller matrix in backward direction
    matrix forward;		///< Mueller matrix in forward direction
};

struct BeamInfo
{
    bool order;
    bool isBad = false;
    double area;
    double projLenght;
    double opticalLengths[3];
    Point3f beamBasis;
    Point3d center;
    Point3d projectedCenter;
    Point3f normal;
    Point3d normald;
    Point3d horAxis;
    Point3d verAxis;
    Point3d lenIndices;
};

/// Precomputed polarization data for fast RotateJones (computed once per beam)
struct BeamPolData
{
    Point3d NTd;    // normal × beamBasis
    Point3d NPd;    // normal × polarizationBasis
    Point3d nxDT;   // normal × (beamDir × beamBasis)
    Point3d nxDP;   // normal × (beamDir × polarizationBasis)
};

/// Precomputed edge data for fast diffraction integral (computed once per beam)
struct BeamEdgeData
{
    static const int MAX_EDGES = 32;
    double x[MAX_EDGES];       ///< vertex x in aperture 2D
    double y[MAX_EDGES];       ///< vertex y in aperture 2D
    // Precomputed per-edge data (GOAD-style EdgeData)
    double dx[MAX_EDGES];      ///< x[next] - x[i]
    double dy[MAX_EDGES];      ///< y[next] - y[i]
    double slope_yx[MAX_EDGES]; ///< (y[i]-y[next])/(x[i]-x[next]) = -dy/dx
    double slope_xy[MAX_EDGES]; ///< (x[i]-x[next])/(y[i]-y[next]) = -dx/dy
    double intercept_y[MAX_EDGES]; ///< y[i] - slope_yx[i]*x[i]
    double intercept_x[MAX_EDGES]; ///< x[i] - slope_xy[i]*y[i]
    bool edge_valid_x[MAX_EDGES];  ///< |dx| > eps (usable for absB>absA branch)
    bool edge_valid_y[MAX_EDGES];  ///< |dy| > eps (usable for absA>=absB branch)
    int nVertices;
    bool valid;
};

class Handler
{
public:
    Handler(Particle *particle, Light *incidentLight, int nTheta,
            double wavelength);

    virtual double ComputeTotalScatteringEnergy() {}
    virtual void HandleBeams(std::vector<Beam> &beams, double sinZenith);
    virtual void SetTracks(Tracks *tracks);
    Tracks *GetTracks() const;
    void SetScattering(Scattering *scattering);
    virtual void WriteMatricesToFile(std::string &destName, double nrg);
    virtual void WriteTotalMatricesToFile(const std::string &destName);
    void SetAbsorptionAccounting(bool value);
    virtual void SetScatteringSphere(const ScatteringRange &grid);

    void SetNormIndex(double value);
    void SetSinZenith(double value);

    Light *m_incidentLight;

    double m_sinZenith;
    Tracks *m_tracks;

    int nTheta;
    int m_nBadBeams;
    bool m_isBadBeam;
    bool isCoh = true;
    ScatteringRange m_sphere;
    double normIndexGamma;
    double m_outputEnergy = 0;

    std::ofstream *betaFile;
    int m_fixedItr = -1;

protected:
    double BeamCrossSection(const Beam &beam) const;

    void ApplyAbsorption(Beam &beam);

    /**
     * @brief Calculate the diffraction of beam in given direction
     * @param beam beam of light
     * @param direction direction for calculating diffraction
     * @return Fresnel coefficient
     */
    complex DiffractIncline(const BeamInfo &info, const Beam &beam,
                            const Point3d &direction) const;

    void PrecomputeEdgeData(const BeamInfo &info, const Beam &beam,
                            BeamEdgeData &edgeData) const;
    complex DiffractInclineFast(const BeamInfo &info, const BeamEdgeData &edgeData,
                                const Point3d &beamDir, const Point3d &direction) const;

    complex DiffractInclineAbs(const BeamInfo &info, const Beam &beam,
                               const Point3d &direction) const;


    Point3d ChangeCoordinateSystem(const Point3d& hor, const Point3d& ver,
                                   const Point3d& normal,
                                   const Point3d& point) const;

    Point3d ChangeCoordinateSystem(const Point3d& normal,
                                   const Point3d &point) const;

    void ComputeCoordinateSystemAxes(const Point3d& normal,
                                     Point3d &hor, Point3d &ver) const;

    void ComputeLengthIndices(const Beam &beam, BeamInfo &info);

protected:
    Scattering *m_scattering;

    Particle *m_particle;
    double m_wavelength; // must be double type!!!
    bool m_hasAbsorption;
    double m_normIndex;
    std::ofstream m_logFile;
    double m_cAbs;

    complex m_ri;
    double m_riIm;

    double m_waveIndex;
    double m_wi2;

    complex m_complWave;
    complex m_invComplWave;
    double m_absMag;

    double m_eps1;
    double m_eps2;
    double m_eps3;

public:
    /// Beam importance cutoff: skip diffraction if |J|²×area < this threshold.
    /// Default 0 = auto (eps² × totalBeamEnergy per orientation).
    /// Set via --beam_cutoff CLI to override.
    double m_beamCutoff = 0;
    /// Target accuracy for auto beam cutoff (set from --auto eps)
    double m_targetEps = 0; // was 0.01, set to 0 to match MBS-raw (no beam cutoff)
    bool m_legacySign = false; ///< Use old (+invComplWave) sign for forward direction
    void SetBeamCutoff(double val) { m_beamCutoff = val; }
    /// Set cutoff relative to geometric cross-section: threshold = eps × C_geo
    void SetBeamCutoffRelative(double eps, double C_geo) { m_beamCutoff = eps * C_geo; }

private:
    void ExtropolateOpticalLenght(Beam &beam, const std::vector<int> &tr);
};
