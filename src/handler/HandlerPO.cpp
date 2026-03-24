#include "HandlerPO.h"
#include "HandlerPO_fast.h"
#include "HandlerPO_avx.h"

#include "Mueller.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <chrono>

// ---- Timing accumulators (static, printed at program exit) ----
static double g_time_rotateJones_us    = 0;
static double g_time_diffractFast_us   = 0;
static double g_time_computeFnJones_us = 0;
static double g_time_matmul_us         = 0;
static double g_time_mueller_us        = 0;
static double g_time_addToMueller_us   = 0;
static double g_time_beamLoop_us       = 0;
static double g_time_handleBeams_us    = 0;
static long long g_count_directions    = 0;
static long long g_count_beams         = 0;
static int g_handleBeams_calls         = 0;

static struct TimingPrinter {
    ~TimingPrinter() {
        if (g_handleBeams_calls == 0) return;
        std::cerr << "\n===== HandlerPO TIMING BREAKDOWN =====\n";
        std::cerr << "HandleBeams calls: " << g_handleBeams_calls << "\n";
        std::cerr << "Total beams processed: " << g_count_beams << "\n";
        std::cerr << "Total direction evaluations: " << g_count_directions << "\n";
        std::cerr << std::fixed << std::setprecision(3);
        std::cerr << "HandleBeams total:    " << g_time_handleBeams_us/1e6 << " s\n";
        std::cerr << "  Beam loop total:    " << g_time_beamLoop_us/1e6 << " s\n";
        std::cerr << "  RotateJones:        " << g_time_rotateJones_us/1e6 << " s\n";
        std::cerr << "  DiffractInclineFast:" << g_time_diffractFast_us/1e6 << " s\n";
        std::cerr << "  ComputeFnJones:     " << g_time_computeFnJones_us/1e6 << " s\n";
        std::cerr << "  MatMul (fres*rot*fn):" << g_time_matmul_us/1e6 << " s\n";
        std::cerr << "  Mueller conversion: " << g_time_mueller_us/1e6 << " s\n";
        std::cerr << "  AddToMueller:       " << g_time_addToMueller_us/1e6 << " s\n";
        double accounted = g_time_rotateJones_us + g_time_diffractFast_us
                         + g_time_computeFnJones_us + g_time_matmul_us
                         + g_time_mueller_us + g_time_addToMueller_us;
        std::cerr << "  Accounted:          " << accounted/1e6 << " s\n";
        std::cerr << "  Unaccounted:        " << (g_time_beamLoop_us - accounted)/1e6 << " s\n";
        if (g_count_directions > 0) {
            std::cerr << "  Per direction (ns): " << (g_time_beamLoop_us*1000.0)/g_count_directions << "\n";
        }
        std::cerr << "======================================\n";
    }
} g_timingPrinter;

HandlerPO::HandlerPO(Particle *particle, Light *incidentLight, int nTheta,
                     double wavelength)
    : Handler(particle, incidentLight, nTheta, wavelength)
{
    m_Lp = new matrix(4, 4);

    (*m_Lp)[0][0] = 1;
    (*m_Lp)[0][1] = 0;
    (*m_Lp)[0][2] = 0;
    (*m_Lp)[0][3] = 0;

    (*m_Lp)[1][0] = 0;
    (*m_Lp)[1][3] = 0;

    (*m_Lp)[2][0] = 0;
    (*m_Lp)[2][3] = 0;

    (*m_Lp)[3][0] = 0;
    (*m_Lp)[3][1] = 0;
    (*m_Lp)[3][2] = 0;
    (*m_Lp)[3][3] = 1;

    m_Ln = new matrix(4, 4);

    (*m_Ln) = (*m_Lp);
}

void HandlerPO::CleanJ()
{
    m_diffractedMatrices.clear();
    Arr2DC tmp(m_sphere.nAzimuth + 1, m_sphere.nZenith + 1, 2, 2);
    tmp.ClearArr();

    for (unsigned q = 0; q < m_tracks->size() + 1; q++)
    {
        m_diffractedMatrices.push_back(tmp);
    }
}

void HandlerPO::WriteMatricesToFile(std::string &destName, double nrg)
{
//    if (!m_tracks->shouldComputeTracksOnly)
    {
        std::ofstream outFile(destName, std::ios::app);

        outFile /*<< std::to_string(m_sphere.radius) << ' '*/
                << std::to_string(m_sphere.nZenith) << ' '
                << std::to_string(m_sphere.nAzimuth+1);

        for (int t = 0; t <= m_sphere.nZenith; ++t)
        {
            double tt = RadToDeg(m_sphere.GetZenith(t));

            for (int p = 0; p <= m_sphere.nAzimuth; ++p)
            {
                double fi = -((double)p)*m_sphere.azinuthStep;
                double degPhi = RadToDeg(-fi);
                outFile << std::endl << tt << " " << degPhi << " " << nrg << " ";

                matrix m = M(p, t);
                outFile << m;
            }
        }
    }
}

void HandlerPO::WriteGroupMatrices(Arr2D &matrices, const std::string &name)
{
    auto &Lp = *m_Lp;
    auto &Ln = *m_Ln;

    std::ofstream outFile(name, std::ios::out);

    outFile /*<< std::to_string(m_sphere.radius) << ' '*/
            << std::to_string(m_sphere.nZenith) << ' '
            << std::to_string(m_sphere.nAzimuth+1);

    matrix sum(4, 4);

    for (int t = m_sphere.nZenith; t >= 0; --t)
    {
        sum.Fill(0.0);
        double tt = RadToDeg(m_sphere.zenithEnd - m_sphere.GetZenith(t));

        for (int p = 0; p <= m_sphere.nAzimuth; ++p)
        {
            double fi = -((double)p)*m_sphere.azinuthStep;
            matrix m = matrices(p, t);

            Lp[1][1] = cos(2*fi);
            Lp[1][2] = sin(2*fi);
            Lp[2][1] = -Lp[1][2];
            Lp[2][2] = Lp[1][1];

            Ln[1][2] = -Lp[1][2];
            Ln[2][1] = -Lp[2][1];

            if (t == 0)
            {
                sum += Lp*m*Lp;
            }
            else if (t == m_sphere.nZenith-1)
            {
                sum += Ln*m*Lp; // OPT: вынести Ln в отдельный случай
            }
            else
            {
                sum += m*Lp;
            }
        }

        outFile << std::endl << tt << " ";
        outFile << sum/m_sphere.nAzimuth;
    }
}

void HandlerPO::WriteJonesToFile(const std::string &destName)
{
    std::string jonesFile = destName + "_jones.dat";
    std::ofstream outFile(jonesFile, std::ios::out);
    outFile << std::scientific << std::setprecision(15);

    outFile << "# theta(deg) phi(deg)"
            << " J00_re J00_im J01_re J01_im"
            << " J10_re J10_im J11_re J11_im" << std::endl;

    for (int t = 0; t <= m_sphere.nZenith; ++t)
    {
        double theta_deg = RadToDeg(m_sphere.GetZenith(t));

        for (int p = 0; p <= m_sphere.nAzimuth; ++p)
        {
            double phi_deg = RadToDeg(p * m_sphere.azinuthStep);

            // Sum Jones matrices from all groups (coherent sum)
            complex j00(0, 0), j01(0, 0), j10(0, 0), j11(0, 0);

            for (size_t q = 0; q < m_diffractedMatrices.size(); ++q)
            {
                matrixC jq = m_diffractedMatrices[q](p, t);
                j00 += jq[0][0];
                j01 += jq[0][1];
                j10 += jq[1][0];
                j11 += jq[1][1];
            }

            outFile << theta_deg << " " << phi_deg
                    << " " << real(j00) << " " << imag(j00)
                    << " " << real(j01) << " " << imag(j01)
                    << " " << real(j10) << " " << imag(j10)
                    << " " << real(j11) << " " << imag(j11)
                    << std::endl;
        }
    }
}

void HandlerPO::WriteTotalMatricesToFile(const std::string &destName)
{
    if (m_tracks->shouldOutputGroups)
    {
        for (size_t i = 0; i < m_diffractedMatrices.size(); ++i)
        {
            if ((*m_tracks)[i].size != 0)
            {
                const std::string subname = (*m_tracks)[i].CreateGroupName();
                const std::string &filename = destName + '_' +  subname;
                WriteGroupMatrices(m_groupMatrices[i], filename);
            }
        }

        WriteGroupMatrices(M, destName + "_total.dat");
    }
}

// double HandlerPO::ComputeTotalScatteringEnergy()
// {
//     double D_tot = 0;

//     for (int t = 0; t <= m_sphere.nZenith; ++t)
//     {
//         for (int p = 0; p <= m_sphere.nAzimuth; ++p)
//         {
//             matrix m = M(p, t);
//             D_tot += m[0][0];
//         }
//     }

//     return D_tot;
// }

void HandlerPO::RotateJones(const Beam &beam, const BeamInfo &info,
                            const Vector3d &vf, const Vector3d &direction,
                            matrixC &matrix) const
{
    auto &dir = direction;

    Vector3d vt = CrossProductD(vf, dir);
    vt = vt/LengthD(vt);

    Vector3f NT = CrossProduct(info.normal, info.beamBasis);
    Vector3f NE = CrossProduct(info.normal, beam.polarizationBasis);

    Vector3d NTd = Vector3d(NT.cx, NT.cy, NT.cz);
    Vector3d NPd = Vector3d(NE.cx, NE.cy, NE.cz);

    Point3f DT = CrossProduct(beam.direction, info.beamBasis);
    Point3f DP = CrossProduct(beam.direction, beam.polarizationBasis);

    Point3d DTd = Point3d(DT.cx, DT.cy, DT.cz);
    Point3d DPd = Point3d(DP.cx, DP.cy, DP.cz);

    const Point3d nd = info.normal;

    // BAC-CAB: dir × (dir × X) = dir(dir·X) - X  (dir is unit vector)
    // So: dir×NT - dir×(dir×(n×DT)) = dir×NT + (n×DT) - dir*(dir·(n×DT))
    Point3d nxDT = CrossProductD(nd, DTd);
    Point3d cpT = CrossProductD(dir, NTd) + nxDT - dir * DotProductD(dir, nxDT);

    Point3d nxDP = CrossProductD(nd, DPd);
    Point3d cpP = CrossProductD(dir, NPd) + nxDP - dir * DotProductD(dir, nxDP);

    matrix[0][0] = DotProductD(cpT, vt)/2.0;
    matrix[0][1] = DotProductD(cpP, vt)/2.0;
    matrix[1][0] = DotProductD(cpT, vf)/2.0;
    matrix[1][1] = DotProductD(cpP, vf)/2.0;
}

void HandlerPO::PrecomputePolData(const Beam &beam, const BeamInfo &info,
                                   BeamPolData &pd)
{
    Vector3f NT = CrossProduct(info.normal, info.beamBasis);
    Vector3f NE = CrossProduct(info.normal, beam.polarizationBasis);
    pd.NTd = Vector3d(NT.cx, NT.cy, NT.cz);
    pd.NPd = Vector3d(NE.cx, NE.cy, NE.cz);

    Point3f DT = CrossProduct(beam.direction, info.beamBasis);
    Point3f DP = CrossProduct(beam.direction, beam.polarizationBasis);
    Point3d DTd(DT.cx, DT.cy, DT.cz);
    Point3d DPd(DP.cx, DP.cy, DP.cz);

    Point3d nd(info.normal.cx, info.normal.cy, info.normal.cz);
    pd.nxDT = CrossProductD(nd, DTd);
    pd.nxDP = CrossProductD(nd, DPd);
}

void HandlerPO::RotateJonesFast(const BeamPolData &pd,
                                 const Vector3d &vf, const Vector3d &dir,
                                 matrixC &matrix)
{
    Vector3d vt = CrossProductD(vf, dir);
    double vtLen = LengthD(vt);
    if (vtLen > 1e-15) vt = vt / vtLen;

    // cpT = dir×NTd + nxDT - dir*(dir·nxDT)
    Point3d cpT = CrossProductD(dir, pd.NTd) + pd.nxDT - dir * DotProductD(dir, pd.nxDT);
    Point3d cpP = CrossProductD(dir, pd.NPd) + pd.nxDP - dir * DotProductD(dir, pd.nxDP);

    matrix[0][0] = complex(DotProductD(cpT, vt) * 0.5, 0);
    matrix[0][1] = complex(DotProductD(cpP, vt) * 0.5, 0);
    matrix[1][0] = complex(DotProductD(cpT, vf) * 0.5, 0);
    matrix[1][1] = complex(DotProductD(cpP, vf) * 0.5, 0);
}

void HandlerPO::KarczewskiJones(const Beam &beam, const BeamInfo &info,
                                 const Vector3d &vf_out, const Vector3d &direction,
                                 matrixC &matrix) const
{
    // =========================================================================
    // GOAD-compatible Karczewski pipeline (replicates diff2.rs)
    // =========================================================================
    //
    // IMPORTANT: This implementation computes the Karczewski matrix in aperture
    // frame but does NOT fully replicate GOAD's 5-step coordinate pipeline.
    // As a result:
    //   - M11 is IDENTICAL to RotateJones (proven: ||R|| = ||K||, same Frobenius
    //     norm, so M11 = (1/2)||A||^2_F is the same regardless of formalism).
    //   - M33, M34, M44 DIFFER from GOAD at phi != 0 because the matrix STRUCTURE
    //     differs (off-diagonal elements reach ±0.87 for RotateJones vs ±0.001
    //     for Karczewski).
    //   - Fixing M33/M34/M44 to match GOAD would require replicating the full
    //     coordinate transformation pipeline, which is extremely complex.
    //
    // STATUS: Experimental. Use for M33/M34/M44 research only.
    // For M11: results are identical to default RotateJones.
    // =========================================================================
    //
    // GOAD pipeline:
    //   1. Build rot3: lab -> aperture frame (face in xy, beam e_perp along +y)
    //   2. Transform prop and k_obs into aperture frame
    //   3. Compute Karczewski matrix + e_perp (m vector) in aperture frame
    //   4. Compute rot4: rotation from Karczewski basis to scattering plane basis
    //   5. Compute prerotation: rotation from incident beam frame to scattering plane
    //   6. result = rot4 × Karczewski × prerotation  (ampl is handled separately)
    //
    // In MBS, beam.J (the amplitude matrix) is multiplied in ApplyDiffraction,
    // so here we only return: rot4 × Karczewski × prerotation
    // =========================================================================

    // --- Step 1: Build aperture coordinate system ---
    // z_ap = face normal (pointing in beam propagation direction)
    Point3d z_ap = info.normald;

    // Project beam's polarizationBasis (e_perp) onto the face plane -> y_ap
    Point3f bp = beam.polarizationBasis;
    Point3d e_perp_beam(bp.cx, bp.cy, bp.cz);
    double dot_en = DotProductD(e_perp_beam, z_ap);
    Point3d y_ap(e_perp_beam.x - dot_en*z_ap.x,
                 e_perp_beam.y - dot_en*z_ap.y,
                 e_perp_beam.z - dot_en*z_ap.z);
    double y_len = LengthD(y_ap);
    if (y_len < 1e-12) {
        // Fallback: beam e_perp is parallel to normal, use MBS default axes
        y_ap = info.verAxis;
        y_len = LengthD(y_ap);
    }
    y_ap = y_ap / y_len;

    // x_ap completes right-handed system: x = y × z
    Point3d x_ap = CrossProductD(y_ap, z_ap);
    double x_len = LengthD(x_ap);
    if (x_len > 1e-12) x_ap = x_ap / x_len;

    // rot3 matrix: transforms lab vector v to aperture frame
    // rot3 * v = (v·x_ap, v·y_ap, v·z_ap)
    // Row 0 = x_ap, Row 1 = y_ap, Row 2 = z_ap

    // --- Step 2: Transform directions into aperture frame ---
    Point3f bd = beam.direction;
    Point3d prop(bd.cx*x_ap.x + bd.cy*x_ap.y + bd.cz*x_ap.z,
                 bd.cx*y_ap.x + bd.cy*y_ap.y + bd.cz*y_ap.z,
                 bd.cx*z_ap.x + bd.cy*z_ap.y + bd.cz*z_ap.z);

    Point3d k(direction.x*x_ap.x + direction.y*x_ap.y + direction.z*x_ap.z,
              direction.x*y_ap.x + direction.y*y_ap.y + direction.z*y_ap.z,
              direction.x*z_ap.x + direction.y*z_ap.y + direction.z*z_ap.z);

    // --- Step 3: Karczewski matrix in aperture frame ---
    double Kx = prop.x, Ky = prop.y, Kz = prop.z;
    double kx = k.x, ky = k.y, kz = k.z;

    double one_minus_ky2 = 1.0 - ky*ky;
    if (one_minus_ky2 < 1e-12) one_minus_ky2 = 1e-12;
    double r_k = sqrt(one_minus_ky2);

    double one_minus_Ky2 = 1.0 - Ky*Ky;
    if (one_minus_Ky2 < 1e-12) one_minus_Ky2 = 1e-12;

    double frac = sqrt(one_minus_ky2 / one_minus_Ky2);
    if (fabs(frac) < 1e-12) frac = 1e-12;

    double a1m = -Kz * frac;
    double b2m = -kz / frac;
    double a1e = b2m;
    double b2e = a1m;
    double b1m = -kx*ky/frac + Kx*Ky*frac;
    double a2e = -b1m;

    double D00 = 0.5*(a1m + a1e);
    double D01 = 0.5*b1m;
    double D10 = 0.5*a2e;
    double D11 = 0.5*(b2m + b2e);

    // Karczewski e_perp (m vector) in aperture frame
    Point3d m_vec(-kx*ky/r_k, r_k, -ky*kz/r_k);

    // --- Step 4: rot4 in aperture frame ---
    // GOAD: hc = rot3 * scattering_e_perp (vf_out in lab)
    // Transform vf_out (scattering e_perp, lab frame) into aperture frame
    Point3d hc(vf_out.x*x_ap.x + vf_out.y*x_ap.y + vf_out.z*x_ap.z,
               vf_out.x*y_ap.x + vf_out.y*y_ap.y + vf_out.z*y_ap.z,
               vf_out.x*z_ap.x + vf_out.y*z_ap.y + vf_out.z*z_ap.z);

    // GOAD rotation_matrix(e_perp_in=m_vec, e_perp_out=hc, prop=k):
    //   dot1 = hc · m_vec
    //   evo2 = k × m_vec, normalized
    //   dot2 = hc · evo2
    //   result = [[dot1, dot2], [-dot2, dot1]] / sqrt(|det|)
    double dot1 = DotProductD(hc, m_vec);
    Point3d evo2 = CrossProductD(k, m_vec);
    double evo2_len = LengthD(evo2);
    if (evo2_len > 1e-12) evo2 = evo2 / evo2_len;
    double dot2 = DotProductD(hc, evo2);

    double det4 = dot1*dot1 + dot2*dot2;
    double inv_sqrt_det4 = (det4 > 1e-12) ? 1.0/sqrt(det4) : 1.0;

    double r4_00 = dot1*inv_sqrt_det4, r4_01 = dot2*inv_sqrt_det4;
    double r4_10 = -dot2*inv_sqrt_det4, r4_11 = dot1*inv_sqrt_det4;

    // --- Step 5: prerotation ---
    // GOAD: prerotation = rotation_matrix(incident.e_perp, scattering_e_perp, incident.prop)^T
    // incident.e_perp = m_incidentLight->polarizationBasis (lab frame)
    // scattering_e_perp = vf_out (lab frame)
    // incident.prop = -m_incidentLight->direction (GOAD uses prop pointing along beam travel)
    //
    // In MBS, beam.polarizationBasis has already been rotated from incident
    // to the beam's own frame. The Jones matrix beam.J encodes this.
    // But KarczewskiJones is a replacement for RotateJones which does the
    // coordinate transformation from beam frame to scattering frame.
    // The prerotation in GOAD handles the same thing: rotating from the
    // incident beam's polarization frame into the scattering plane frame.
    //
    // MBS incident beam: polarizationBasis and direction are in lab frame
    Point3d inc_e_perp(m_incidentLight->polarizationBasis.cx,
                       m_incidentLight->polarizationBasis.cy,
                       m_incidentLight->polarizationBasis.cz);
    Point3d inc_prop(-m_incidentLight->direction.cx,
                     -m_incidentLight->direction.cy,
                     -m_incidentLight->direction.cz);

    // rotation_matrix(e_perp_in=inc_e_perp, e_perp_out=vf_out, prop=inc_prop)
    double pre_dot1 = DotProductD(vf_out, inc_e_perp);
    Point3d pre_evo2 = CrossProductD(inc_prop, inc_e_perp);
    double pre_evo2_len = LengthD(pre_evo2);
    if (pre_evo2_len > 1e-12) pre_evo2 = pre_evo2 / pre_evo2_len;
    double pre_dot2 = DotProductD(vf_out, pre_evo2);

    double pre_det = pre_dot1*pre_dot1 + pre_dot2*pre_dot2;
    double pre_norm = (pre_det > 1e-12) ? 1.0/sqrt(pre_det) : 1.0;

    // rotation_matrix result (before transpose)
    double pr_00 = pre_dot1*pre_norm, pr_01 = pre_dot2*pre_norm;
    double pr_10 = -pre_dot2*pre_norm, pr_11 = pre_dot1*pre_norm;

    // Transpose for prerotation
    double pre_00 = pr_00, pre_01 = pr_10;
    double pre_10 = pr_01, pre_11 = pr_11;

    // --- Step 6: Combine: rot4 × Karczewski × prerotation ---
    // First: temp = rot4 × Karczewski
    double t00 = r4_00*D00 + r4_01*D10;
    double t01 = r4_00*D01 + r4_01*D11;
    double t10 = r4_10*D00 + r4_11*D10;
    double t11 = r4_10*D01 + r4_11*D11;

    // Then: result = temp × prerotation
    matrix[0][0] = complex(t00*pre_00 + t01*pre_10, 0);
    matrix[0][1] = complex(t00*pre_01 + t01*pre_11, 0);
    matrix[1][0] = complex(t10*pre_00 + t11*pre_10, 0);
    matrix[1][1] = complex(t10*pre_01 + t11*pre_11, 0);
}

matrixC HandlerPO::ApplyDiffractionFast(const Beam &beam, const BeamInfo &info,
                                        const BeamEdgeData &edgeData,
                                        const Point3d &beamDirD,
                                        const Vector3d &direction,
                                        const Vector3d &vf)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    matrixC fnJones = (beam.lastFacetId != __INT_MAX__) ?
                ComputeFnJones(beam.J, info, direction) :
                beam.J * fast_exp_im(m_waveIndex*beam.opticalPath);
    auto t1 = std::chrono::high_resolution_clock::now();

    matrixC jones_rot(2, 2);
    RotateJones(beam, info, vf, direction, jones_rot);
    auto t2 = std::chrono::high_resolution_clock::now();

    complex fresnel = DiffractInclineFast(info, edgeData, beamDirD, direction);
    auto t3 = std::chrono::high_resolution_clock::now();

    if (isnan(real(fresnel)))
        return matrixC(2, 2);

    matrixC result = fresnel*jones_rot*fnJones;
    auto t4 = std::chrono::high_resolution_clock::now();

    g_time_computeFnJones_us += std::chrono::duration<double, std::micro>(t1 - t0).count();
    g_time_rotateJones_us    += std::chrono::duration<double, std::micro>(t2 - t1).count();
    g_time_diffractFast_us   += std::chrono::duration<double, std::micro>(t3 - t2).count();
    g_time_matmul_us         += std::chrono::duration<double, std::micro>(t4 - t3).count();

    return result;
}

matrixC HandlerPO::ApplyDiffractionFast2(const Beam &beam, const BeamInfo &info,
                                         const BeamEdgeData &edgeData,
                                         const Point3d &beamDirD,
                                         const matrixC &J_phased,
                                         bool isExternal,
                                         const Vector3d &direction,
                                         const Vector3d &vf)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    // J_phased already contains beam.J * fast_exp_im(k * projLength) [internal]
    // or beam.J * fast_exp_im(k * opticalPath) [external].
    // For internal beams, we only need the direction-dependent phase:
    //   exp_im(-k * dot(direction, center))
    // For external beams, J_phased is already complete (no direction dependence).
    complex dirPhase;
    if (!isExternal) {
        double dp = DotProductD(direction, info.center);
        dirPhase = fast_exp_im(-m_waveIndex * dp);
    } else {
        dirPhase = complex(1.0, 0.0);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    matrixC jones_rot(2, 2);
    RotateJones(beam, info, vf, direction, jones_rot);

    complex fresnel = DiffractInclineFast(info, edgeData, beamDirD, direction);

    if (isnan(real(fresnel)))
        return matrixC(2, 2);

    complex combined_phase = fresnel * dirPhase;
    matrixC scaled_rot = jones_rot * combined_phase;
    matrixC result = scaled_rot * J_phased;

    return result;
}

// Fastest version: precomputed polData + edgeData + J_phased
matrixC HandlerPO::ApplyDiffractionFast3(const BeamPolData &polData,
                                          const BeamInfo &info,
                                          const BeamEdgeData &edgeData,
                                          const Point3d &beamDirD,
                                          const matrixC &J_phased,
                                          bool isExternal,
                                          const Vector3d &direction,
                                          const Vector3d &vf)
{
    // Direction-dependent phase
    complex dirPhase = isExternal ? complex(1.0, 0.0)
                                  : fast_exp_im(-m_waveIndex * DotProductD(direction, info.center));

    // Polarization rotation (uses precomputed beam-independent data)
    matrixC jones_rot(2, 2);
    RotateJonesFast(polData, vf, direction, jones_rot);

    // Diffraction integral (uses precomputed edge vertices)
    complex fresnel = DiffractInclineFast(info, edgeData, beamDirD, direction);

    if (isnan(real(fresnel)))
        return matrixC(2, 2);

    complex combined_phase = fresnel * dirPhase;
    matrixC scaled_rot = jones_rot * combined_phase;
    return scaled_rot * J_phased;
}

matrixC HandlerPO::ApplyDiffraction(const Beam &beam, const BeamInfo &info,
                                    const Vector3d &direction,
                                    const Vector3d &vf)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    matrixC fnJones = (beam.lastFacetId != __INT_MAX__) ?
                ComputeFnJones(beam.J, info, direction) :
                beam.J * fast_exp_im(m_waveIndex*beam.opticalPath);
    auto t1 = std::chrono::high_resolution_clock::now();

    matrixC jones_rot(2, 2);
    if (useKarczewski)
        KarczewskiJones(beam, info, vf, direction, jones_rot);
    else
        RotateJones(beam, info, vf, direction, jones_rot);
    auto t2 = std::chrono::high_resolution_clock::now();

    complex fresnel = DiffractIncline(info, beam, direction);
    auto t3 = std::chrono::high_resolution_clock::now();
//	complex fresnel = (m_hasAbsorption && beam.lastFacetId != INT_MAX && beam.nActs > 0)
//			? DiffractInclineAbs(info, beam, direction)
//			: DiffractIncline(info, beam, direction);

    if (isnan(real(fresnel)))
    {
        return matrixC(2, 2);
    }

    matrixC difracted = fresnel*jones_rot*fnJones;
    auto t4 = std::chrono::high_resolution_clock::now();

    g_time_computeFnJones_us += std::chrono::duration<double, std::micro>(t1 - t0).count();
    g_time_rotateJones_us    += std::chrono::duration<double, std::micro>(t2 - t1).count();
    g_time_diffractFast_us   += std::chrono::duration<double, std::micro>(t3 - t2).count();
    g_time_matmul_us         += std::chrono::duration<double, std::micro>(t4 - t3).count();
#ifdef _DEBUG // DEB
    complex ddd[4];
    ddd[0] = jones_rot[0][0];
    ddd[1] = jones_rot[0][1];
    ddd[2] = jones_rot[1][0];
    ddd[3] = jones_rot[1][1];

    complex qqq[4];
    qqq[0] = fnJones[0][0];
    qqq[1] = fnJones[0][1];
    qqq[2] = fnJones[1][0];
    qqq[3] = fnJones[1][1];

    complex bbb[4];
    bbb[0] = difracted[0][0];
    bbb[1] = difracted[0][1];
    bbb[2] = difracted[1][0];
    bbb[3] = difracted[1][1];
#endif
    return difracted;
}

matrixC HandlerPO::ComputeFnJones(const Matrix2x2c &matrix, const BeamInfo &info,
                                  const Vector3d &direction)
{
    double dp = DotProductD(direction, info.center);
    double arg = m_waveIndex*(info.projLenght - dp);
    return matrix*fast_exp_im(arg);
}

void HandlerPO::AddToMueller()
{
    for (size_t q = 0; q < m_diffractedMatrices.size(); ++q)
    {
        for (int t = 0; t <= m_sphere.nZenith; ++t)
        {
            for (int p = 0; p <= m_sphere.nAzimuth; ++p)
            {
                matrix m = Mueller(m_diffractedMatrices[q](p, t));
                m *= m_normIndex;
                M.insert(p, t, m);
                m_groupMatrices[q].insert(p, t, m);
            }
        }
    }
}

BeamInfo HandlerPO::ComputeBeamInfo(Beam &beam)
{
    BeamInfo info;
    info.normal = beam.Normal();
    info.normald = Point3d(info.normal.cx, info.normal.cy, info.normal.cz);

    info.order = DotProduct(info.normal, beam.direction) > 0;

    if (!info.order)
    {
        //info.normal = -info.normal;
//        info.normald = -info.normald;
    }

    m_isBadBeam = false;

    ComputeCoordinateSystemAxes(info.normald, info.horAxis, info.verAxis);

    info.center = beam.Center();
    info.projectedCenter = ChangeCoordinateSystem(info.horAxis, info.verAxis,
                                                  info.normald, info.center);

    if (m_hasAbsorption && beam.lastFacetId != __INT_MAX__ && beam.nActs > 1)
    {
        ComputeOpticalLengths(beam, info);
        ComputeLengthIndices(beam, info);
    }

    info.area = beam.Area();

    info.projLenght = beam.opticalPath + DotProductD(info.center, beam.direction);
//#ifdef _DEBUG // DEB
//	std::vector<int> tr;
//	Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

//	double sss = m_scattering->ComputeInternalOpticalPath(
//				beam, beam.Center(), tr);
//#endif
    info.beamBasis = CrossProduct(beam.polarizationBasis, beam.direction);
    info.beamBasis = info.beamBasis/Length(info.beamBasis); // basis of beam

    return info;
}

void HandlerPO::SetBackScatteringConus(double radAngle)
{
    isBackScatteringConusEnabled = true;
    backScatteringConus = cos(radAngle);
}

void HandlerPO::SetScatteringSphere(const ScatteringRange &grid)
{
    m_sphere = grid;
    M = Arr2D(m_sphere.nAzimuth+1, m_sphere.nZenith+1, 4, 4);
    M_noshadow = Arr2D(m_sphere.nAzimuth+1, m_sphere.nZenith+1, 4, 4);

    m_sphere.ComputeSphereDirections(*m_incidentLight);
}

void HandlerPO::ComputeOpticalLengths(const Beam &beam, BeamInfo &info)
{
    std::vector<int> tr;
    Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

    for (int i = 0; i < 3; ++i)
    {
        info.opticalLengths[i] = m_scattering->ComputeInternalOpticalPath(
                    beam, beam.arr[i], tr);
    }

//#ifdef _DEBUG // DEB
//	if (info.opticalLengths[0] > 200)
//		int fff = 0;
//#endif

//	double path = m_scattering->ComputeInternalOpticalPath(beam, beam.Center(), tr);

//	if (path > DBL_EPSILON)
//	{
//		double abs = exp(m_cAbs*path);
//		beam.J *= abs;
//	}

    //	ExtropolateOpticalLenght(beam, tr);
}

void HandlerPO::HandleBeams(std::vector<Beam> &beams, double sinZenith)
{
    auto hb_start = std::chrono::high_resolution_clock::now();
    ++g_handleBeams_calls;

    m_sinZenith = sinZenith;
    double sum = 0;
    int groupId;
    CleanJ();

    // --- BEAM LOGGING (first call only) ---
    static int handleBeamsCallCount = 0;
    ++handleBeamsCallCount;
    if (handleBeamsCallCount == 1)
    {
        std::ofstream logFile("beam_log.dat");
        logFile << "# beam_num  nActs  area  opticalPath  |Jones|  cx  cy  cz  nVertices\n";
        logFile << std::scientific << std::setprecision(8);
        int beamCounter = 0;
        for (const Beam &b : beams)
        {
            double jonesNorm = std::sqrt(b.J.Norm());
            logFile << beamCounter
                    << "  " << b.nActs
                    << "  " << b.Area()
                    << "  " << b.opticalPath
                    << "  " << jonesNorm
                    << "  " << b.direction.cx
                    << "  " << b.direction.cy
                    << "  " << b.direction.cz
                    << "  " << b.nVertices
                    << "\n";
            ++beamCounter;
        }
        logFile.close();
        std::cerr << "[beam_log] Logged " << beamCounter << " beams to beam_log.dat\n";
    }
    // --- END BEAM LOGGING ---

    auto beam_loop_start = std::chrono::high_resolution_clock::now();

    // Precompute sin(theta), cos(theta) for ALL theta bins (once for all beams)
    int nZen_global = m_sphere.nZenith;
    std::vector<double> sin_theta_arr(nZen_global+1), cos_theta_arr(nZen_global+1);
    for (int j = 0; j <= nZen_global; ++j) {
        double theta_rad = m_sphere.GetZenith(j);
        fast_sincos(theta_rad, sin_theta_arr[j], cos_theta_arr[j]);
    }
    // Precompute sin(phi), cos(phi) for all phi bins
    int nAz_global = m_sphere.nAzimuth;
    std::vector<double> sin_phi_arr(nAz_global+1), cos_phi_arr(nAz_global+1);
    for (int i = 0; i <= nAz_global; ++i) {
        double phi_rad = i * m_sphere.azinuthStep;
        fast_sincos(phi_rad, sin_phi_arr[i], cos_phi_arr[i]);
    }

    for (Beam &beam : beams)
    {
        //++g_count_beams;
        if (isBackScatteringConusEnabled && beam.direction.cz < backScatteringConus)
        {
            continue;
        }

        groupId = 0;
#ifdef _DEBUG // DEB
        // if (beam.id == 423)
        //     int dsdsdsd = 0;
//		std::vector<int> tr;
//		Tracks::RecoverTrack(beam, m_particle->nFacets, tr);
//		if (tr.size() == 2 && tr[0] == 2 && tr[1] == 4)
//			int fff = 0;
#endif
//		if (m_tracks->shouldComputeTracksOnly)
//		{
//			groupId = m_tracks->FindGroupByTrackId(beam.id);

            if (groupId < 0)
            {
//				groupId = 0;
//				continue;
            }
//		}

        beam.polarizationBasis = beam.RotateSpherical(
                    -m_incidentLight->direction,
                    m_incidentLight->polarizationBasis);

        BeamInfo info = ComputeBeamInfo(beam);

        if (m_isBadBeam)
        {
            continue;
        }

        if (m_hasAbsorption && beam.lastFacetId != __INT_MAX__ && beam.lastFacetId != -1)
        {
            ApplyAbsorption(beam);
        }

        if (beam.lastFacetId != __INT_MAX__)
        {
            matrix m_ = Mueller(beam.J);
            m_outputEnergy += BeamCrossSection(beam)*m_[0][0]*m_sinZenith;
        }

        // Skip diffraction for negligible beams (importance cutoff)
        // Physical criterion: beam contributes < ε of geometric cross-section
        // contribution = |J|² × area (proportional to beam's scattering power)
        // threshold = C_geo × ε_rel  (where C_geo = incoming energy ≈ projected area)
        if (beam.lastFacetId != __INT_MAX__)
        {
            double contribution = beam.J.Norm() * info.area;
            double threshold = m_beamCutoff;  // set via SetBeamCutoff()
            if (contribution < threshold)
                continue;
        }

        // Precompute ONCE per beam
        BeamEdgeData edgeData;
        PrecomputeEdgeData(info, beam, edgeData);
        Point3d beamDirD(beam.direction.cx, beam.direction.cy, beam.direction.cz);

        BeamPolData polData;
        PrecomputePolData(beam, info, polData);

        bool isExternal = (beam.lastFacetId == __INT_MAX__);
        matrixC J_phased = isExternal
            ? beam.J * fast_exp_im(m_waveIndex * beam.opticalPath)
            : beam.J * fast_exp_im(m_waveIndex * info.projLenght);

        // Extract scalars for inlined computation
        double horAx = info.horAxis.x, horAy = info.horAxis.y, horAz = info.horAxis.z;
        double verAx = info.verAxis.x, verAy = info.verAxis.y, verAz = info.verAxis.z;
        double normx = info.normald.x, normy = info.normald.y, normz = info.normald.z;
        double bdx = beamDirD.x, bdy = beamDirD.y, bdz = beamDirD.z;
        double cenx = info.center.x, ceny = info.center.y, cenz = info.center.z;
        double beam_area = info.area;
        // PolData scalars
        double pNTx = polData.NTd.x, pNTy = polData.NTd.y, pNTz = polData.NTd.z;
        double pNPx = polData.NPd.x, pNPy = polData.NPd.y, pNPz = polData.NPd.z;
        double pnxDTx = polData.nxDT.x, pnxDTy = polData.nxDT.y, pnxDTz = polData.nxDT.z;
        double pnxDPx = polData.nxDP.x, pnxDPy = polData.nxDP.y, pnxDPz = polData.nxDP.z;
        // J_phased elements (double components for inline math)
        complex jp00 = J_phased[0][0], jp01 = J_phased[0][1];
        complex jp10 = J_phased[1][0], jp11 = J_phased[1][1];
        double jp00r=real(jp00),jp00i=imag(jp00),jp01r=real(jp01),jp01i=imag(jp01);
        double jp10r=real(jp10),jp10i=imag(jp10),jp11r=real(jp11),jp11i=imag(jp11);

#ifdef _DEBUG // DEB
        sum += beam.Area();
        for (int i = 0; i <= m_sphere.nAzimuth; ++i)
        {
            for (int j = 0; j <= m_sphere.nZenith; ++j)
            {
                if (/*i == m_sphere.nAzimuth &&*/ j == m_sphere.nZenith)
                    int ddd = 0;
#else

        // Precompute vf in SOA layout for better cache access
        int nAz_total = m_sphere.nAzimuth;
        int nZen_total = m_sphere.nZenith;

        for (int i = 0; i <= nAz_total; ++i)
        {
            int nZen = nZen_total;

            // --- Theta-coefficients: precompute per-phi ---
            double cp = cos_phi_arr[i], sp = sin_phi_arr[i];

            ThetaCoeffs tc;
            if (edgeData.valid) {
                precompute_theta_coeffs(
                    edgeData.x, edgeData.y, edgeData.nVertices,
                    horAx, horAy, horAz, verAx, verAy, verAz,
                    bdx, bdy, bdz, cenx, ceny, cenz,
                    cp, sp, m_waveIndex,
                    pNTx, pNTy, pNTz, pNPx, pNPy, pNPz,
                    pnxDTx, pnxDTy, pnxDTz, pnxDPx, pnxDPy, pnxDPz,
                    tc);
            }

            for (int j = 0; j <= nZen; ++j)
            {
#endif
                complex d00(0,0), d01(0,0), d10(0,0), d11(0,0);

                if (edgeData.valid)
                {
                    double sin_t = sin_theta_arr[j];
                    double cos_t = cos_theta_arr[j];

                    // Construct direction from angles (MBS convention: z = -cos(theta))
                    double dx = sin_t * cp;
                    double dy = sin_t * sp;
                    double dz = -cos_t;

                    // vf from sphere (direction-dependent)
                    Point3d &vf = m_sphere.vf[i][j];
                    double vfx = vf.x, vfy = vf.y, vfz = vf.z;

                    // A, B from theta-coefficients
                    double A = sin_t * tc.a_sin + cos_t * tc.a_cos + tc.a0;
                    double B = sin_t * tc.b_sin + cos_t * tc.b_cos + tc.b0;

                    double absA = fabs(A);
                    double absB = fabs(B);

                    complex fresnel;
                    if (absA < m_eps2 && absB < m_eps2)
                    {
                        fresnel = -m_invComplWave * beam_area;
                    }
                    else
                    {
                        // Vertex phases from theta-coefficients
                        int nv = tc.nv;
                        double vc[32], vs[32];
                        double phases[32];
                        for (int v = 0; v < nv; ++v)
                            phases[v] = sin_t * tc.psin[v] + cos_t * tc.pcos[v] + tc.p0[v];

                        // AVX-512 / AVX2 / scalar sincos
                        int vv = 0;
                        for (; vv + 7 < nv; vv += 8)
                            fast_sincos_8x(&phases[vv], &vs[vv], &vc[vv]);
                        for (; vv + 3 < nv; vv += 4)
                            fast_sincos_4x(&phases[vv], &vs[vv], &vc[vv]);
                        for (; vv < nv; ++vv)
                            fast_sincos(phases[vv], vs[vv], vc[vv]);

                        double sr = 0, si = 0;
                        if (absB > absA)
                        {
                            for (int e = 0; e < nv; ++e)
                            {
                                if (!edgeData.edge_valid_x[e]) continue;
                                int enext = (e + 1 < nv) ? e + 1 : 0;
                                double Ci = A + edgeData.slope_yx[e] * B;
                                double absCi = fabs(Ci);
                                double inv_Ci = (absCi > m_eps1) ? (1.0 / Ci) : 0.0;
                                sr += (vc[enext] - vc[e]) * inv_Ci;
                                si += (vs[enext] - vs[e]) * inv_Ci;
                                if (__builtin_expect(absCi <= m_eps1, 0)) {
                                    double p1x = edgeData.x[e], p2x = edgeData.x[enext];
                                    double tr = -m_wi2*Ci*(p2x*p2x-p1x*p1x)*0.5;
                                    double ti = m_waveIndex*(p2x-p1x);
                                    sr += vc[e]*tr - vs[e]*ti;
                                    si += vc[e]*ti + vs[e]*tr;
                                }
                            }
                            double inv_B = 1.0 / B;
                            sr *= inv_B; si *= inv_B;
                        }
                        else
                        {
                            for (int e = 0; e < nv; ++e)
                            {
                                if (!edgeData.edge_valid_y[e]) continue;
                                int enext = (e + 1 < nv) ? e + 1 : 0;
                                double Ei = A * edgeData.slope_xy[e] + B;
                                double absEi = fabs(Ei);
                                double inv_Ei = (absEi > m_eps1) ? (1.0 / Ei) : 0.0;
                                sr += (vc[enext] - vc[e]) * inv_Ei;
                                si += (vs[enext] - vs[e]) * inv_Ei;
                                if (__builtin_expect(absEi <= m_eps1, 0)) {
                                    double p1y = edgeData.y[e], p2y = edgeData.y[enext];
                                    double tr = -m_wi2*Ei*(p2y*p2y-p1y*p1y)*0.5;
                                    double ti = m_waveIndex*(p2y-p1y);
                                    sr += vc[e]*tr - vs[e]*ti;
                                    si += vc[e]*ti + vs[e]*tr;
                                }
                            }
                            double inv_nA = -1.0 / A;
                            sr *= inv_nA; si *= inv_nA;
                        }

                        fresnel = m_complWave * complex(sr, si);
                    }

                    if (!isnan(real(fresnel)))
                    {
                        double dpr,dpi;
                        if (!isExternal) {
                            double dpArg = -m_waveIndex*(sin_t*tc.dp_sin + cos_t*tc.dp_cos);
                            fast_sincos(dpArg,dpi,dpr);
                        } else { dpr=1.0; dpi=0.0; }

                        double r00, r01, r10, r11;
                        rotate_jones_inline(
                            pNTx, pNTy, pNTz, pNPx, pNPy, pNPz,
                            pnxDTx, pnxDTy, pnxDTz, pnxDPx, pnxDPy, pnxDPz,
                            vfx, vfy, vfz, dx, dy, dz,
                            r00, r01, r10, r11);

                        double fr = real(fresnel), fi = imag(fresnel);
                        double cpr = fr*dpr - fi*dpi;
                        double cpi = fr*dpi + fi*dpr;
                        double sr00r = cpr*r00, sr00i = cpi*r00;
                        double sr01r = cpr*r01, sr01i = cpi*r01;
                        double sr10r = cpr*r10, sr10i = cpi*r10;
                        double sr11r = cpr*r11, sr11i = cpi*r11;
                        d00=complex(sr00r*jp00r-sr00i*jp00i+sr01r*jp10r-sr01i*jp10i,
                                    sr00r*jp00i+sr00i*jp00r+sr01r*jp10i+sr01i*jp10r);
                        d01=complex(sr00r*jp01r-sr00i*jp01i+sr01r*jp11r-sr01i*jp11i,
                                    sr00r*jp01i+sr00i*jp01r+sr01r*jp11i+sr01i*jp11r);
                        d10=complex(sr10r*jp00r-sr10i*jp00i+sr11r*jp10r-sr11i*jp10i,
                                    sr10r*jp00i+sr10i*jp00r+sr11r*jp10i+sr11i*jp10r);
                        d11=complex(sr10r*jp01r-sr10i*jp01i+sr11r*jp11r-sr11i*jp11i,
                                    sr10r*jp01i+sr10i*jp01r+sr11r*jp11i+sr11i*jp11r);

                    }
                }
                else
                {
                    Point3d &dir = m_sphere.directions[i][j];
                    Point3d &vf = m_sphere.vf[i][j];
                    matrixC tmp = ApplyDiffraction(beam, info, dir, vf);
                    d00 = tmp[0][0]; d01 = tmp[0][1];
                    d10 = tmp[1][0]; d11 = tmp[1][1];
                }

                if (!isCoh)
                {
                    matrixC diffractedMatrix(2, 2);
                    diffractedMatrix[0][0] = d00; diffractedMatrix[0][1] = d01;
                    diffractedMatrix[1][0] = d10; diffractedMatrix[1][1] = d11;
                    matrix m = Mueller(diffractedMatrix);
                    m *= m_sinZenith;
                    M.insert(i, j, m);
                }
                else
                {
                    m_diffractedMatrices[groupId].insert_2x2(i, j, d00, d01, d10, d11);
                }
            }
        }
#ifdef _DEBUG // DEB
        complex ddd = m_diffractedMatrices[0](0,0)[0][0];
        int fff = 0;
#endif
    }

    auto beam_loop_end = std::chrono::high_resolution_clock::now();
    g_time_beamLoop_us += std::chrono::duration<double, std::micro>(beam_loop_end - beam_loop_start).count();

    if (isCoh)
    {
        auto atm0 = std::chrono::high_resolution_clock::now();
        AddToMueller();
        auto atm1 = std::chrono::high_resolution_clock::now();
        g_time_addToMueller_us += std::chrono::duration<double, std::micro>(atm1 - atm0).count();
    }

#ifdef _DEBUG // DEB
//    double mm = M(28, 20)[0][0];
    int fff = 0;
#endif

//    if (m_tracks->shouldComputeTracksOnly)
//    {
//        AddToMueller();
//    }

    auto hb_end = std::chrono::high_resolution_clock::now();
    g_time_handleBeams_us += std::chrono::duration<double, std::micro>(hb_end - hb_start).count();
}

void HandlerPO::SetTracks(Tracks *tracks)
{
    m_tracks = tracks;

    for (int i = 0; i < m_tracks->size() + 1; ++i)
    {
        Arr2D tmp = Arr2D(m_sphere.nAzimuth+1, m_sphere.nZenith+1, 4, 4);
        m_groupMatrices.push_back(tmp);
    }
}

// =============================================================================
// PrepareBeams: sequential preprocessing of beams from one orientation.
// Extracts all scalar data needed for the direction loop.
// =============================================================================
void HandlerPO::PrepareBeams(std::vector<Beam> &beams, double sinZenith,
                              PreparedOrientation &out)
{
    out.beams.clear();
    out.sinZenith = sinZenith;

    // Pass 1: compute max |J|² and max area for two-threshold cutoff
    double maxJnorm = 0, maxArea = 0;
    for (Beam &beam : beams)
    {
        if (beam.lastFacetId != __INT_MAX__)
        {
            double jn = beam.J.Norm();  // |J|² (Frobenius norm²)
            double ar = beam.Area();
            if (jn > maxJnorm) maxJnorm = jn;
            if (ar > maxArea) maxArea = ar;
        }
    }
    // Two-threshold cutoff: skip beam only if BOTH |J|²/max AND area/max < eps.
    // Protects: large-area beams (forward peak) and strong narrow beams.
    // Both ratios are dimensionless (0..1), independent of particle size.
    double jThreshold = maxJnorm * m_targetEps;
    double areaThreshold = maxArea * m_targetEps;
    int skippedBeams = 0;

    // Pass 2: prepare beams, skip negligible ones
    for (Beam &beam : beams)
    {
        if (isBackScatteringConusEnabled && beam.direction.cz < backScatteringConus)
            continue;

        beam.polarizationBasis = beam.RotateSpherical(
                    -m_incidentLight->direction,
                    m_incidentLight->polarizationBasis);

        m_isBadBeam = false;
        BeamInfo info = ComputeBeamInfo(beam);

        if (m_isBadBeam)
            continue;

        if (m_hasAbsorption && beam.lastFacetId != __INT_MAX__ && beam.lastFacetId != -1)
            ApplyAbsorption(beam);

        if (beam.lastFacetId != __INT_MAX__)
        {
            matrix m_ = Mueller(beam.J);
            m_outputEnergy += BeamCrossSection(beam)*m_[0][0]*sinZenith;

            // Skip beam only if BOTH |J|² and area are small
            // Protects: large-area beams (forward peak) and strong narrow beams
            double jn = beam.J.Norm();
            double ar = info.area;
            if (jn < jThreshold && ar < areaThreshold)
            {
                skippedBeams++;
                continue;
            }
        }

        PreparedBeam pb;
        pb.info = info;

        // Precompute edge and pol data
        PrecomputeEdgeData(info, beam, pb.edgeData);
        PrecomputePolData(beam, info, pb.polData);

        pb.isExternal = (beam.lastFacetId == __INT_MAX__);
        matrixC J_phased = pb.isExternal
            ? beam.J * fast_exp_im(m_waveIndex * beam.opticalPath)
            : beam.J * fast_exp_im(m_waveIndex * info.projLenght);

        // Extract scalars
        pb.horAx = info.horAxis.x; pb.horAy = info.horAxis.y; pb.horAz = info.horAxis.z;
        pb.verAx = info.verAxis.x; pb.verAy = info.verAxis.y; pb.verAz = info.verAxis.z;
        pb.normx = info.normald.x; pb.normy = info.normald.y; pb.normz = info.normald.z;
        Point3d beamDirD(beam.direction.cx, beam.direction.cy, beam.direction.cz);
        pb.bdx = beamDirD.x; pb.bdy = beamDirD.y; pb.bdz = beamDirD.z;
        pb.cenx = info.center.x; pb.ceny = info.center.y; pb.cenz = info.center.z;
        pb.beam_area = info.area;
        pb.pNTx = pb.polData.NTd.x; pb.pNTy = pb.polData.NTd.y; pb.pNTz = pb.polData.NTd.z;
        pb.pNPx = pb.polData.NPd.x; pb.pNPy = pb.polData.NPd.y; pb.pNPz = pb.polData.NPd.z;
        pb.pnxDTx = pb.polData.nxDT.x; pb.pnxDTy = pb.polData.nxDT.y; pb.pnxDTz = pb.polData.nxDT.z;
        pb.pnxDPx = pb.polData.nxDP.x; pb.pnxDPy = pb.polData.nxDP.y; pb.pnxDPz = pb.polData.nxDP.z;

        complex jp00 = J_phased[0][0], jp01 = J_phased[0][1];
        complex jp10 = J_phased[1][0], jp11 = J_phased[1][1];
        pb.jp00r = real(jp00); pb.jp00i = imag(jp00);
        pb.jp01r = real(jp01); pb.jp01i = imag(jp01);
        pb.jp10r = real(jp10); pb.jp10i = imag(jp10);
        pb.jp11r = real(jp11); pb.jp11i = imag(jp11);

        // Store original beam for fallback path
        pb.origBeam = beam;

        out.beams.push_back(pb);
    }

    // Log beam cutoff statistics (first call only)
    bool logged = false;
    if (!logged && skippedBeams > 0) {
        std::cerr << "Beam cutoff: " << skippedBeams << "/" << (skippedBeams + (int)out.beams.size())
                  << " beams skipped (|J|²<" << std::scientific << std::setprecision(1)
                  << jThreshold << " AND area<" << areaThreshold << ")" << std::endl;
        logged = true;
    }
}

// =============================================================================
// HandleBeamsToLocal: Process prepared beams into a LOCAL Mueller accumulator.
// Thread-safe: only reads immutable handler data (m_sphere, wave constants)
// and writes to the provided localM / localJ arrays.
// =============================================================================
// =============================================================================
// DiffractControlPoints: fast diffraction for nPoints theta values,
// phi-averaged (same as full HandleBeamsToLocal but only at control thetas).
// Cost: nPoints × N_phi × beams — ~100× cheaper than full N_theta × N_phi.
// =============================================================================
void HandlerPO::DiffractControlPoints(const PreparedOrientation &prepared,
                                       const int *thetaIndices, int nPoints,
                                       double *m11_out)
{
    for (int k = 0; k < nPoints; ++k) m11_out[k] = 0;

    int nAz = m_sphere.nAzimuth;

    for (const PreparedBeam &pb : prepared.beams)
    {
        const BeamEdgeData &edgeData = pb.edgeData;
        if (!edgeData.valid) continue;

        double bdx = pb.bdx, bdy = pb.bdy, bdz = pb.bdz;
        double horAx = pb.horAx, horAy = pb.horAy, horAz = pb.horAz;
        double verAx = pb.verAx, verAy = pb.verAy, verAz = pb.verAz;
        double cenx = pb.cenx, ceny = pb.ceny, cenz = pb.cenz;
        double beam_area = pb.beam_area;
        double pNTx = pb.pNTx, pNTy = pb.pNTy, pNTz = pb.pNTz;
        double pNPx = pb.pNPx, pNPy = pb.pNPy, pNPz = pb.pNPz;
        double pnxDTx = pb.pnxDTx, pnxDTy = pb.pnxDTy, pnxDTz = pb.pnxDTz;
        double pnxDPx = pb.pnxDPx, pnxDPy = pb.pnxDPy, pnxDPz = pb.pnxDPz;
        double jp00r = pb.jp00r, jp00i = pb.jp00i;
        double jp01r = pb.jp01r, jp01i = pb.jp01i;
        double jp10r = pb.jp10r, jp10i = pb.jp10i;
        double jp11r = pb.jp11r, jp11i = pb.jp11i;
        bool isExternal = pb.isExternal;
        int nv = edgeData.nVertices;

        // Loop over ALL phi bins (same as HandleBeamsToLocal)
        for (int iPhi = 0; iPhi <= nAz; ++iPhi)
        {
            double cp = cos(iPhi * m_sphere.azinuthStep);
            double sp = sin(iPhi * m_sphere.azinuthStep);

            ThetaCoeffs tc;
            precompute_theta_coeffs(
                edgeData.x, edgeData.y, nv,
                horAx, horAy, horAz, verAx, verAy, verAz,
                bdx, bdy, bdz, cenx, ceny, cenz,
                cp, sp, m_waveIndex,
                pNTx, pNTy, pNTz, pNPx, pNPy, pNPz,
                pnxDTx, pnxDTy, pnxDTz, pnxDPx, pnxDPy, pnxDPz,
                tc);

            // Only nPoints theta values (not full nZenith)
            for (int k = 0; k < nPoints; ++k)
            {
                int j = thetaIndices[k];
                double theta_rad = m_sphere.GetZenith(j);
                double sin_t, cos_t;
                fast_sincos(theta_rad, sin_t, cos_t);

                double dx = sin_t*cp, dy = sin_t*sp, dz = -cos_t;

                double A = sin_t*tc.a_sin + cos_t*tc.a_cos + tc.a0;
                double B = sin_t*tc.b_sin + cos_t*tc.b_cos + tc.b0;
                double absA = fabs(A), absB = fabs(B);

                complex fresnel;
                if (absA < m_eps2 && absB < m_eps2) {
                    fresnel = -m_invComplWave * beam_area;
                } else {
                    double phases[32], vc[32], vs[32];
                    for (int v = 0; v < nv; ++v)
                        phases[v] = sin_t*tc.psin[v] + cos_t*tc.pcos[v] + tc.p0[v];
                    for (int v = 0; v < nv; ++v)
                        fast_sincos(phases[v], vs[v], vc[v]);

                    double sr = 0, si = 0;
                    if (absB > absA) {
                        for (int e = 0; e < nv; ++e) {
                            if (!edgeData.edge_valid_x[e]) continue;
                            int en = (e+1<nv)?e+1:0;
                            double Ci = A + edgeData.slope_yx[e]*B;
                            double inv = (fabs(Ci)>m_eps1)?1.0/Ci:0.0;
                            sr += (vc[en]-vc[e])*inv; si += (vs[en]-vs[e])*inv;
                        }
                        sr /= B; si /= B;
                    } else {
                        for (int e = 0; e < nv; ++e) {
                            if (!edgeData.edge_valid_y[e]) continue;
                            int en = (e+1<nv)?e+1:0;
                            double Ei = A*edgeData.slope_xy[e]+B;
                            double inv = (fabs(Ei)>m_eps1)?1.0/Ei:0.0;
                            sr += (vc[en]-vc[e])*inv; si += (vs[en]-vs[e])*inv;
                        }
                        double inv_nA = -1.0/A; sr *= inv_nA; si *= inv_nA;
                    }
                    fresnel = m_complWave * complex(sr, si);
                }
                if (isnan(real(fresnel))) continue;

                double dpr=1.0, dpi=0.0;
                if (!isExternal) {
                    double dpArg = -m_waveIndex*(sin_t*tc.dp_sin+cos_t*tc.dp_cos);
                    fast_sincos(dpArg, dpi, dpr);
                }

                Point3d &vf = m_sphere.vf[iPhi][j];
                double r00,r01,r10,r11;
                rotate_jones_inline(pNTx,pNTy,pNTz,pNPx,pNPy,pNPz,
                    pnxDTx,pnxDTy,pnxDTz,pnxDPx,pnxDPy,pnxDPz,
                    vf.x,vf.y,vf.z,dx,dy,dz,r00,r01,r10,r11);

                double fr=real(fresnel),fi=imag(fresnel);
                double cpr=fr*dpr-fi*dpi, cpi=fr*dpi+fi*dpr;
                double sr00r=cpr*r00,sr00i=cpi*r00,sr01r=cpr*r01,sr01i=cpi*r01;
                double sr10r=cpr*r10,sr10i=cpi*r10,sr11r=cpr*r11,sr11i=cpi*r11;
                double d00r=sr00r*jp00r-sr00i*jp00i+sr01r*jp10r-sr01i*jp10i;
                double d00i=sr00r*jp00i+sr00i*jp00r+sr01r*jp10i+sr01i*jp10r;
                double d11r=sr10r*jp01r-sr10i*jp01i+sr11r*jp11r-sr11i*jp11i;
                double d11i=sr10r*jp01i+sr10i*jp01r+sr11r*jp11i+sr11i*jp11r;

                // Incoherent M11 (phi-averaged): |d00|²+|d11|²
                m11_out[k] += (d00r*d00r+d00i*d00i + d11r*d11r+d11i*d11i)
                              * prepared.sinZenith / (nAz + 1);
            }
        }
    }
}

void HandlerPO::HandleBeamsToLocal(const PreparedOrientation &prepared,
                                    Arr2D &localM,
                                    std::vector<Arr2DC> &localJ,
                                    std::vector<Arr2DC> *localJ_noshadow)
{
    double sinZenith = prepared.sinZenith;
    int nAz_global = m_sphere.nAzimuth;
    int nZen_global = m_sphere.nZenith;

    // Precompute sin/cos for theta and phi (same as HandleBeams)
    std::vector<double> sin_theta_arr(nZen_global+1), cos_theta_arr(nZen_global+1);
    for (int j = 0; j <= nZen_global; ++j) {
        double theta_rad = m_sphere.GetZenith(j);
        fast_sincos(theta_rad, sin_theta_arr[j], cos_theta_arr[j]);
    }
    std::vector<double> sin_phi_arr(nAz_global+1), cos_phi_arr(nAz_global+1);
    for (int i = 0; i <= nAz_global; ++i) {
        double phi_rad = i * m_sphere.azinuthStep;
        fast_sincos(phi_rad, sin_phi_arr[i], cos_phi_arr[i]);
    }

    for (const PreparedBeam &pb : prepared.beams)
    {
        const BeamEdgeData &edgeData = pb.edgeData;
        double bdx = pb.bdx, bdy = pb.bdy, bdz = pb.bdz;
        double horAx = pb.horAx, horAy = pb.horAy, horAz = pb.horAz;
        double verAx = pb.verAx, verAy = pb.verAy, verAz = pb.verAz;
        double cenx = pb.cenx, ceny = pb.ceny, cenz = pb.cenz;
        double beam_area = pb.beam_area;
        double pNTx = pb.pNTx, pNTy = pb.pNTy, pNTz = pb.pNTz;
        double pNPx = pb.pNPx, pNPy = pb.pNPy, pNPz = pb.pNPz;
        double pnxDTx = pb.pnxDTx, pnxDTy = pb.pnxDTy, pnxDTz = pb.pnxDTz;
        double pnxDPx = pb.pnxDPx, pnxDPy = pb.pnxDPy, pnxDPz = pb.pnxDPz;
        double jp00r = pb.jp00r, jp00i = pb.jp00i;
        double jp01r = pb.jp01r, jp01i = pb.jp01i;
        double jp10r = pb.jp10r, jp10i = pb.jp10i;
        double jp11r = pb.jp11r, jp11i = pb.jp11i;
        bool isExternal = pb.isExternal;

        for (int i = 0; i <= nAz_global; ++i)
        {
            int nZen = nZen_global;
            double cp = cos_phi_arr[i], sp = sin_phi_arr[i];

            ThetaCoeffs tc;
            if (edgeData.valid) {
                precompute_theta_coeffs(
                    edgeData.x, edgeData.y, edgeData.nVertices,
                    horAx, horAy, horAz, verAx, verAy, verAz,
                    bdx, bdy, bdz, cenx, ceny, cenz,
                    cp, sp, m_waveIndex,
                    pNTx, pNTy, pNTz, pNPx, pNPy, pNPz,
                    pnxDTx, pnxDTy, pnxDTz, pnxDPx, pnxDPy, pnxDPz,
                    tc);
            }

            // Pre-batch vertex sincos for ALL thetas at once (AVX-512/AVX2)
            int nv = edgeData.valid ? tc.nv : 0;
            std::vector<double> all_vc((nZen+1)*nv), all_vs((nZen+1)*nv);
            std::vector<double> dir_dpr(nZen+1), dir_dpi(nZen+1);
            if (edgeData.valid && nv > 0) {
                std::vector<double> all_phases((nZen+1)*nv);
                int total = 0;
                for (int j = 0; j <= nZen; ++j) {
                    double sin_t = sin_theta_arr[j], cos_t = cos_theta_arr[j];
                    int base = j * nv;
                    for (int v = 0; v < nv; ++v)
                        all_phases[base + v] = sin_t * tc.psin[v] + cos_t * tc.pcos[v] + tc.p0[v];
                    total = base + nv;
                }
                int pp = 0;
                for (; pp + 7 < total; pp += 8)
                    fast_sincos_8x(&all_phases[pp], &all_vs[pp], &all_vc[pp]);
                for (; pp + 3 < total; pp += 4)
                    fast_sincos_4x(&all_phases[pp], &all_vs[pp], &all_vc[pp]);
                for (; pp < total; ++pp)
                    fast_sincos(all_phases[pp], all_vs[pp], all_vc[pp]);
                if (!isExternal) {
                    std::vector<double> dp_phases(nZen+1);
                    for (int j = 0; j <= nZen; ++j)
                        dp_phases[j] = -m_waveIndex*(sin_theta_arr[j]*tc.dp_sin + cos_theta_arr[j]*tc.dp_cos);
                    int jj = 0;
                    for (; jj + 7 <= nZen; jj += 8)
                        fast_sincos_8x(&dp_phases[jj], &dir_dpi[jj], &dir_dpr[jj]);
                    for (; jj + 3 <= nZen; jj += 4)
                        fast_sincos_4x(&dp_phases[jj], &dir_dpi[jj], &dir_dpr[jj]);
                    for (; jj <= nZen; ++jj)
                        fast_sincos(dp_phases[jj], dir_dpi[jj], dir_dpr[jj]);
                } else {
                    for (int j = 0; j <= nZen; ++j) { dir_dpr[j] = 1.0; dir_dpi[j] = 0.0; }
                }
            }

            for (int j = 0; j <= nZen; ++j)
            {
                complex d00(0,0), d01(0,0), d10(0,0), d11(0,0);

                if (edgeData.valid)
                {
                    double sin_t = sin_theta_arr[j];
                    double cos_t = cos_theta_arr[j];

                    double dx = sin_t * cp;
                    double dy = sin_t * sp;
                    double dz = -cos_t;

                    Point3d &vf = m_sphere.vf[i][j];
                    double vfx = vf.x, vfy = vf.y, vfz = vf.z;

                    double A = sin_t * tc.a_sin + cos_t * tc.a_cos + tc.a0;
                    double B = sin_t * tc.b_sin + cos_t * tc.b_cos + tc.b0;

                    double absA = fabs(A);
                    double absB = fabs(B);

                    complex fresnel;
                    if (absA < m_eps2 && absB < m_eps2)
                    {
                        fresnel = -m_invComplWave * beam_area;
                    }
                    else
                    {
                        int base = j * nv;
                        double *vc = &all_vc[base];
                        double *vs = &all_vs[base];

                        double sr = 0, si = 0;
                        if (absB > absA)
                        {
                            for (int e = 0; e < nv; ++e)
                            {
                                if (!edgeData.edge_valid_x[e]) continue;
                                int enext = (e + 1 < nv) ? e + 1 : 0;
                                double Ci = A + edgeData.slope_yx[e] * B;
                                double absCi = fabs(Ci);
                                double inv_Ci = (absCi > m_eps1) ? (1.0 / Ci) : 0.0;
                                sr += (vc[enext] - vc[e]) * inv_Ci;
                                si += (vs[enext] - vs[e]) * inv_Ci;
                                if (__builtin_expect(absCi <= m_eps1, 0)) {
                                    double p1x = edgeData.x[e], p2x = edgeData.x[enext];
                                    double tr = -m_wi2*Ci*(p2x*p2x-p1x*p1x)*0.5;
                                    double ti = m_waveIndex*(p2x-p1x);
                                    sr += vc[e]*tr - vs[e]*ti;
                                    si += vc[e]*ti + vs[e]*tr;
                                }
                            }
                            double inv_B = 1.0 / B;
                            sr *= inv_B; si *= inv_B;
                        }
                        else
                        {
                            for (int e = 0; e < nv; ++e)
                            {
                                if (!edgeData.edge_valid_y[e]) continue;
                                int enext = (e + 1 < nv) ? e + 1 : 0;
                                double Ei = A * edgeData.slope_xy[e] + B;
                                double absEi = fabs(Ei);
                                double inv_Ei = (absEi > m_eps1) ? (1.0 / Ei) : 0.0;
                                sr += (vc[enext] - vc[e]) * inv_Ei;
                                si += (vs[enext] - vs[e]) * inv_Ei;
                                if (__builtin_expect(absEi <= m_eps1, 0)) {
                                    double p1y = edgeData.y[e], p2y = edgeData.y[enext];
                                    double tr = -m_wi2*Ei*(p2y*p2y-p1y*p1y)*0.5;
                                    double ti = m_waveIndex*(p2y-p1y);
                                    sr += vc[e]*tr - vs[e]*ti;
                                    si += vc[e]*ti + vs[e]*tr;
                                }
                            }
                            double inv_nA = -1.0 / A;
                            sr *= inv_nA; si *= inv_nA;
                        }

                        fresnel = m_complWave * complex(sr, si);
                    }

                    if (!isnan(real(fresnel)))
                    {
                        double dpr = dir_dpr[j], dpi = dir_dpi[j];

                        double r00, r01, r10, r11;
                        rotate_jones_inline(
                            pNTx, pNTy, pNTz, pNPx, pNPy, pNPz,
                            pnxDTx, pnxDTy, pnxDTz, pnxDPx, pnxDPy, pnxDPz,
                            vfx, vfy, vfz, dx, dy, dz,
                            r00, r01, r10, r11);

                        double fr = real(fresnel), fi = imag(fresnel);
                        double cpr = fr*dpr - fi*dpi;
                        double cpi = fr*dpi + fi*dpr;
                        double sr00r = cpr*r00, sr00i = cpi*r00;
                        double sr01r = cpr*r01, sr01i = cpi*r01;
                        double sr10r = cpr*r10, sr10i = cpi*r10;
                        double sr11r = cpr*r11, sr11i = cpi*r11;
                        d00=complex(sr00r*jp00r-sr00i*jp00i+sr01r*jp10r-sr01i*jp10i,
                                    sr00r*jp00i+sr00i*jp00r+sr01r*jp10i+sr01i*jp10r);
                        d01=complex(sr00r*jp01r-sr00i*jp01i+sr01r*jp11r-sr01i*jp11i,
                                    sr00r*jp01i+sr00i*jp01r+sr01r*jp11i+sr01i*jp11r);
                        d10=complex(sr10r*jp00r-sr10i*jp00i+sr11r*jp10r-sr11i*jp10i,
                                    sr10r*jp00i+sr10i*jp00r+sr11r*jp10i+sr11i*jp10r);
                        d11=complex(sr10r*jp01r-sr10i*jp01i+sr11r*jp11r-sr11i*jp11i,
                                    sr10r*jp01i+sr10i*jp01r+sr11r*jp11i+sr11i*jp11r);

                    }
                }
                else
                {
                    // Fallback path for beams without valid edge data
                    Point3d &dir = m_sphere.directions[i][j];
                    Point3d &vf = m_sphere.vf[i][j];
                    matrixC tmp = ApplyDiffraction(pb.origBeam, pb.info, dir, vf);
                    d00 = tmp[0][0]; d01 = tmp[0][1];
                    d10 = tmp[1][0]; d11 = tmp[1][1];
                }

                if (!isCoh)
                {
                    matrixC diffractedMatrix(2, 2);
                    diffractedMatrix[0][0] = d00; diffractedMatrix[0][1] = d01;
                    diffractedMatrix[1][0] = d10; diffractedMatrix[1][1] = d11;
                    matrix m = Mueller(diffractedMatrix);
                    m *= sinZenith;
                    localM.insert(i, j, m);
                }
                else
                {
                    // For coherent mode, accumulate into localJ[0]
                    localJ[0].insert_2x2(i, j, d00, d01, d10, d11);
                    // Also accumulate without shadow beam (for separate output)
                    if (localJ_noshadow && !isExternal)
                        (*localJ_noshadow)[0].insert_2x2(i, j, d00, d01, d10, d11);
                }
            }
        }
    }
}

// =============================================================================
// AddToMuellerLocal: Convert Jones to Mueller for local arrays (static, thread-safe)
// =============================================================================
void HandlerPO::AddToMuellerLocal(const std::vector<Arr2DC> &localJ,
                                   double normIndex, Arr2D &localM,
                                   int nAz, int nZen)
{
    for (size_t q = 0; q < localJ.size(); ++q)
    {
        for (int t = 0; t <= nZen; ++t)
        {
            for (int p = 0; p <= nAz; ++p)
            {
                matrix m = Mueller(localJ[q](p, t));
                m *= normIndex;
                localM.insert(p, t, m);
            }
        }
    }
}

// =============================================================================
// Phase 1: Cache beams from one orientation
// =============================================================================
void HandlerPO::CacheBeams(std::vector<Beam> &beams, double weight,
                            double D_ref, double incomingEnergy,
                            OrientationBeams &out)
{
    out.beams.clear();
    out.weight = weight;
    out.incomingEnergy = incomingEnergy;
    double inv_D = 1.0 / D_ref;
    double inv_D2 = inv_D * inv_D;

    for (Beam &beam : beams)
    {
        if (isBackScatteringConusEnabled && beam.direction.cz < backScatteringConus)
            continue;

        beam.polarizationBasis = beam.RotateSpherical(
                    -m_incidentLight->direction,
                    m_incidentLight->polarizationBasis);

        BeamInfo info = ComputeBeamInfo(beam);

        if (m_isBadBeam)
            continue;

        if (m_hasAbsorption && beam.lastFacetId != __INT_MAX__ && beam.lastFacetId != -1)
            ApplyAbsorption(beam);

        // Precompute edge data
        BeamEdgeData edgeData;
        PrecomputeEdgeData(info, beam, edgeData);

        if (!edgeData.valid)
            continue; // skip beams we can't handle with fast path

        // Build cached beam
        CachedBeam cb;
        cb.nVertices = edgeData.nVertices;
        cb.edgeDataValid = edgeData.valid;

        // Normalize vertices and intercepts by D_ref
        for (int i = 0; i < cb.nVertices; ++i)
        {
            cb.vx_norm[i] = edgeData.x[i] * inv_D;
            cb.vy_norm[i] = edgeData.y[i] * inv_D;
            cb.slope_yx[i] = edgeData.slope_yx[i]; // dimensionless
            cb.slope_xy[i] = edgeData.slope_xy[i]; // dimensionless
            cb.intercept_y_norm[i] = edgeData.intercept_y[i] * inv_D; // unused in fast path
            cb.intercept_x_norm[i] = edgeData.intercept_x[i] * inv_D; // unused in fast path
            cb.edge_valid_x[i] = edgeData.edge_valid_x[i];
            cb.edge_valid_y[i] = edgeData.edge_valid_y[i];
        }

        // Store Jones matrix
        cb.J = beam.J;

        // Normalized lengths
        cb.opticalPath_norm = beam.opticalPath * inv_D;
        cb.projLength_norm = info.projLenght * inv_D;
        cb.area_norm = info.area * inv_D2;

        cb.isExternal = (beam.lastFacetId == __INT_MAX__);

        // Polarization data (precomputed, direction-independent)
        PrecomputePolData(beam, info, cb.polData);

        // Beam direction (unit vector)
        cb.dirx = beam.direction.cx;
        cb.diry = beam.direction.cy;
        cb.dirz = beam.direction.cz;

        // Aperture axes (unit vectors)
        cb.horAx = info.horAxis.x; cb.horAy = info.horAxis.y; cb.horAz = info.horAxis.z;
        cb.verAx = info.verAxis.x; cb.verAy = info.verAxis.y; cb.verAz = info.verAxis.z;
        cb.normx = info.normald.x; cb.normy = info.normald.y; cb.normz = info.normald.z;

        // Normalized center
        cb.cenx_norm = info.center.x * inv_D;
        cb.ceny_norm = info.center.y * inv_D;
        cb.cenz_norm = info.center.z * inv_D;

        out.beams.push_back(cb);
    }
}

// =============================================================================
// Phase 2: Compute Mueller matrices for multiple sizes from cache
// =============================================================================
void HandlerPO::ComputeFromCache(const BeamCache &cache,
                                  const std::vector<double> &x_sizes,
                                  std::vector<Arr2D> &results_M,
                                  std::vector<double> &results_energy)
{
    int nSizes = x_sizes.size();
    int nAz = m_sphere.nAzimuth;
    int nZen = m_sphere.nZenith;

    // Allocate output Mueller arrays
    results_M.clear();
    results_M.reserve(nSizes);
    for (int s = 0; s < nSizes; ++s)
    {
        results_M.push_back(Arr2D(nAz + 1, nZen + 1, 4, 4));
        results_M[s].ClearArr();
    }

    results_energy.clear();
    results_energy.resize(nSizes, 0.0);

    // Precompute per-size constants
    // x = pi*D/lambda, so D = x*lambda/pi
    // waveIndex = 2*pi/lambda
    // k*D = waveIndex * D = 2*pi/lambda * x*lambda/pi = 2*x
    // But we need to work in terms of D/D_ref scaling
    // For cached beam vertices: actual_vertex = vx_norm * D_ref
    // phase = k * (A * actual_vx + B * actual_vy) where A,B depend on direction (not size)
    //       = k * D_ref * (A * vx_norm + B * vy_norm)
    //       = (2*pi/lambda) * (x*lambda/pi) * (D_ref/(x*lambda/pi)) * D_ref * (A*vx_norm+B*vy_norm)
    // Wait, let me think more carefully.
    //
    // The cached vertices are vx_norm = vx_actual / D_ref (where D_ref = x_ref * lambda / pi).
    // For a new size x, the actual particle is scaled by D/D_ref = x/x_ref.
    // The actual vertices become vx_actual_new = vx_norm * D_ref * (D/D_ref) = vx_norm * D
    //
    // k = 2*pi/lambda (fixed)
    // phase_vertex = k * (A * vx_actual_new + B * vy_actual_new)
    //              = k * D * (A * vx_norm + B * vy_norm)
    //              = (2*pi/lambda) * (x*lambda/pi) * (A * vx_norm + B * vy_norm)
    //              = 2*x * (A * vx_norm + B * vy_norm)
    //
    // Similarly:
    // dirPhase = exp(-i*k*D * (dir . center_norm))
    //          = exp(-i * 2*x * (dir . center_norm))
    //
    // J_phased = J * exp(i*k*D * projLength_norm) [internal]
    //          = J * exp(i * 2*x * projLength_norm) [internal]
    // or       = J * exp(i * 2*x * opticalPath_norm) [external]
    //
    // area_new = area_norm * D^2
    //
    // complWave = (-i * lambda) / (2*pi)^2  -- this scales with lambda, but lambda is fixed
    // invComplWave = i / lambda
    // These don't change with x (wavelength is fixed).
    //
    // The edge sum: sum (exp_next - exp_i) * inv_Ci / B_or_negA
    // The Ci = A + slope*B -- A and B are direction-dependent but NOT size-dependent
    //   because A = dot(k_k0, horAxis) where k_k0 = -dir + beamDir
    //   Wait: A is actually the projection of k_k0 onto aperture plane.
    //   In the original code: A,B come from ChangeCoordinateSystem which projects
    //   the difference vector. The key point is that A,B are dimensionless
    //   (they're components of a unit-ish vector difference), NOT multiplied by k.
    //   The phases are k*A*vx + k*B*vy.
    //
    //   But with scaling, the actual A,B don't change (they depend only on directions),
    //   and the actual vertex positions scale with D. So phase = k*(D/D_ref)*D_ref*(A*vx_norm+B*vy_norm)
    //   = k*D*(A*vx_norm+B*vy_norm) = 2*x*(A*vx_norm+B*vy_norm).
    //
    // The final diffraction integral result: complWave * edge_sum
    //   where edge_sum has units of (1/A or 1/B) * exp differences
    //   In the forward direction (A,B->0): result = invComplWave * area
    //   area scales as D^2. So the diffraction amplitude scales as D^2.
    //
    //   More precisely: the edge_sum involves exp(k*D*stuff)/Ci where Ci is
    //   independent of D. The exp's oscillate faster with D but their
    //   differences (exp_next-exp_i) don't have a simple D scaling.
    //   So we must recompute the full diffraction integral for each x.

    // Precompute 2*x for each size
    std::vector<double> two_x(nSizes);
    for (int s = 0; s < nSizes; ++s)
        two_x[s] = 2.0 * x_sizes[s];

    // D values for each size: D = x * lambda / pi
    std::vector<double> D_vals(nSizes);
    for (int s = 0; s < nSizes; ++s)
        D_vals[s] = x_sizes[s] * m_wavelength / M_PI;

    // k*D for each size
    std::vector<double> kD(nSizes);
    for (int s = 0; s < nSizes; ++s)
        kD[s] = m_waveIndex * D_vals[s]; // = 2*x

    // Hoist wave constants outside all loops (they never change)
    double cwr = real(m_complWave), cwi = imag(m_complWave);
    double icwr = real(m_invComplWave), icwi = imag(m_invComplWave);

    // Precompute sin/cos theta and phi ONCE for all beams
    std::vector<double> sin_theta_cache(nZen+1), cos_theta_cache(nZen+1);
    for (int j = 0; j <= nZen; ++j) {
        double theta_rad = m_sphere.GetZenith(j);
        fast_sincos(theta_rad, sin_theta_cache[j], cos_theta_cache[j]);
    }
    std::vector<double> cos_phi_cache(nAz+1), sin_phi_cache(nAz+1);
    for (int i = 0; i <= nAz; ++i) {
        double phi_rad = i * m_sphere.azinuthStep;
        fast_sincos(phi_rad, sin_phi_cache[i], cos_phi_cache[i]);
    }

    // Process all orientations and beams
    for (const auto &orient : cache.orientations)
    {
        double weight = orient.weight;

        for (const auto &cb : orient.beams)
        {
            if (!cb.edgeDataValid)
                continue;

            double pNTx = cb.polData.NTd.x, pNTy = cb.polData.NTd.y, pNTz = cb.polData.NTd.z;
            double pNPx = cb.polData.NPd.x, pNPy = cb.polData.NPd.y, pNPz = cb.polData.NPd.z;
            double pnxDTx = cb.polData.nxDT.x, pnxDTy = cb.polData.nxDT.y, pnxDTz = cb.polData.nxDT.z;
            double pnxDPx = cb.polData.nxDP.x, pnxDPy = cb.polData.nxDP.y, pnxDPz = cb.polData.nxDP.z;

            double bdx = cb.dirx, bdy = cb.diry, bdz = cb.dirz;
            double horAx = cb.horAx, horAy = cb.horAy, horAz = cb.horAz;
            double verAx = cb.verAx, verAy = cb.verAy, verAz = cb.verAz;
            int nv = cb.nVertices;

            // Precompute J_phased for each size
            // For internal: J * exp(i * kD * projLength_norm)
            // For external: J * exp(i * kD * opticalPath_norm)
            // J_phased elements — stack array, no heap allocation
            struct JPhased {
                double jp00r, jp00i, jp01r, jp01i;
                double jp10r, jp10i, jp11r, jp11i;
            };
            JPhased jp_vec[32]; // max 32 sizes on stack

            // Hoist J component extraction out of size loop
            double j11r = real(cb.J.m11), j11i = imag(cb.J.m11);
            double j12r = real(cb.J.m12), j12i = imag(cb.J.m12);
            double j21r = real(cb.J.m21), j21i = imag(cb.J.m21);
            double j22r = real(cb.J.m22), j22i = imag(cb.J.m22);

            // Phase base: same for all sizes, differs only by kD[s] scale
            double phase_base = cb.isExternal ? cb.opticalPath_norm : cb.projLength_norm;

            for (int s = 0; s < nSizes; ++s)
            {
                double phase_arg = kD[s] * phase_base;
                double pe_r, pe_i;
                fast_sincos(phase_arg, pe_i, pe_r);

                // J_phased = J * exp(i*phase_arg), store doubles directly
                jp_vec[s].jp00r = pe_r * j11r - pe_i * j11i;
                jp_vec[s].jp00i = pe_r * j11i + pe_i * j11r;
                jp_vec[s].jp01r = pe_r * j12r - pe_i * j12i;
                jp_vec[s].jp01i = pe_r * j12i + pe_i * j12r;
                jp_vec[s].jp10r = pe_r * j21r - pe_i * j21i;
                jp_vec[s].jp10i = pe_r * j21i + pe_i * j21r;
                jp_vec[s].jp11r = pe_r * j22r - pe_i * j22i;
                jp_vec[s].jp11i = pe_r * j22i + pe_i * j22r;
            }

            // Loop: direction → size (size loop is innermost)
            for (int i = 0; i <= nAz; ++i)
            {
                double cp = cos_phi_cache[i], sp = sin_phi_cache[i];

                // Theta-coefficients: precompute A,B decomposition per phi
                double neg_cp = -cp, neg_sp = -sp;
                double a_sin = neg_cp*horAx + neg_sp*horAy;
                double a_cos = +horAz;
                double a0 = bdx*horAx + bdy*horAy + bdz*horAz;
                double b_sin = neg_cp*verAx + neg_sp*verAy;
                double b_cos = +verAz;
                double b0 = bdx*verAx + bdy*verAy + bdz*verAz;

                // Per-vertex base_phase coefficients (size-independent)
                double bp_sin[32], bp_cos[32], bp_0[32];
                for (int v = 0; v < nv; ++v) {
                    bp_sin[v] = a_sin*cb.vx_norm[v] + b_sin*cb.vy_norm[v];
                    bp_cos[v] = a_cos*cb.vx_norm[v] + b_cos*cb.vy_norm[v];
                    bp_0[v]   = a0   *cb.vx_norm[v] + b0   *cb.vy_norm[v];
                }

                // dirPhase geometric part
                double dp_sin_coeff = cp*cb.cenx_norm + sp*cb.ceny_norm;
                double dp_cos_coeff = -cb.cenz_norm;

                for (int j = 0; j <= nZen; ++j)
                {
                    double sin_t = sin_theta_cache[j], cos_t = cos_theta_cache[j];
                    Point3d &vf = m_sphere.vf[i][j];

                    double dx = sin_t*cp, dy = sin_t*sp, dz = -cos_t;
                    double vfx = vf.x, vfy = vf.y, vfz = vf.z;

                    // A, B from theta-coefficients (2 FMA each)
                    double A = sin_t*a_sin + cos_t*a_cos + a0;
                    double B = sin_t*b_sin + cos_t*b_cos + b0;

                    double absA = fabs(A);
                    double absB = fabs(B);

                    // Compute rotate_jones (size-independent)
                    double r00, r01, r10, r11;
                    rotate_jones_inline(
                        pNTx, pNTy, pNTz, pNPx, pNPy, pNPz,
                        pnxDTx, pnxDTy, pnxDTz, pnxDPx, pnxDPy, pnxDPz,
                        vfx, vfy, vfz, dx, dy, dz,
                        r00, r01, r10, r11);

                    // dirPhase from theta-coefficients
                    double dir_dot_cen = sin_t*dp_sin_coeff + cos_t*dp_cos_coeff;

                    // Precompute inv_Ci for each edge (size-independent)
                    double inv_Ci_arr[32];
                    bool use_B_branch = (absB > absA);

                    if (use_B_branch)
                    {
                        for (int e = 0; e < nv; ++e)
                        {
                            if (!cb.edge_valid_x[e]) { inv_Ci_arr[e] = 0; continue; }
                            double Ci = A + cb.slope_yx[e] * B;
                            double absCi = fabs(Ci);
                            inv_Ci_arr[e] = (absCi > m_eps1) ? (1.0 / Ci) : 0.0;
                        }
                    }
                    else
                    {
                        for (int e = 0; e < nv; ++e)
                        {
                            if (!cb.edge_valid_y[e]) { inv_Ci_arr[e] = 0; continue; }
                            double Ei = A * cb.slope_xy[e] + B;
                            double absEi = fabs(Ei);
                            inv_Ci_arr[e] = (absEi > m_eps1) ? (1.0 / Ei) : 0.0;
                        }
                    }

                    // base_phase from theta-coefficients (2 FMA per vertex instead of full A*vx+B*vy)
                    double base_phase[32];
                    for (int v = 0; v < nv; ++v)
                        base_phase[v] = sin_t*bp_sin[v] + cos_t*bp_cos[v] + bp_0[v];

                    // Precompute edge next-vertex indices (size-independent)
                    int enext_arr[32];
                    for (int e = 0; e < nv; ++e)
                        enext_arr[e] = (e + 1 < nv) ? e + 1 : 0;

                    // Precompute branch divisor (size-independent, hoisted out of size loop)
                    double branch_inv = use_B_branch ? (1.0 / B) : (-1.0 / A);

                    // -------------------------------------------------------
                    // Batch sincos: compute ALL vertex sincos for ALL sizes
                    // before entering the size loop. Uses AVX to batch across
                    // sizes: for each vertex, compute sincos(kD[0..3]*bp)
                    // simultaneously with fast_sincos_4x.
                    // Layout: vc_all[s][v], vs_all[s][v]
                    // -------------------------------------------------------
                    bool is_forward = (absA < m_eps2 && absB < m_eps2);

                    double vc_all[32][32], vs_all[32][32]; // [size][vertex]

                    if (!is_forward)
                    {
                        // For each vertex, batch the nSizes phases via AVX
                        // phase[s] = kD[s] * base_phase[v]
                        for (int v = 0; v < nv; ++v)
                        {
                            double bp = base_phase[v];
                            // Build array of phases across sizes for this vertex
                            double ph_sizes[32];
                            for (int s = 0; s < nSizes; ++s)
                                ph_sizes[s] = kD[s] * bp;

                            double sv_sizes[32], cv_sizes[32];
                            int ss = 0;
                            for (; ss + 7 < nSizes; ss += 8)
                                fast_sincos_8x(&ph_sizes[ss], &sv_sizes[ss], &cv_sizes[ss]);
                            for (; ss + 3 < nSizes; ss += 4)
                                fast_sincos_4x(&ph_sizes[ss], &sv_sizes[ss], &cv_sizes[ss]);
                            for (; ss < nSizes; ++ss)
                                fast_sincos(ph_sizes[ss], sv_sizes[ss], cv_sizes[ss]);

                            for (int s = 0; s < nSizes; ++s)
                            {
                                vc_all[s][v] = cv_sizes[s];
                                vs_all[s][v] = sv_sizes[s];
                            }
                        }
                    }

                    // Precompute dirPhase for all sizes (batch across sizes)
                    double dpr_all[32], dpi_all[32];
                    if (!cb.isExternal)
                    {
                        double dp_phases[32];
                        for (int s = 0; s < nSizes; ++s)
                            dp_phases[s] = -kD[s] * dir_dot_cen;

                        int ss = 0;
                        for (; ss + 7 < nSizes; ss += 8)
                            fast_sincos_8x(&dp_phases[ss], &dpi_all[ss], &dpr_all[ss]);
                        for (; ss + 3 < nSizes; ss += 4)
                            fast_sincos_4x(&dp_phases[ss], &dpi_all[ss], &dpr_all[ss]);
                        for (; ss < nSizes; ++ss)
                            fast_sincos(dp_phases[ss], dpi_all[ss], dpr_all[ss]);
                    }
                    else
                    {
                        for (int s = 0; s < nSizes; ++s)
                        { dpr_all[s] = 1.0; dpi_all[s] = 0.0; }
                    }

                    // Now loop over sizes (the tight inner loop)
                    for (int s = 0; s < nSizes; ++s)
                    {
                        double D_s = D_vals[s];

                        double fr, fi; // fresnel result

                        // Forward direction: area * invComplWave
                        if (is_forward)
                        {
                            double area_s = cb.area_norm * D_s * D_s;
                            fr = -icwr * area_s;
                            fi = -icwi * area_s;
                        }
                        else
                        {
                            // Use precomputed sincos arrays
                            const double *vc = vc_all[s];
                            const double *vs_arr = vs_all[s];

                            // Edge sum
                            double esr = 0, esi = 0;
                            if (use_B_branch)
                            {
                                for (int e = 0; e < nv; ++e)
                                {
                                    if (!cb.edge_valid_x[e]) continue;
                                    int en = enext_arr[e];
                                    double ic = inv_Ci_arr[e];
                                    esr += (vc[en] - vc[e]) * ic;
                                    esi += (vs_arr[en] - vs_arr[e]) * ic;
                                    if (__builtin_expect(ic == 0.0 && cb.edge_valid_x[e], 0)) {
                                        double Ci = A + cb.slope_yx[e] * B;
                                        double p1x = cb.vx_norm[e] * D_s;
                                        double p2x = cb.vx_norm[en] * D_s;
                                        double tr = -m_wi2*Ci*(p2x*p2x-p1x*p1x)*0.5;
                                        double ti = m_waveIndex*(p2x-p1x);
                                        esr += vc[e]*tr - vs_arr[e]*ti;
                                        esi += vc[e]*ti + vs_arr[e]*tr;
                                    }
                                }
                                esr *= branch_inv;
                                esi *= branch_inv;
                            }
                            else
                            {
                                for (int e = 0; e < nv; ++e)
                                {
                                    if (!cb.edge_valid_y[e]) continue;
                                    int en = enext_arr[e];
                                    double ic = inv_Ci_arr[e];
                                    esr += (vc[en] - vc[e]) * ic;
                                    esi += (vs_arr[en] - vs_arr[e]) * ic;
                                    if (__builtin_expect(ic == 0.0 && cb.edge_valid_y[e], 0)) {
                                        double Ei = A * cb.slope_xy[e] + B;
                                        double p1y = cb.vy_norm[e] * D_s;
                                        double p2y = cb.vy_norm[en] * D_s;
                                        double tr = -m_wi2*Ei*(p2y*p2y-p1y*p1y)*0.5;
                                        double ti = m_waveIndex*(p2y-p1y);
                                        esr += vc[e]*tr - vs_arr[e]*ti;
                                        esi += vc[e]*ti + vs_arr[e]*tr;
                                    }
                                }
                                esr *= branch_inv;
                                esi *= branch_inv;
                            }

                            // fresnel = complWave * s
                            fr = cwr*esr - cwi*esi;
                            fi = cwr*esi + cwi*esr;

                            if (isnan(fr)) continue;
                        }

                        // dirPhase (precomputed above)
                        double dpr = dpr_all[s], dpi = dpi_all[s];

                        // combined = fresnel * dirPhase
                        double cpr = fr*dpr - fi*dpi;
                        double cpi = fr*dpi + fi*dpr;

                        // scaled_rot = combined * jones_rot
                        double sr00r = cpr*r00, sr00i = cpi*r00;
                        double sr01r = cpr*r01, sr01i = cpi*r01;
                        double sr10r = cpr*r10, sr10i = cpi*r10;
                        double sr11r = cpr*r11, sr11i = cpi*r11;

                        // result = scaled_rot * J_phased
                        const JPhased &jp = jp_vec[s];
                        double d00r=sr00r*jp.jp00r-sr00i*jp.jp00i+sr01r*jp.jp10r-sr01i*jp.jp10i;
                        double d00i=sr00r*jp.jp00i+sr00i*jp.jp00r+sr01r*jp.jp10i+sr01i*jp.jp10r;
                        double d01r=sr00r*jp.jp01r-sr00i*jp.jp01i+sr01r*jp.jp11r-sr01i*jp.jp11i;
                        double d01i=sr00r*jp.jp01i+sr00i*jp.jp01r+sr01r*jp.jp11i+sr01i*jp.jp11r;
                        double d10r=sr10r*jp.jp00r-sr10i*jp.jp00i+sr11r*jp.jp10r-sr11i*jp.jp10i;
                        double d10i=sr10r*jp.jp00i+sr10i*jp.jp00r+sr11r*jp.jp10i+sr11i*jp.jp10r;
                        double d11r=sr10r*jp.jp01r-sr10i*jp.jp01i+sr11r*jp.jp11r-sr11i*jp.jp11i;
                        double d11i=sr10r*jp.jp01i+sr10i*jp.jp01r+sr11r*jp.jp11i+sr11i*jp.jp11r;

                        // Inline Mueller computation: avoid matrixC/matrix heap allocation
                        double a11 = d00r*d00r+d00i*d00i, a12 = d01r*d01r+d01i*d01i;
                        double a21 = d10r*d10r+d10i*d10i, a22 = d11r*d11r+d11i*d11i;
                        double A1p = a11+a21, A2p = a12+a22;
                        double A1m = a11-a21, A2m = a12-a22;
                        // C1 = d00*conj(d01), C2 = d11*conj(d10)
                        double c1r = d00r*d01r+d00i*d01i, c1i = d00i*d01r-d00r*d01i;
                        double c2r = d11r*d10r+d11i*d10i, c2i = d11i*d10r-d11r*d10i;
                        // C3 = d00*conj(d10), C4 = d11*conj(d01)
                        double c3r = d00r*d10r+d00i*d10i, c3i = d00i*d10r-d00r*d10i;
                        double c4r = d11r*d01r+d11i*d01i, c4i = d11i*d01r-d11r*d01i;
                        // C5 = d00*conj(d11), C6 = d01*conj(d10)
                        double c5r = d00r*d11r+d00i*d11i, c5i = d00i*d11r-d00r*d11i;
                        double c6r = d01r*d10r+d01i*d10i, c6i = d01i*d10r-d01r*d10i;

                        double w2 = weight * 0.5;
                        // Accumulate into results_M[s](i,j) += Mueller * weight
                        results_M[s](i, j, 0, 0) += (A1p+A2p)*w2;
                        results_M[s](i, j, 0, 1) += (A1p-A2p)*w2;
                        results_M[s](i, j, 0, 2) += (-c1r-c2r)*weight;
                        results_M[s](i, j, 0, 3) += (c2i-c1i)*weight;
                        results_M[s](i, j, 1, 0) += (A1m+A2m)*w2;
                        results_M[s](i, j, 1, 1) += (A1m-A2m)*w2;
                        results_M[s](i, j, 1, 2) += (c2r-c1r)*weight;
                        results_M[s](i, j, 1, 3) += (-c1i-c2i)*weight;
                        results_M[s](i, j, 2, 0) += (-c3r-c4r)*weight;
                        results_M[s](i, j, 2, 1) += (c4r-c3r)*weight;
                        results_M[s](i, j, 2, 2) += (c5r+c6r)*weight;
                        results_M[s](i, j, 2, 3) += (c5i-c6i)*weight;
                        results_M[s](i, j, 3, 0) += (c3i-c4i)*weight;
                        results_M[s](i, j, 3, 1) += (c4i+c3i)*weight;
                        results_M[s](i, j, 3, 2) += (-c5i-c6i)*weight;
                        results_M[s](i, j, 3, 3) += (c5r-c6r)*weight;
                    } // end size loop
                } // end zenith
            } // end azimuth
        } // end beam
    } // end orientation

    // Compute energy: scale incoming energy by (D/D_ref)^2
    for (int s = 0; s < nSizes; ++s)
    {
        double scale2 = (D_vals[s] / cache.D_ref) * (D_vals[s] / cache.D_ref);
        double energy = 0;
        for (const auto &orient : cache.orientations)
            energy += orient.incomingEnergy;
        results_energy[s] = energy * scale2;
    }
}
