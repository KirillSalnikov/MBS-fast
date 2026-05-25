#include "HandlerPOTotal.h"

#include "Mueller.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace
{
std::string LogNameForResult(const std::string &destName)
{
    const std::string suffix = "_noshadow";
    if (destName.size() >= suffix.size()
            && destName.compare(destName.size() - suffix.size(), suffix.size(), suffix) == 0)
    {
        return destName.substr(0, destName.size() - suffix.size()) + "_log.txt";
    }
    return destName + "_log.txt";
}

void AppendTextLog(const std::string &destName, const std::string &text)
{
    std::ofstream out(LogNameForResult(destName), std::ios::app);
    if (out.is_open())
        out << text;
}

void ApplyForwardPoleSymmetry(matrix &m)
{
    const double M00 = m[0][0];
    const double M11 = 0.5 * (m[1][1] + m[2][2]);
    const double M33 = m[3][3];
    const double M03 = m[0][3];

    m.Fill(0.0);
    m[0][0] = M00;
    m[1][1] = M11;
    m[2][2] = M11;
    m[3][3] = M33;
    m[0][3] = M03;
    m[3][0] = M03;
}

void ApplyBackwardPoleSymmetry(matrix &m)
{
    const double M00 = m[0][0];
    const double M11 = 0.5 * (m[1][1] - m[2][2]);
    const double M33 = m[3][3];
    const double M03 = m[0][3];

    m.Fill(0.0);
    m[0][0] = M00;
    m[1][1] = M11;
    m[2][2] = -M11;
    m[3][3] = M33;
    m[0][3] = M03;
    m[3][0] = M03;
}
}

HandlerPOTotal::HandlerPOTotal(Particle *particle, Light *incidentLight, int nTheta,
                               double wavelength)
    : HandlerPO(particle, incidentLight, nTheta, wavelength)
{
    betaMueller = new matrix(4, 4);
}

void HandlerPOTotal::WriteMatricesToFile(std::string &destName, double nrg)
{
    std::ofstream outFile(destName + ".dat", std::ios::out);

    if (!outFile.is_open())
    {
        // int fff = 0;
    }

    outFile << std::setprecision(10);
    outFile << "ScAngle 2pi*dcos "\
            "M11 M12 M13 M14 "\
            "M21 M22 M23 M24 "\
            "M31 M32 M33 M34 "\
            "M41 M42 M43 M44";
    matrix Msum(4, 4);

    auto &Lp = *m_Lp;

    int &nZen = m_sphere.nZenith;
    int &nAz = m_sphere.nAzimuth;

    // Accumulate C_sca = integral of M11 * 2pi*dcos over zenith bins
    double C_sca = 0.0;

    for (int iZen = 0; iZen <= nZen; ++iZen)
//    for (int iZen = nZen; iZen >= 0; --iZen)
    {
        Msum.Fill(0.0);
        double radZen = m_sphere.GetZenith(iZen);
        const bool isForwardPole = radZen < __FLT_EPSILON__;
        const bool isBackwardPole = radZen > M_PI-__FLT_EPSILON__;
//        double tt = RadToDeg(m_sphere.zenithEnd) - RadToDeg((t*dT));

        for (int iAz = 0; iAz < nAz; ++iAz)
        {
            double radAz = -iAz*m_sphere.azinuthStep;
            matrix m = M(iAz, iZen);

            Lp[1][1] = cos(2*radAz);
            Lp[1][2] = sin(2*radAz);
            Lp[2][1] = -Lp[1][2];
            Lp[2][2] = Lp[1][1];

            if (isForwardPole || isBackwardPole)
            {
                Msum += m;
            }
            else
            {
                Msum += m*Lp;
            }
        }

        double _2Pi_dcos = m_sphere.Compute2PiDcos(iZen);

        Msum /= m_sphere.nAzimuth;
        if (isBackwardPole)
        {
            ApplyBackwardPoleSymmetry(Msum);
        }
        else if (isForwardPole)
        {
            ApplyForwardPoleSymmetry(Msum);
        }

        // Accumulate C_sca from azimuth-averaged M11
        C_sca += Msum[0][0] * _2Pi_dcos;

        outFile << std::endl << RadToDeg(radZen) << ' ' << _2Pi_dcos << ' ' /*<< nrg << ' '*/;
        outFile << Msum;
    }

    outFile.close();

    // Compute efficiencies. nrg is the orientation-averaged projected area.
    if (nrg > 0)
    {
        const double C_abs_raw = m_hasAbsorption ? (nrg - m_outputEnergy) : 0.0;
        const double absTol = std::max(1.0, nrg) * 1e-10;
        const double C_abs = (std::fabs(C_abs_raw) < absTol) ? 0.0 : C_abs_raw;
        const double C_ext_legacy = C_sca + C_abs;
        const double C_ext = m_hasExtinctionOt
            ? m_extinctionCrossSectionOt : C_ext_legacy;
        const double Q_sca = C_sca / nrg;
        const double Q_abs = C_abs / nrg;
        const double Q_ext = C_ext / nrg;
        const double Q_ext_legacy = C_ext_legacy / nrg;
        const std::string label =
            (destName.find("noshadow") != std::string::npos) ? "no-shadow" : "full";
        if (label == "full")
            m_integralSummary.clear();

        std::ostringstream log;
        log << std::fixed << std::setprecision(4);
        log << "\n===== SCATTERING EFFICIENCY: " << label << " =====\n";
        log << "C_sca = " << C_sca << "\n";
        log << "A_proj (incoming energy) = " << nrg << "\n";
        log << "Q_sca (" << label << ") = C_sca / A_proj = " << Q_sca << "\n";
        if (label == "full")
        {
            log << "Outcoming energy = " << m_outputEnergy << "\n";
            log << "C_abs = A_proj - outcoming energy = " << C_abs << "\n";
            if (m_hasExtinctionOt)
            {
                log << "C_ext_OT = optical theorem forward amplitude = "
                    << m_extinctionCrossSectionOt << "\n";
                log << "Q_ext_OT = C_ext_OT / A_proj = " << Q_ext << "\n";
                log << "C_ext_legacy = C_sca + C_abs_GO = "
                    << C_ext_legacy << "\n";
                log << "Q_ext_legacy = C_ext_legacy / A_proj = "
                    << Q_ext_legacy << "\n";
            }
            else
            {
                log << "C_ext = C_sca + C_abs = " << C_ext << "\n";
            }
            log << "Q_abs = C_abs / A_proj = " << Q_abs << "\n";
            log << "Q_ext = C_ext / A_proj = " << Q_ext << "\n";
            log << "EFFICIENCY_SUMMARY "
                << "Qext=" << Q_ext << ' '
                << "Cext=" << C_ext << ' '
                << "Qext_legacy=" << Q_ext_legacy << ' '
                << "Cext_legacy=" << C_ext_legacy << ' '
                << "Cext_OT=" << (m_hasExtinctionOt
                    ? m_extinctionCrossSectionOt : 0.0) << ' '
                << "Qabs=" << Q_abs << ' '
                << "Qsca=" << Q_sca << ' '
                << "Csca=" << C_sca << ' '
                << "Cabs=" << C_abs << "\n";
        }
        if (Q_sca > 2.5)
        {
            log << "WARNING: Q_sca = " << Q_sca
                << " > 2. PO overestimates scattering at this size parameter.\n"
                << "  Physical limit (extinction paradox): Q_ext -> 2 for large x.\n"
                << "  This is a known PO limitation.\n";
        }
        log << "=========================================\n";
        std::cerr << log.str();
        m_integralSummary += log.str();
        AppendTextLog(destName, log.str());
    }
}

void HandlerPOTotal::AddToMueller()
{
#ifdef _DEBUG // DEB
    double sum = 0;
#endif
    for (size_t q = 0; q < m_diffractedMatrices.size(); ++q)
    {
        auto &diffM = m_diffractedMatrices[q];

        for (int t = 0; t <= m_sphere.nZenith; ++t)
        {
            for (int p = 0; p < m_sphere.nAzimuth; ++p)
            {
                matrix m = Mueller(diffM(p, t));
#ifdef _DEBUG // DEB
                double fff[4];
                fff[0] = m[0][0];
#endif
                if (q == 0 && t == 0 && p == 0)
                {
                    *betaMueller += m*normIndexGamma;
                }

                m *= m_sinZenith;

#ifdef _DEBUG // DEB
                complex ddd[4];
                ddd[0] = diffM(p, t)[0][0];
                double d = m[0][0];
//                if (t == 160)
                {
                    sum += d;
//                    m_logFile << p << ' ' << t << ' ' << sum << std::endl;
                }
#endif
                M.insert(p, t, m);
            }
        }
    }
#ifdef _DEBUG // DEB
    double fffefwe = M(0,0)[0][0];
    int ff = 0;
#endif
}

void HandlerPOTotal::OutputContribution(double angle, double energy)
{
    *(betaFile) << RadToDeg(angle) << ' ' << energy << ' ';
    *(betaFile) << *(betaMueller) << std::endl;
    betaMueller->Fill(0.f);
}
