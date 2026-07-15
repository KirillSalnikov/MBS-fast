#include "HandlerPOTotal.h"

#include "IntegralCharacteristics.h"
#include "Mueller.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace
{
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
        throw std::runtime_error(
            "cannot open Mueller output '" + destName
            + ".dat'.\n  Fix: verify output permissions and free disk space.");

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

    double C_sca_integral = 0.0;
    std::vector<double> thetaRows;
    std::vector<double> m11Rows;
    thetaRows.reserve(nZen + 1);
    m11Rows.reserve(nZen + 1);

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

        C_sca_integral += Msum[0][0] * _2Pi_dcos;
        thetaRows.push_back(radZen);
        m11Rows.push_back(Msum[0][0]);

        outFile << std::endl << RadToDeg(radZen) << ' ' << _2Pi_dcos << ' ' /*<< nrg << ' '*/;
        outFile << Msum;
    }

    outFile.flush();
    if (!outFile)
        throw std::runtime_error(
            "failed while writing Mueller output '" + destName
            + ".dat'.\n  Fix: verify free disk space and filesystem health.");
    outFile.close();

    const std::string label =
        (destName.find("noshadow") != std::string::npos) ? "no-shadow" : "full";
    const double cScaError = EstimateAngularIntegralRelativeError(
        thetaRows, m11Rows, C_sca_integral);
    const IntegralCharacteristics characteristics = ComputeIntegralCharacteristics(
        IntegralMethod::PhysicalOptics, nrg, C_sca_integral,
        m_extinctionCrossSectionOt, m_hasExtinctionOt, m_hasAbsorption,
        cScaError);
    const std::string log = FormatIntegralCharacteristicsLog(
        characteristics, label);
    WriteIntegralCharacteristicsTsv(destName, label, characteristics);
    AppendIntegralCharacteristicsLog(destName, log);
    std::cerr << log;
    if (label == "full")
        m_integralSummary = log;
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
