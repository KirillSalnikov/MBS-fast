#include "HandlerTotalGO.h"

#include "IntegralCharacteristics.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <fstream>
#include <iomanip>
#include <sstream>

HandlerTotalGO::HandlerTotalGO(Particle *particle, Light *incidentLight, int nTheta,
                               double wavelength)
    : HandlerGO(particle, incidentLight, nTheta, wavelength)
{
}

void HandlerTotalGO::HandleBeams(std::vector<Beam> &beams, double sinZenith)
{
    m_sinZenith = sinZenith;

    for (Beam &beam : beams)
    {
#ifdef _DEBUG // DEB
//		m_logFile << beam.id << std::endl;
#endif
        beam.RotateSpherical(-m_incidentLight->direction,
                             m_incidentLight->polarizationBasis);
        // absorbtion
        if (m_hasAbsorption && beam.nActs > 0)
        {
            ApplyAbsorption(beam);
        }

        const double z = acos(std::max(-1.0, std::min(
            1.0, static_cast<double>(beam.direction.cz))));
        matrix m = ComputeMueller(RadToDeg(z), beam);

        m_totalContrib.AddMueller(z, m);
    }
}

void HandlerTotalGO::WriteMatricesToFile(std::string &destName, double nrg)
{
    AverageOverAlpha(true, m_normIndex, m_totalContrib, destName);
    WriteToFile(m_totalContrib, m_normIndex, destName + "_all");

    if(fullautoOutput){
        // Internal fullauto product: finite-bin densities followed by two
        // integrated delta matrices. Column 2 is zero only for delta rows.
        std::ifstream bins((destName+"_all.dat").c_str());
        std::ofstream combined((destName+".dat").c_str());
        if(!bins||!combined)throw std::runtime_error("cannot write GO fullauto product");
        combined<<bins.rdbuf();combined<<std::setprecision(17);
        auto pole=[&](const matrix& raw,bool forward){
            matrix mean(4,4);mean.Fill(0);
            mean[0][0]=raw[0][0];mean[0][3]=raw[0][3];mean[3][0]=raw[3][0];mean[3][3]=raw[3][3];
            const double diagonal=.5*(raw[1][1]+(forward?raw[2][2]:-raw[2][2]));
            const double cross=.5*(raw[1][2]+(forward?-raw[2][1]:raw[2][1]));
            mean[1][1]=diagonal;mean[2][2]=forward?diagonal:-diagonal;
            mean[1][2]=cross;mean[2][1]=forward?-cross:cross;
            combined<<'\n'<<(forward?0:180)<<" 0";
            for(int i=0;i<4;++i)for(int j=0;j<4;++j)combined<<' '<<mean[i][j]*m_normIndex;
        };
        pole(m_totalContrib.forward,true);pole(m_totalContrib.back,false);combined<<'\n';combined.flush();
        if(!combined)throw std::runtime_error("GO fullauto product write failed");
    }

    const IntegralCharacteristics characteristics = ComputeIntegralCharacteristics(
        IntegralMethod::GeometricalOptics, nrg, m_outputEnergy, 0.0,
        false, m_hasAbsorption, std::numeric_limits<double>::quiet_NaN());
    const std::string log = FormatIntegralCharacteristicsLog(
        characteristics, "full");
    WriteIntegralCharacteristicsTsv(destName, "full", characteristics);
    AppendIntegralCharacteristicsLog(destName, log);
    std::cerr << log;
    m_integralSummary = log;
}

void HandlerTotalGO::SetScatteringSphere(const ScatteringRange &grid)
{
    m_sphere = grid;
    m_totalContrib.SetStep(nTheta, grid.zenithEnd - grid.zenithStart);

    for (auto contr : m_tracksContrib)
    {
        contr.SetStep(nTheta, grid.zenithEnd - grid.zenithStart);
    };
}
