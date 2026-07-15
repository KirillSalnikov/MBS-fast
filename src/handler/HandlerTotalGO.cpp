#include "HandlerTotalGO.h"

#include "IntegralCharacteristics.h"

#include <iostream>
#include <limits>

HandlerTotalGO::HandlerTotalGO(Particle *particle, Light *incidentLight, int nTheta,
                               float wavelength)
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

        const float &z = acos(beam.direction.cz);
        matrix m = ComputeMueller(RadToDeg(z), beam);

        m_totalContrib.AddMueller(z, m);
    }
}

void HandlerTotalGO::WriteMatricesToFile(std::string &destName, double nrg)
{
    AverageOverAlpha(true, m_normIndex, m_totalContrib, destName);
    WriteToFile(m_totalContrib, m_normIndex, destName + "_all");

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
