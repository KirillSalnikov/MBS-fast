#include "TracerGO.h"
#include "HandlerGO.h"
#include "HandlerTotalGO.h"
#include "Sobol.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <atomic>
#include <exception>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace
{
bool OldautoBetaMidpointEnabledGO()
{
    const char *env = std::getenv("MBS_OLDAUTO_BETA_MIDPOINT");
    if (!env || !*env)
        return true;
    return !(env[0] == '0' && env[1] == '\0');
}

bool OldautoGammaStaggerEnabledGO()
{
    const char *env = std::getenv("MBS_OLDAUTO_GAMMA_STAGGER");
    if (!env || !*env)
        return false;
    return env[0] == '1' && env[1] == '\0';
}

double OldautoGammaAngleGO(const AngleRange &gammaRange, int nGamma,
                           int gammaIndex, int betaIndex, bool stagger)
{
    double unit = gammaIndex + 0.5;
    if (stagger && nGamma > 1)
    {
        const double golden = 0.6180339887498948482;
        double shift = std::fmod((betaIndex + 0.5)*golden, 1.0);
        unit += shift;
        unit -= std::floor(unit/nGamma)*nGamma;
    }
    return gammaRange.min + unit*gammaRange.step;
}

int OldautoBetaCountGO(const AngleRange &betaRange, bool midpoint)
{
    return midpoint ? betaRange.number : betaRange.number + 1;
}

double OldautoBetaAngleGO(const AngleRange &betaRange, int betaIndex,
                          bool midpoint)
{
    return betaRange.min + (betaIndex + (midpoint ? 0.5 : 0.0))
        * betaRange.step;
}

double OldautoTraceBetaGO(double beta, const AngleRange &betaRange)
{
    if (betaRange.step <= 0.0)
        return beta;
    if (fabs(beta) <= FLT_EPSILON)
        return beta + 0.5*betaRange.step;
    if (fabs(beta - M_PI) <= FLT_EPSILON)
        return beta - 0.5*betaRange.step;
    return beta;
}

unsigned int MonteSeedGO()
{
    const char *env = std::getenv("MBS_MONTE_SEED");
    if (env && *env)
    {
        char *end = nullptr;
        unsigned long seed = std::strtoul(env, &end, 10);
        if (end && *end == '\0')
            return (unsigned int)seed;
    }

    unsigned int lo, hi;
    asm("rdtsc" : "=a"(lo), "=d"(hi));
    return lo ^ (hi << 16);
}

int GoThreadCount()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}
}

TracerGO::TracerGO(Particle *particle, int reflNum, const std::string &resultFileName)
	: Tracer(particle, reflNum, resultFileName)
{
}

void TracerGO::TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange)
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif

	vector<Beam> outBeams;
	double beta, gamma;

	CalcTimer timer;
	OutputStartTime(timer);

    const bool betaMidpoint = OldautoBetaMidpointEnabledGO();
    const bool gammaStagger = OldautoGammaStaggerEnabledGO();
    const int nBeta = OldautoBetaCountGO(betaRange, betaMidpoint);
    const int nGamma = gammaRange.number;
    long long orNum = (long long)nBeta * nGamma;
    int betaNorm = (m_symmetry.beta < M_PI_2+FLT_EPSILON && m_symmetry.beta > M_PI_2-FLT_EPSILON) ? 1 : 2;
    double normIndex = gammaRange.number * betaNorm;
    // double norm = CalcNorm(orNum);
    m_handler->SetNormIndex(normIndex);
    if (dynamic_cast<HandlerTotalGO*>(m_handler) != nullptr)
    {
        m_handler->SetNormIndex(1.0);
    }

    HandlerGO *handlerGO = dynamic_cast<HandlerGO*>(m_handler);
    HandlerTotalGO *handlerTotal = dynamic_cast<HandlerTotalGO*>(m_handler);
    const bool parallelTotal = handlerGO && handlerTotal
        && GoThreadCount() > 1 && orNum > 1;

    double cs_beta = 0.0;
    long long count = 0;

    if (parallelTotal)
    {
        m_scattering->PrepareForParallelTrace();
        double incomingEnergySum = 0.0;
        std::atomic<bool> parallelFailed(false);
        std::exception_ptr parallelError;

#pragma omp parallel reduction(+:incomingEnergySum)
        {
            Particle localParticle = *m_particle;
            Scattering *localScatter =
                m_scattering->CloneFor(&localParticle, &m_incidentLight);
            HandlerTotalGO localHandler(&localParticle, &m_incidentLight,
                                        handlerGO->nTheta, m_scattering->m_wave);
            localHandler.ConfigureForThreadLocal(*handlerGO, localScatter);
            localHandler.SetScatteringSphere(handlerGO->m_sphere);
            std::vector<Beam> localBeams;

#pragma omp for schedule(dynamic, 1)
            for (long long idx = 0; idx < orNum; ++idx)
            {
                if (parallelFailed.load(std::memory_order_relaxed))
                    continue;
                try
                {
                const int i = (int)(idx/nGamma);
                const int j = (int)(idx%nGamma);
                const double beta = OldautoBetaAngleGO(betaRange, i,
                                                       betaMidpoint);
                double csLocal = 0.0;
                CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normIndex,
                           csLocal);
                const double gamma = OldautoGammaAngleGO(gammaRange, nGamma, j,
                                                         i, gammaStagger);
                const double traceBeta = OldautoTraceBetaGO(beta, betaRange);
                localParticle.Rotate(traceBeta, gamma, 0);
                localScatter->ScatterLight(0, 0, localBeams);
                localHandler.HandleBeams(localBeams, csLocal);

#ifdef _CHECK_ENERGY_BALANCE
                incomingEnergySum += localScatter->GetIncedentEnergy()*csLocal;
#endif
                const int beamCount = (int)localBeams.size();
                localBeams.clear();

                long long done;
#pragma omp atomic capture
                done = ++count;

#pragma omp critical(go_progress)
                {
                    if (m_logTime == 0)
                    {
                        OutputProgress(orNum, done,
                                       std::lround(RadToDeg(beta)),
                                       std::lround(RadToDeg(gamma)), timer,
                                       beamCount);
                    }
                    else
                    {
                        OutputProgress(orNum, done, i, j, timer, beamCount);
                    }
                }
                }
                catch (...)
                {
                    parallelFailed.store(true, std::memory_order_relaxed);
#pragma omp critical(go_parallel_exception)
                    {
                        if (!parallelError)
                            parallelError = std::current_exception();
                    }
                }
            }

#pragma omp critical(go_merge)
            {
                handlerGO->MergeTotalContributionFrom(localHandler);
            }

            delete localScatter;
        }

#ifdef _CHECK_ENERGY_BALANCE
        m_incomingEnergy += incomingEnergySum;
#endif
        if (parallelError)
            std::rethrow_exception(parallelError);
    }
    else
    {
        for (int i = 0; i < nBeta; ++i)
        {
            beta = OldautoBetaAngleGO(betaRange, i, betaMidpoint);
            CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normIndex, cs_beta);

            for (int j = 0; j < nGamma; ++j)
            {
                gamma = OldautoGammaAngleGO(gammaRange, nGamma, j, i,
                                            gammaStagger);

                const double traceBeta = OldautoTraceBetaGO(beta, betaRange);
                m_particle->Rotate(traceBeta, gamma, 0);
                m_scattering->ScatterLight(0, 0, outBeams);
                m_handler->HandleBeams(outBeams, cs_beta);

#ifdef _CHECK_ENERGY_BALANCE
                m_incomingEnergy += m_scattering->GetIncedentEnergy()*cs_beta;
#endif
                if (m_logTime == 0)
                {
                    OutputProgress(orNum, ++count,
                                   std::lround(RadToDeg(beta)),
                                   std::lround(RadToDeg(gamma)), timer,
                                   outBeams.size());
                }
                else
                {
                    OutputProgress(orNum, ++count, i,j, timer, outBeams.size());
                }

                outBeams.clear();
            }

        }
    }

    // m_incomingEnergy *= normIndex;
    m_handler->m_outputEnergy = ((HandlerGO*)m_handler)->ComputeTotalScatteringEnergy();
    m_handler->WriteMatricesToFile(m_resultDirName, 1000);
    OutputSummary(orNum, timer);
}

void TracerGO::TraceFixed(const double &beta, const double &gamma)
{
	double b = DegToRad(beta);
	double g = DegToRad(gamma);

	vector<Beam> outBeams;
	m_particle->Rotate(b, g, 0);
	m_scattering->ScatterLight(0, 0, outBeams);
    m_handler->HandleBeams(outBeams, 1);
	outBeams.clear();

//	double D_tot = CalcTotalScatteringEnergy();

    m_handler->WriteMatricesToFile(m_resultDirName, 1000);
    //	WriteStatisticsToFileGO(1, D_tot, 1, timer); // TODO: раскомментить
}

void TracerGO::TraceMonteCarlo(const AngleRange &betaRange, const AngleRange &gammaRange, int nOrientations)
{
#ifdef _CHECK_ENERGY_BALANCE
    m_incomingEnergy = 0;
    m_outcomingEnergy = 0;
#endif
    vector<Beam> outBeams;
    double beta, gamma;

    CalcTimer timer;
    OutputStartTime(timer);

    string fulldir = m_resultDirName + "/";

    unsigned int seed = MonteSeedGO();
    srand(seed);
    cout << endl << "Monte seed = " << seed << endl;
    long long count = 0;
    const double cosMin = cos(betaRange.max);
    const double dCos = 1.0 - cosMin;
    std::vector<double> betas(nOrientations);
    std::vector<double> gammas(nOrientations);

    for (int i = 0; i < nOrientations; ++i)
    {
        betas[i] = acos(1.0 - RandomDouble(0, 1)*dCos);
        gammas[i] = RandomDouble(0, 1)*gammaRange.max;
    }

    HandlerGO *handlerGO = dynamic_cast<HandlerGO*>(m_handler);
    HandlerTotalGO *handlerTotal = dynamic_cast<HandlerTotalGO*>(m_handler);
    const bool parallelTotal = handlerGO && handlerTotal
        && GoThreadCount() > 1 && nOrientations > 1;

    if (parallelTotal)
    {
        m_scattering->PrepareForParallelTrace();
        double incomingEnergySum = 0.0;
        std::atomic<bool> parallelFailed(false);
        std::exception_ptr parallelError;

#pragma omp parallel reduction(+:incomingEnergySum)
        {
            Particle localParticle = *m_particle;
            Scattering *localScatter =
                m_scattering->CloneFor(&localParticle, &m_incidentLight);
            HandlerTotalGO localHandler(&localParticle, &m_incidentLight,
                                        handlerGO->nTheta, m_scattering->m_wave);
            localHandler.ConfigureForThreadLocal(*handlerGO, localScatter);
            localHandler.SetScatteringSphere(handlerGO->m_sphere);
            std::vector<Beam> localBeams;

#pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < nOrientations; ++i)
            {
                if (parallelFailed.load(std::memory_order_relaxed))
                    continue;
                try
                {
                    localParticle.Rotate(betas[i], gammas[i], 0);
                    localScatter->ScatterLight(0, 0, localBeams);
                    localHandler.HandleBeams(localBeams, 1.0);
#ifdef _CHECK_ENERGY_BALANCE
                    incomingEnergySum += localScatter->GetIncedentEnergy();
#endif
                }
                catch (...)
                {
                    parallelFailed.store(true, std::memory_order_relaxed);
#pragma omp critical(go_parallel_exception)
                    {
                        if (!parallelError)
                            parallelError = std::current_exception();
                    }
                }

                const int beamCount = (int)localBeams.size();
                localBeams.clear();

                long long done;
#pragma omp atomic capture
                done = ++count;

#pragma omp critical(go_progress)
                {
                    OutputProgress(nOrientations, done, i, i, timer,
                                   beamCount);
                }
            }

#pragma omp critical(go_merge)
            {
                handlerGO->MergeTotalContributionFrom(localHandler);
            }

            delete localScatter;
        }

#ifdef _CHECK_ENERGY_BALANCE
        m_incomingEnergy += incomingEnergySum;
#endif
        if (parallelError)
            std::rethrow_exception(parallelError);
    }
    else
    {
        for (int i = 0; i < nOrientations; ++i)
        {
            beta = betas[i];
            gamma = gammas[i];

            m_particle->Rotate(beta, gamma, 0);
            m_scattering->ScatterLight(0, 0, outBeams);
            m_handler->HandleBeams(outBeams, 1.0);
#ifdef _CHECK_ENERGY_BALANCE
            m_incomingEnergy += m_scattering->GetIncedentEnergy();
#endif

            const int beamCount = (int)outBeams.size();
            outBeams.clear();

            OutputProgress(nOrientations, ++count, i, i, timer, beamCount);
        }
    }

    double norm = 1.0/nOrientations;
    m_handler->SetNormIndex(norm);

    m_outcomingEnergy = ((HandlerGO*)m_handler)->ComputeTotalScatteringEnergy()*norm;
    m_handler->m_outputEnergy = m_outcomingEnergy;
#ifdef _CHECK_ENERGY_BALANCE
    m_incomingEnergy *= norm;
#endif

    // m_resultDirName already has full path (dir/basename) from main.cpp

    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
    OutputSummary(nOrientations, timer);

}

void TracerGO::TraceSobol(int nOrientations, unsigned int seed,
                          double betaMax, double gammaMax)
{
#ifdef _CHECK_ENERGY_BALANCE
    m_incomingEnergy = 0;
    m_outcomingEnergy = 0;
#endif
    vector<Beam> outBeams;

    CalcTimer timer;
    OutputStartTime(timer);

    Sobol2D sobol(seed);
    std::vector<double> su, sv;
    sobol.generate(nOrientations, su, sv);

    const double cosMin = cos(betaMax);
    const double dCos = 1.0 - cosMin;

    long long count = 0;
    std::vector<double> betas(nOrientations);
    std::vector<double> gammas(nOrientations);
    for (int i = 0; i < nOrientations; ++i)
    {
        betas[i] = acos(1.0 - su[i]*dCos);
        gammas[i] = sv[i]*gammaMax;
    }

    HandlerGO *handlerGO = dynamic_cast<HandlerGO*>(m_handler);
    HandlerTotalGO *handlerTotal = dynamic_cast<HandlerTotalGO*>(m_handler);
    const bool parallelTotal = handlerGO && handlerTotal
        && GoThreadCount() > 1 && nOrientations > 1;

    if (parallelTotal)
    {
        m_scattering->PrepareForParallelTrace();
        double incomingEnergySum = 0.0;
        std::atomic<bool> parallelFailed(false);
        std::exception_ptr parallelError;

#pragma omp parallel reduction(+:incomingEnergySum)
        {
            Particle localParticle = *m_particle;
            Scattering *localScatter =
                m_scattering->CloneFor(&localParticle, &m_incidentLight);
            HandlerTotalGO localHandler(&localParticle, &m_incidentLight,
                                        handlerGO->nTheta, m_scattering->m_wave);
            localHandler.ConfigureForThreadLocal(*handlerGO, localScatter);
            localHandler.SetScatteringSphere(handlerGO->m_sphere);
            std::vector<Beam> localBeams;

#pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < nOrientations; ++i)
            {
                if (parallelFailed.load(std::memory_order_relaxed))
                    continue;
                try
                {
                    localParticle.Rotate(betas[i], gammas[i], 0);
                    localScatter->ScatterLight(0, 0, localBeams);
                    localHandler.HandleBeams(localBeams, 1.0);
#ifdef _CHECK_ENERGY_BALANCE
                    incomingEnergySum += localScatter->GetIncedentEnergy();
#endif
                }
                catch (...)
                {
                    parallelFailed.store(true, std::memory_order_relaxed);
#pragma omp critical(go_parallel_exception)
                    {
                        if (!parallelError)
                            parallelError = std::current_exception();
                    }
                }

                const int beamCount = (int)localBeams.size();
                localBeams.clear();

                long long done;
#pragma omp atomic capture
                done = ++count;

#pragma omp critical(go_progress)
                {
                    OutputProgress(nOrientations, done, i, i, timer,
                                   beamCount);
                }
            }

#pragma omp critical(go_merge)
            {
                handlerGO->MergeTotalContributionFrom(localHandler);
            }

            delete localScatter;
        }

#ifdef _CHECK_ENERGY_BALANCE
        m_incomingEnergy += incomingEnergySum;
#endif
        if (parallelError)
            std::rethrow_exception(parallelError);
    }
    else
    {
        for (int i = 0; i < nOrientations; ++i)
        {
            const double beta = betas[i];
            const double gamma = gammas[i];

            m_particle->Rotate(beta, gamma, 0);
            m_scattering->ScatterLight(0, 0, outBeams);
            m_handler->HandleBeams(outBeams, 1.0);
#ifdef _CHECK_ENERGY_BALANCE
            m_incomingEnergy += m_scattering->GetIncedentEnergy();
#endif

            const int beamCount = (int)outBeams.size();
            outBeams.clear();
            OutputProgress(nOrientations, ++count, i, i, timer, beamCount);
        }
    }

    const double norm = 1.0/nOrientations;
    m_handler->SetNormIndex(norm);

    m_outcomingEnergy = ((HandlerGO*)m_handler)->ComputeTotalScatteringEnergy()*norm;
    m_handler->m_outputEnergy = m_outcomingEnergy;
#ifdef _CHECK_ENERGY_BALANCE
    m_incomingEnergy *= norm;
#endif

    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
    OutputSummary(nOrientations, timer);
}

double TracerGO::CalcNorm(long long orNum)
{
	double &symBeta = m_symmetry.beta;
	double dBeta = -(cos(symBeta) - cos(0));
	return symBeta/(orNum*dBeta);
}

void TracerGO::OutputSummary(int orNumber, CalcTimer &timer)
{
	string startTime = ctime(&m_startTime);
	string totalTime = timer.Elapsed();
	time_t end = timer.Stop();
	string endTime = ctime(&end);

	m_summary += "\nStart of calculation = " + startTime
			+ "End of calculation   = " + endTime
			+ "\nTotal time of calculation = " + totalTime
			+ "\nTotal number of body orientation = " + to_string(orNumber)
            /*+ "\nTotal scattering energy = " + to_string(D_tot)*/;

#ifdef _CHECK_ENERGY_BALANCE
    double passedEnergy = (m_handler->m_outputEnergy/m_incomingEnergy)*100;

    m_summary += "\nTotal incoming energy = " + to_string(m_incomingEnergy)
            + "\nTotal outcoming energy = " + to_string(m_handler->m_outputEnergy)
                 + " (S/4 = " + to_string(m_particle->Area()/4)
            + ")\nEnergy passed = " + to_string(passedEnergy) + '%';
#endif

	// out << "\nAveraged cross section = " << incomingEnergy*NRM;
	ofstream out(m_resultDirName + "_out.txt", ios::out);
	if (!out.is_open())
	{
		throw std::runtime_error(
			"cannot open GO summary output.\n"
			"  Fix: verify output permissions and free disk space.");
	}
	out << m_summary;
	out.flush();
	if (!out)
	{
		throw std::runtime_error(
			"failed while writing GO summary output.\n"
			"  Fix: verify free disk space and filesystem health.");
	}

	cout << m_summary;
}
