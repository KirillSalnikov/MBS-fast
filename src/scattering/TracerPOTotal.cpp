#include "TracerPOTotal.h"
#include "HandlerPOTotal.h"
#include "HandlerPO.h"
#include "BeamCache.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <chrono>
#include <iomanip>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

TracerPOTotal::TracerPOTotal(Particle *particle, int nActs,
                             const string &resultFileName)
    : TracerPO(particle, nActs, resultFileName)
{
}

void TracerPOTotal::TraceRandom(const AngleRange &betaRange,
                                const AngleRange &gammaRange)
{
#ifdef _CHECK_ENERGY_BALANCE
    m_incomingEnergy = 0;
    m_outcomingEnergy = 0;
#endif

    CalcTimer timer;
    long long count = 0;
    long long nOrientations = (betaRange.number) * (gammaRange.number);

#ifdef _DEBUG  /* DEB */
    ofstream outFile(m_resultDirName + "log.dat", ios::out);
    if (!outFile.is_open())
    {
        std::cerr << "Error! File \"" << m_resultDirName
                  << "\" was not opened. " << __FUNCTION__;
        throw std::exception();
    }
#endif


    vector<Beam> outBeams;
    double beta, gamma;

    int betaNorm = (m_symmetry.beta < M_PI_2+FLT_EPSILON && m_symmetry.beta > M_PI_2-FLT_EPSILON) ? 1 : 2;
    double normIndex = gammaRange.number * betaNorm;
    m_handler->SetNormIndex(normIndex);

    std::string dir = CreateFolder(m_resultDirName);
#ifdef _WIN32
    m_resultDirName += '\\' + m_resultDirName;
#else
    m_resultDirName = dir + m_resultDirName;
#endif

    m_handler->betaFile = new std::ofstream(m_resultDirName + "_beta.dat", ios::out);

    *(m_handler->betaFile) << "Beta CS M11 M12 M13 M14 "\
        "M21 M22 M23 M24 "\
        "M31 M32 M33 M34 "\
        "M41 M42 M43 M44";

//    ++nOrientations;
    timer.Start();
    OutputStartTime(timer);

#ifdef _DEBUG  /* DEB */
    for (int i = 0; i <= betaRange.number; ++i)
    {
//        std::cout  << "i: "<< i << std::endl;
#else
    for (int i = 0; i <= betaRange.number; ++i)
    {
#endif
        beta = betaRange.min + i*betaRange.step;

//		double dcos = (i == 0 || i == betaRange.number)
//				? (1.0-cos(0.5*betaRange.step))/normIndex
//				: (cos((i-0.5)*betaRange.step) -
//				   cos((i+0.5)*betaRange.step))/normIndex;

        double dcos;
        CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normIndex, dcos);
        m_handler->SetSinZenith(dcos);

#ifdef _DEBUG // DEB
        for (int j = 0; j < gammaRange.number; ++j)
        {
//            std::cout  << "j: "<< j << std::endl;
//			beta = DegToRad(15); gamma = DegToRad(0);
#else
        for (int j = 0; j < gammaRange.number; ++j)
        {
#endif
            gamma = gammaRange.min + j*gammaRange.step;
            m_particle->Rotate(/*M_PI-*/beta, /*M_PI+*/gamma, 0);

            if (!shadowOff)
            {
                m_scattering->FormShadowBeam(outBeams);
            }

#ifdef _DEBUG  /* DEB */
            // if (i == 4 && j ==0)
            //     int ffff = 0;
#endif
            bool ok = m_scattering->ScatterLight(0, 0, outBeams);

            if (ok)
            {
#ifdef _DEBUG  /* DEB */
//             vector<Beam> be = outBeams;
//            for (int k = 1; k < be.size(); ++k) {
//                if (be[k].nActs==2) {
//                    outBeams.push_back(be[k]);
//                }
//            }
#endif
                m_handler->HandleBeams(outBeams, dcos);
//#ifdef _DEBUG  /* DEB */
//            double sum = 0;
//            for (int k = 0; k < outBeams.size(); ++k) {
//                double aaa = outBeams[k].Area();
////                ofstream outFile(m_resultDirName + to_string(k) + "_vertices.dat", ios::out);
//                sum += aaa;
//                std::vector<int> tr;
//                m_handler->m_tracks->RecoverTrack(outBeams[k], m_particle->nFacets, tr);
//                outFile << outBeams[k].id << " ";
//                for (int l = 0; l < tr.size(); ++l) {
//                    outFile << tr[l] << " ";
//                }
//               outFile << aaa << std::endl;
////                outFile << outBeams[k].arr[0] << endl << endl;
////                outFile.close();
//            }
//            outFile.close();
//            double mm = ((HandlerPO*)m_handler)->M(0, 0)[0][0];
//            std::cout << mm << " " << j << std::endl;
//#else
//#endif
            }
            else
            {
                std::cout << std::endl << "Orientation (" << i << ", " << j << ") has been skipped!!!" << std::endl;
            }

            m_incomingEnergy += m_scattering->GetIncedentEnergy()*dcos;

            OutputProgress(nOrientations, count, i, j, timer, outBeams.size());
            outBeams.clear();
//            std::cout << "sdfewtwr";
            ++count;
        }

        static_cast<HandlerPOTotal*>(m_handler)->OutputContribution(beta, m_incomingEnergy);
    }

    static_cast<HandlerPOTotal*>(m_handler)->betaFile->close();

//#ifdef _DEBUG  /* DEB */
//    double mm = ((HandlerPO*)m_handler)->M(0, 0)[0][0];
//    outFile.close();
//#endif

    EraseConsoleLine(60);
    std::cout << "100%" << std::endl;

    // m_handler->m_outputEnergy = m_handler->ComputeTotalScatteringEnergy();
    m_handler->WriteTotalMatricesToFile(m_resultDirName);
    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
//#ifndef _DEBUG

    OutputStatisticsPO(timer, nOrientations, m_resultDirName);
//#endif
}

void TracerPOTotal::TraceMonteCarlo(const AngleRange &betaRange,
                                    const AngleRange &gammaRange,
                                    int nOrientations)
{
    CalcTimer timer;
    long long count = 0;

    ofstream outFile(m_resultDirName + ".dat", ios::out);

    if (!outFile.is_open())
    {
        std::cerr << "Error! File \"" << m_resultDirName << "\" was not opened. "
                  << __FUNCTION__;

        throw std::exception();
    }

    vector<Beam> outBeams;
    double beta, gamma;
    timer.Start();

    long long nTacts;
    asm("rdtsc" : "=A"(nTacts));
    srand(nTacts);
//    srand(static_cast<unsigned>(time(0)));

    for (int i = 0; i < nOrientations; ++i)
    {
        beta = RandomDouble(0, 1)*betaRange.max;
        gamma = RandomDouble(0, 1)*gammaRange.max;

        m_particle->Rotate(beta, gamma, 0);
        m_scattering->ScatterLight(beta, gamma, outBeams);

        m_handler->HandleBeams(outBeams, sin(beta));
        outBeams.clear();

        ++count;
        OutputProgress(nOrientations, count, i, i, timer, outBeams.size());
    }

    m_handler->WriteTotalMatricesToFile(m_resultDirName);

    std::string dir = CreateFolder(m_resultDirName);
    m_resultDirName = dir + m_resultDirName + '\\' + m_resultDirName;
    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
    OutputStatisticsPO(timer, nOrientations, dir);
    outFile.close();
}

void TracerPOTotal::TraceFromFile(const std::string &orientFile)
{
    // Read orientations from file
    std::ifstream inFile(orientFile);
    if (!inFile.is_open())
    {
        std::cerr << "Error! Cannot open orientation file: " << orientFile << std::endl;
        throw std::exception();
    }

    std::vector<std::pair<double,double>> orientations;
    std::string line;
    while (std::getline(inFile, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        double b, g;
        if (iss >> b >> g)
            orientations.push_back({b, g});
    }
    inFile.close();

    int nOrientations = orientations.size();
    if (nOrientations == 0)
    {
        std::cerr << "Error! No orientations in file: " << orientFile << std::endl;
        throw std::exception();
    }

    CalcTimer timer;
    timer.Start();
    OutputStartTime(timer);

    double weight = 1.0 / nOrientations;
    m_handler->SetNormIndex(1);

    std::string dir = CreateFolder(m_resultDirName);
#ifdef _WIN32
    m_resultDirName += '\\' + m_resultDirName;
#else
    m_resultDirName = dir + m_resultDirName;
#endif

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO)
    {
        std::cerr << "Error! Handler is not HandlerPO in TraceFromFile" << std::endl;
        throw std::exception();
    }

    // =========================================================================
    // Phase 1 (sequential): Trace beams for all orientations, preprocess them
    // =========================================================================
    auto t_phase1_start = std::chrono::high_resolution_clock::now();

    std::vector<PreparedOrientation> allPrepared(nOrientations);
    std::vector<Beam> outBeams;
    long long count = 0;

    for (int i = 0; i < nOrientations; ++i)
    {
        m_particle->Rotate(orientations[i].first, orientations[i].second, 0);

        if (!shadowOff)
            m_scattering->FormShadowBeam(outBeams);

        bool ok = m_scattering->ScatterLight(0, 0, outBeams);

        if (ok)
            handlerPO->PrepareBeams(outBeams, weight, allPrepared[i]);
        else
            allPrepared[i].sinZenith = weight;

        m_incomingEnergy += m_scattering->GetIncedentEnergy() * weight;
        OutputProgress(nOrientations, count, i, 0, timer, outBeams.size());
        outBeams.clear();
        ++count;
    }

    auto t_phase1_end = std::chrono::high_resolution_clock::now();
    double phase1_sec = std::chrono::duration<double>(t_phase1_end - t_phase1_start).count();

    EraseConsoleLine(60);
    std::cout << "Phase 1 (tracing + preprocessing): " << std::fixed
              << std::setprecision(2) << phase1_sec << " s" << std::endl;

    // =========================================================================
    // Phase 2 (parallel): Process beams into Mueller matrices using OpenMP
    // Each thread accumulates into its own private Mueller array, then reduce.
    // =========================================================================
    auto t_phase2_start = std::chrono::high_resolution_clock::now();

    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    #pragma omp parallel
    {
        // Each thread gets its own zero-initialized Mueller accumulator
        Arr2D localM(nAz + 1, nZen + 1, 4, 4);
        localM.ClearArr();

        // For coherent mode, each thread also needs local Jones matrices
        std::vector<Arr2DC> localJ;
        if (handlerPO->isCoh)
        {
            Arr2DC tmp(nAz + 1, nZen + 1, 2, 2);
            tmp.ClearArr();
            localJ.push_back(tmp);
        }

        #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < nOrientations; ++i)
        {
            if (!allPrepared[i].beams.empty())
            {
                handlerPO->HandleBeamsToLocal(allPrepared[i], localM, localJ);
            }

            // For coherent mode: convert Jones->Mueller after each orientation's beams
            if (handlerPO->isCoh && !localJ.empty())
            {
                HandlerPO::AddToMuellerLocal(localJ, 1.0, localM, nAz, nZen);
                // Clear Jones for next orientation
                localJ[0].ClearArr();
            }
        }

        // Reduce: merge thread-local M into handler's global M
        #pragma omp critical
        {
            for (int p = 0; p <= nAz; ++p)
            {
                for (int t = 0; t <= nZen; ++t)
                {
                    matrix m = localM(p, t);
                    handlerPO->M.insert(p, t, m);
                }
            }
        }
    } // end omp parallel

    auto t_phase2_end = std::chrono::high_resolution_clock::now();
    double phase2_sec = std::chrono::duration<double>(t_phase2_end - t_phase2_start).count();

    std::cout << "Phase 2 (diffraction, OpenMP): " << std::fixed
              << std::setprecision(2) << phase2_sec << " s" << std::endl;
    std::cout << "Total: " << phase1_sec + phase2_sec << " s" << std::endl;

    // Free prepared beam memory before writing results
    allPrepared.clear();
    allPrepared.shrink_to_fit();

    m_handler->WriteTotalMatricesToFile(m_resultDirName);
    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
    OutputStatisticsPO(timer, nOrientations, m_resultDirName);
}

void TracerPOTotal::TraceFromFileMultiSize(const std::string &orientFile,
                                            const std::vector<double> &x_sizes,
                                            double x_ref)
{
    // Phase 1: Read orientations
    std::ifstream inFile(orientFile);
    if (!inFile.is_open())
    {
        std::cerr << "Error! Cannot open orientation file: " << orientFile << std::endl;
        throw std::exception();
    }

    std::vector<std::pair<double,double>> orientations;
    std::string line;
    while (std::getline(inFile, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        double b, g;
        if (iss >> b >> g)
            orientations.push_back({b, g});
    }
    inFile.close();

    int nOrientations = orientations.size();
    if (nOrientations == 0)
    {
        std::cerr << "Error! No orientations in file: " << orientFile << std::endl;
        throw std::exception();
    }

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO)
    {
        std::cerr << "Error! Handler is not HandlerPO in TraceFromFileMultiSize" << std::endl;
        throw std::exception();
    }

    // D_ref = maximal dimension of the particle as traced.
    // The particle is already set to x_ref size by the caller.
    double D_ref = m_particle->MaximalDimention();

    CalcTimer timer;
    long long count = 0;
    std::vector<Beam> outBeams;
    timer.Start();
    OutputStartTime(timer);

    double weight = 1.0 / nOrientations;

    // Phase 1: Trace and cache
    BeamCache cache;
    cache.D_ref = D_ref;
    cache.orientations.resize(nOrientations);

    auto t_cache_start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < nOrientations; ++i)
    {
        double beta  = orientations[i].first;
        double gamma = orientations[i].second;

        m_particle->Rotate(beta, gamma, 0);

        if (!shadowOff)
        {
            m_scattering->FormShadowBeam(outBeams);
        }

        bool ok = m_scattering->ScatterLight(0, 0, outBeams);

        double incomingE = m_scattering->GetIncedentEnergy() * weight;

        if (ok)
        {
            handlerPO->CacheBeams(outBeams, weight, D_ref, incomingE,
                                   cache.orientations[i]);
        }
        else
        {
            cache.orientations[i].weight = weight;
            cache.orientations[i].incomingEnergy = incomingE;
            std::cout << std::endl << "Orientation " << i
                      << " (beta=" << beta << ", gamma=" << gamma
                      << ") has been skipped!!!" << std::endl;
        }

        OutputProgress(nOrientations, count, i, 0, timer, outBeams.size());
        outBeams.clear();
        ++count;
    }

    auto t_cache_end = std::chrono::high_resolution_clock::now();
    double cache_sec = std::chrono::duration<double>(t_cache_end - t_cache_start).count();

    EraseConsoleLine(60);
    std::cout << "Phase 1 (tracing): 100%" << std::endl;
    std::cout << "Cached " << cache.totalBeams() << " beams from "
              << nOrientations << " orientations in " << cache_sec << " s" << std::endl;

    // Phase 2: Compute diffraction for all sizes
    auto t_diff_start = std::chrono::high_resolution_clock::now();

    std::vector<Arr2D> results_M;
    std::vector<double> results_energy;
    handlerPO->ComputeFromCache(cache, x_sizes, results_M, results_energy);

    auto t_diff_end = std::chrono::high_resolution_clock::now();
    double diff_sec = std::chrono::duration<double>(t_diff_end - t_diff_start).count();

    std::cout << "Phase 2 (diffraction for " << x_sizes.size() << " sizes): "
              << diff_sec << " s" << std::endl;
    std::cout << "Total: " << cache_sec + diff_sec << " s" << std::endl;

    // Write results for each size
    std::string dir = CreateFolder(m_resultDirName);
#ifdef _WIN32
    std::string baseName = m_resultDirName + '\\' + m_resultDirName;
#else
    std::string baseName = dir + m_resultDirName;
#endif

    // Save original M, write each size's result
    Arr2D &origM = handlerPO->M;

    for (size_t s = 0; s < x_sizes.size(); ++s)
    {
        // Swap in this size's Mueller matrix
        Arr2D savedM = origM;
        origM = results_M[s];
        m_incomingEnergy = results_energy[s];

        // Write with size suffix
        std::string sizeName = baseName + "_x" + std::to_string((int)x_sizes[s]);
        handlerPO->WriteMatricesToFile(sizeName, m_incomingEnergy);

        // Also write in HandlerPOTotal format
        // Write the azimuth-averaged format
        {
            std::ofstream outFile(sizeName + ".dat", std::ios::out);
            outFile << std::setprecision(10);
            outFile << "ScAngle 2pi*dcos "
                    "M11 M12 M13 M14 "
                    "M21 M22 M23 M24 "
                    "M31 M32 M33 M34 "
                    "M41 M42 M43 M44";

            auto &sphere = handlerPO->m_sphere;
            matrix *Lp = handlerPO->m_Lp;
            int nZen = sphere.nZenith;
            int nAz = sphere.nAzimuth;

            for (int iZen = 0; iZen <= nZen; ++iZen)
            {
                matrix Msum(4, 4);
                Msum.Fill(0.0);
                double radZen = sphere.GetZenith(iZen);

                for (int iAz = 0; iAz <= nAz; ++iAz)
                {
                    double radAz = -iAz * sphere.azinuthStep;
                    matrix m = results_M[s](iAz, iZen);

                    (*Lp)[1][1] = cos(2*radAz);
                    (*Lp)[1][2] = sin(2*radAz);
                    (*Lp)[2][1] = -(*Lp)[1][2];
                    (*Lp)[2][2] = (*Lp)[1][1];

                    matrix Ln = *Lp;
                    Ln[1][2] *= -1;
                    Ln[2][1] *= -1;

                    if (radZen > M_PI - __FLT_EPSILON__)
                        Msum += (*Lp) * m * (*Lp);
                    else if (radZen < __FLT_EPSILON__)
                        Msum += Ln * m * (*Lp);
                    else
                        Msum += m * (*Lp);
                }

                double _2Pi_dcos = sphere.Compute2PiDcos(iZen);

                Msum /= nAz;
                outFile << std::endl << RadToDeg(radZen) << ' ' << _2Pi_dcos << ' ';
                outFile << Msum;
            }
            outFile.close();
        }

        // Restore
        origM = savedM;
    }

    std::cout << "Results written to " << baseName << "_x*.dat" << std::endl;
}
