#include "TracerPOTotal.h"
#include "HandlerPOTotal.h"
#include "HandlerPO.h"
#include "BeamCache.h"
#include "Sobol.h"

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
    // Phase 1 (sequential): Trace beams for all orientations, preprocess them.
    //
    // WHY SEQUENTIAL: Particle::Rotate() and ScatterLight() modify shared
    // state (particle geometry, scattering internals). Parallelizing would
    // require N_THREADS copies of Particle and Scattering (no Clone() method
    // exists, and the virtual hierarchy makes it complex to implement).
    //
    // PERFORMANCE: Phase 1 is typically <5% of total time. The dominant
    // Phase 2 (diffraction integrals) IS parallelized with OpenMP.
    // =========================================================================
    // =========================================================================
    // Chunked streaming: process orientations in chunks to limit memory.
    // Each chunk: Phase 1 (sequential trace) → Phase 2 (parallel diffraction).
    // Memory: O(chunkSize * beams_per_orient * sizeof(PreparedBeam)).
    //
    // Chunk size: auto-select based on available RAM.
    // PreparedBeam ~ 3.5 KB, ~100 beams/orient → ~350 KB/orient.
    // Target: use at most 50% of available RAM for beam storage.
    // =========================================================================
    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    // Estimate available memory
    long long availMemBytes = 2LL * 1024 * 1024 * 1024; // default 2 GB
#ifdef __linux__
    {
        std::ifstream meminfo("/proc/meminfo");
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.find("MemAvailable:") == 0) {
                long long kb = 0;
                sscanf(line.c_str(), "MemAvailable: %lld", &kb);
                if (kb > 0) availMemBytes = kb * 1024;
                break;
            }
        }
    }
#endif
    // Reserve memory for thread-local Mueller arrays: nThreads * nAz * nZen * 16 * 8 bytes
    int nThreads = 1;
    #pragma omp parallel
    {
        #pragma omp single
        nThreads = omp_get_num_threads();
    }
    long long muellerMem = (long long)nThreads * (nAz+1) * (nZen+1) * 16 * 8 * 2; // Mueller + Jones
    long long beamBudget = (availMemBytes / 2) - muellerMem; // 50% of RAM for beams
    if (beamBudget < 100LL * 1024 * 1024) beamBudget = 100LL * 1024 * 1024; // min 100 MB

    long long bytesPerOrient = 350LL * 1024; // ~350 KB estimated (100 beams × 3.5 KB)
    int chunkSize = std::max(32, std::min(nOrientations, (int)(beamBudget / bytesPerOrient)));
    // Round up to nice number for Sobol (power of 2 or at least multiple of nThreads)
    if (chunkSize >= nOrientations) chunkSize = nOrientations;

    int nChunks = (nOrientations + chunkSize - 1) / chunkSize;

    std::cerr << "Memory: " << availMemBytes / (1024*1024) << " MB available, "
              << "chunk=" << chunkSize << " orientations (" << nChunks << " chunks), "
              << nThreads << " threads" << std::endl;

    auto t_total_start = std::chrono::high_resolution_clock::now();
    double phase1_total = 0, phase2_total = 0;

    std::vector<Beam> outBeams;
    long long count = 0;

    for (int chunk = 0; chunk < nChunks; ++chunk)
    {
        int iStart = chunk * chunkSize;
        int iEnd = std::min(iStart + chunkSize, nOrientations);
        int thisChunkSize = iEnd - iStart;

        // Phase 1: trace and preprocess this chunk
        auto t_p1 = std::chrono::high_resolution_clock::now();

        std::vector<PreparedOrientation> chunkPrepared(thisChunkSize);

        for (int i = 0; i < thisChunkSize; ++i)
        {
            int idx = iStart + i;
            m_particle->Rotate(orientations[idx].first, orientations[idx].second, 0);

            if (!shadowOff)
                m_scattering->FormShadowBeam(outBeams);

            bool ok = m_scattering->ScatterLight(0, 0, outBeams);

            if (ok)
                handlerPO->PrepareBeams(outBeams, weight, chunkPrepared[i]);
            else
                chunkPrepared[i].sinZenith = weight;

            m_incomingEnergy += m_scattering->GetIncedentEnergy() * weight;
            OutputProgress(nOrientations, count, iStart + i, 0, timer, outBeams.size());
            outBeams.clear();
            ++count;
        }

        auto t_p1_end = std::chrono::high_resolution_clock::now();
        phase1_total += std::chrono::duration<double>(t_p1_end - t_p1).count();

        // Phase 2: parallel diffraction for this chunk
        auto t_p2 = std::chrono::high_resolution_clock::now();

        #pragma omp parallel
        {
            Arr2D localM(nAz + 1, nZen + 1, 4, 4);
            localM.ClearArr();

            std::vector<Arr2DC> localJ;
            if (handlerPO->isCoh)
            {
                Arr2DC tmp(nAz + 1, nZen + 1, 2, 2);
                tmp.ClearArr();
                localJ.push_back(tmp);
            }

            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < thisChunkSize; ++i)
            {
                if (!chunkPrepared[i].beams.empty())
                {
                    handlerPO->HandleBeamsToLocal(chunkPrepared[i], localM, localJ);
                }

                if (handlerPO->isCoh && !localJ.empty())
                {
                    double w = chunkPrepared[i].sinZenith;
                    HandlerPO::AddToMuellerLocal(localJ, w, localM, nAz, nZen);
                    localJ[0].ClearArr();
                }
            }

            #pragma omp critical
            {
                for (int p = 0; p <= nAz; ++p)
                    for (int t = 0; t <= nZen; ++t)
                    {
                        matrix m = localM(p, t);
                        handlerPO->M.insert(p, t, m);
                    }
            }
        } // end omp parallel

        auto t_p2_end = std::chrono::high_resolution_clock::now();
        phase2_total += std::chrono::duration<double>(t_p2_end - t_p2).count();

        // Free chunk memory immediately
        chunkPrepared.clear();
        chunkPrepared.shrink_to_fit();
    }

    EraseConsoleLine(60);
    std::cout << "Phase 1 (tracing + preprocessing): " << std::fixed
              << std::setprecision(2) << phase1_total << " s" << std::endl;

    // (t_phase2_end - t_phase2_start) compatibility: set from totals
    auto t_phase2_start = t_total_start; // dummy for downstream code
    auto t_phase2_end = std::chrono::high_resolution_clock::now();
    double phase2_sec = std::chrono::duration<double>(t_phase2_end - t_phase2_start).count();

    std::cout << "Phase 2 (diffraction, OpenMP): " << std::fixed
              << std::setprecision(2) << phase2_total << " s" << std::endl;
    std::cout << "Total: " << phase1_total + phase2_total << " s" << std::endl;

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

void TracerPOTotal::TraceFromSobol(int nOrient, double betaSym, double gammaSym)
{
    // Generate Sobol orientations mapped to [0, betaSym] x [0, gammaSym]
    // using the correct solid-angle measure:
    //   beta = arccos(1 - (1-cos(betaSym)) * u)  for uniform dOmega
    //   gamma = gammaSym * v
    Sobol2D sobol(42);
    std::vector<double> su, sv;
    sobol.generate(nOrient, su, sv);

    double cosBetaSym = cos(betaSym);
    std::vector<std::pair<double,double>> orientations(nOrient);
    for (int i = 0; i < nOrient; ++i)
    {
        double beta  = acos(1.0 - (1.0 - cosBetaSym) * su[i]);
        double gamma = gammaSym * sv[i];
        orientations[i] = {beta, gamma};
    }

    std::cout << "Sobol: " << nOrient << " orientations, beta_sym="
              << RadToDeg(betaSym) << " deg, gamma_sym="
              << RadToDeg(gammaSym) << " deg" << std::endl;

    CalcTimer timer;
    timer.Start();
    OutputStartTime(timer);

    double weight = 1.0 / nOrient;
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
        std::cerr << "Error! Handler is not HandlerPO in TraceFromSobol" << std::endl;
        throw std::exception();
    }

    // =========================================================================
    // Phase 1 (sequential): Trace beams for all orientations.
    //
    // WHY SEQUENTIAL: Particle::Rotate() modifies shared particle geometry
    // (vertices, normals, facets) and ScatterLight() modifies shared
    // Scattering state. Making these thread-safe would require either:
    //   (a) Particle::Clone() — complex due to virtual hierarchy, or
    //   (b) Critical sections around stateful parts — serializes the work.
    //
    // PERFORMANCE IMPACT: Phase 1 typically takes <5% of total time
    // (e.g., ~0.5s out of ~30s for x=50, 512 orientations). The dominant
    // cost is Phase 2 (diffraction integrals), which IS parallelized.
    // Parallelizing Phase 1 would save <2% of total runtime — not worth
    // the complexity (Amdahl's law).
    // =========================================================================
    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    // Chunked streaming: auto-size chunks based on available RAM
    long long availMB = 2048; // default 2 GB
#ifdef __linux__
    {
        std::ifstream meminfo("/proc/meminfo");
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.find("MemAvailable:") == 0) {
                long long kb = 0;
                sscanf(line.c_str(), "MemAvailable: %lld", &kb);
                if (kb > 0) availMB = kb / 1024;
                break;
            }
        }
    }
#endif
    long long beamBudget = std::max(100LL, availMB / 2); // 50% of RAM, min 100 MB
    int chunkSize = std::max(32, std::min(nOrient, (int)(beamBudget * 1024 / 350)));
    int nChunks = (nOrient + chunkSize - 1) / chunkSize;

    std::cerr << "Memory: " << availMB << " MB available, chunk="
              << chunkSize << " orientations (" << nChunks << " chunks)" << std::endl;

    double phase1_total = 0, phase2_total = 0;
    std::vector<Beam> outBeams;
    long long count = 0;

    for (int chunk = 0; chunk < nChunks; ++chunk)
    {
        int iStart = chunk * chunkSize;
        int iEnd = std::min(iStart + chunkSize, nOrient);
        int thisChunk = iEnd - iStart;

        // Phase 1: trace this chunk
        auto tp1 = std::chrono::high_resolution_clock::now();
        std::vector<PreparedOrientation> chunkPrepared(thisChunk);

        for (int i = 0; i < thisChunk; ++i)
        {
            int idx = iStart + i;
            m_particle->Rotate(orientations[idx].first, orientations[idx].second, 0);
            if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
            bool ok = m_scattering->ScatterLight(0, 0, outBeams);
            if (ok) handlerPO->PrepareBeams(outBeams, weight, chunkPrepared[i]);
            else    chunkPrepared[i].sinZenith = weight;
            m_incomingEnergy += m_scattering->GetIncedentEnergy() * weight;
            OutputProgress(nOrient, count, iStart + i, 0, timer, outBeams.size());
            outBeams.clear();
            ++count;
        }
        phase1_total += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tp1).count();

        // Phase 2: parallel diffraction for this chunk
        auto tp2 = std::chrono::high_resolution_clock::now();

        #pragma omp parallel
        {
            Arr2D localM(nAz + 1, nZen + 1, 4, 4);
            localM.ClearArr();
            std::vector<Arr2DC> localJ;
            if (handlerPO->isCoh) {
                Arr2DC tmp(nAz + 1, nZen + 1, 2, 2);
                tmp.ClearArr();
                localJ.push_back(tmp);
            }

            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < thisChunk; ++i)
            {
                if (!chunkPrepared[i].beams.empty())
                    handlerPO->HandleBeamsToLocal(chunkPrepared[i], localM, localJ);
                if (handlerPO->isCoh && !localJ.empty()) {
                    double w = chunkPrepared[i].sinZenith;
                    HandlerPO::AddToMuellerLocal(localJ, w, localM, nAz, nZen);
                    localJ[0].ClearArr();
                }
            }

            #pragma omp critical
            {
                for (int p = 0; p <= nAz; ++p)
                    for (int t = 0; t <= nZen; ++t) {
                        matrix m = localM(p, t);
                        handlerPO->M.insert(p, t, m);
                    }
            }
        }
        phase2_total += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tp2).count();

        chunkPrepared.clear();
        chunkPrepared.shrink_to_fit();
    }

    EraseConsoleLine(60);
    std::cout << "Phase 1 (tracing + preprocessing): " << std::fixed
              << std::setprecision(2) << phase1_total << " s" << std::endl;
    std::cout << "Phase 2 (diffraction, OpenMP): " << std::fixed
              << std::setprecision(2) << phase2_total << " s" << std::endl;
    std::cout << "Total: " << phase1_total + phase2_total << " s" << std::endl;

    m_handler->WriteTotalMatricesToFile(m_resultDirName);
    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
    OutputStatisticsPO(timer, nOrient, m_resultDirName);
}

void TracerPOTotal::TraceAdaptive(double eps, double betaSym, double gammaSym)
{
    std::cout << "Adaptive mode: target relative change = " << eps << std::endl;
    std::cout << "  Convergence criterion: M11(180°) backscattering" << std::endl;
    std::cout << "  beta_sym=" << RadToDeg(betaSym) << " deg, gamma_sym="
              << RadToDeg(gammaSym) << " deg" << std::endl;

    HandlerPO *hp = dynamic_cast<HandlerPO*>(m_handler);
    if (!hp)
    {
        std::cerr << "Error: handler is not HandlerPO" << std::endl;
        TraceFromSobol(1024, betaSym, gammaSym);
        return;
    }

    double cosBetaSym = cos(betaSym);
    int nZen = hp->m_sphere.nZenith;
    int nAz = hp->m_sphere.nAzimuth;

    // Restart-from-scratch adaptive: each iteration runs full TraceFromSobol.
    // Sobol property ensures first 2^m points are optimal → no wasted work
    // since each larger run includes all smaller runs' physics (same orientations).
    // Total overhead: sum of geometric series = 2×final cost (worst case 2× vs optimal).

    // Start small, scale up. Sobol property: first 2^m points are always optimal.
    // For small x (<50): C_sca converges fast, start at 64.
    // For large x (>500): backscattering needs more, start at 256.
    int nOrient = 64;
    // Max orientations: limited by available RAM (~350 KB per orientation).
    // Default 16384 (~5.5 GB). Auto-scale with available memory.
    long long availMB_ad = 2048;
#ifdef __linux__
    {
        std::ifstream meminfo("/proc/meminfo");
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.find("MemAvailable:") == 0) {
                long long kb = 0;
                sscanf(line.c_str(), "MemAvailable: %lld", &kb);
                if (kb > 0) availMB_ad = kb / 1024;
                break;
            }
        }
    }
#endif
    // Use 50% of RAM, ~350 KB/orient. Min 1024, cap at 131072.
    int maxOrient = std::max(1024, std::min(131072, (int)(availMB_ad / 2 * 1024 / 350)));
    // Round down to power of 2 (Sobol property)
    { int p = 1; while (p * 2 <= maxOrient) p *= 2; maxOrient = p; }
    std::cerr << "Adaptive: max orientations = " << maxOrient
              << " (" << availMB_ad << " MB available)" << std::endl;
    double prevM11_180 = 0;
    double prevCsca = 0;

    for (int iter = 0; iter < 10; ++iter)
    {
        // Full clean + compute (sequential, no OpenMP — avoids parallel bugs)
        hp->M.ClearArr();
        hp->CleanJ();
        m_incomingEnergy = 0;
        hp->m_outputEnergy = 0;

        Sobol2D sobol_gen(42);
        std::vector<double> su, sv;
        sobol_gen.generate(nOrient, su, sv);

        double weight = 1.0 / nOrient;
        std::vector<Beam> outBeams;

        for (int i = 0; i < nOrient; ++i)
        {
            double beta = acos(1.0 - (1.0 - cosBetaSym) * su[i]);
            double gamma = gammaSym * sv[i];
            m_particle->Rotate(beta, gamma, 0);
            if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
            bool ok = m_scattering->ScatterLight(0, 0, outBeams);
            if (ok) m_handler->HandleBeams(outBeams, weight);
            m_incomingEnergy += m_scattering->GetIncedentEnergy() * weight;
            outBeams.clear();
        }

        // Extract M11(180°) from raw Mueller (phi-averaged)
        double M11_180_avg = 0;
        for (int p = 0; p <= nAz; ++p)
            M11_180_avg += hp->M(p, nZen)[0][0];
        M11_180_avg /= (nAz + 1);

        double Csca = m_incomingEnergy;

        double relChange_m11 = (prevM11_180 != 0) ?
            fabs(M11_180_avg - prevM11_180) / fabs(prevM11_180) : 1.0;
        double relChange_csca = (prevCsca > 0) ?
            fabs(Csca - prevCsca) / prevCsca : 1.0;

        std::cout << "  N=" << nOrient
                  << " (+" << nOrient
                  << " new), M11(180)=" << M11_180_avg
                  << ", dM11=" << relChange_m11*100 << "%"
                  << ", Csca=" << Csca
                  << ", dCsca=" << relChange_csca*100 << "%"
                  << std::endl;

        // Convergence: C_sca must converge. M11(180°) is bonus.
        bool csca_ok = (relChange_csca < eps && iter > 0);
        bool m11_ok = (relChange_m11 < eps && iter > 0);
        // Relaxed M11 criterion: backscattering is noisy, accept 5x eps
        bool m11_relaxed = (relChange_m11 < 5.0 * eps && iter > 1);

        if (csca_ok && m11_ok)
        {
            std::cout << "Converged at N=" << nOrient << " (both C_sca and M11)" << std::endl;
            return;
        }
        if (csca_ok && m11_relaxed)
        {
            std::cout << "Converged at N=" << nOrient
                      << " (C_sca ok, M11 within " << relChange_m11*100 << "%)" << std::endl;
            return;
        }
        // Hard stop: C_sca very stable (eps/10) — M11 will converge eventually
        if (relChange_csca < eps * 0.1 && iter >= 4)
        {
            std::cout << "Converged at N=" << nOrient
                      << " (C_sca stable, M11 still " << relChange_m11*100 << "%)" << std::endl;
            return;
        }

        prevM11_180 = M11_180_avg;
        prevCsca = Csca;
        nOrient *= 2;
        if (nOrient > maxOrient) break;
    }

    std::cout << "Max iterations reached, using N=" << nOrient/2 << std::endl;
}
