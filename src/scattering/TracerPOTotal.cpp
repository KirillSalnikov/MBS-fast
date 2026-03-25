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
#include <cstdio>
#include <sys/stat.h>

// ---- Checkpoint save/load for resume after crash ----
static void SaveCheckpoint(const std::string &path, const Arr2D &M, const Arr2D &M_ns,
                            double energy, int completedOrient, int totalOrient,
                            unsigned long paramHash, int nAz, int nZen)
{
    std::ofstream f(path, std::ios::binary);
    if (!f.is_open()) return;
    uint32_t magic = 0x4D425343; // "MBSC"
    f.write((char*)&magic, 4);
    f.write((char*)&paramHash, sizeof(paramHash));
    f.write((char*)&completedOrient, sizeof(completedOrient));
    f.write((char*)&totalOrient, sizeof(totalOrient));
    f.write((char*)&energy, sizeof(energy));
    f.write((char*)&nAz, 4);
    f.write((char*)&nZen, 4);
    // Dump M data
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t < nZen; ++t) {
            matrix m = M(p, t);
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c) {
                    double v = m[r][c];
                    f.write((char*)&v, 8);
                }
        }
    // Dump M_noshadow
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t < nZen; ++t) {
            matrix m = M_ns(p, t);
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c) {
                    double v = m[r][c];
                    f.write((char*)&v, 8);
                }
        }
    f.close();
}

static bool LoadCheckpoint(const std::string &path, Arr2D &M, Arr2D &M_ns,
                            double &energy, int &completedOrient, int totalOrient,
                            unsigned long paramHash)
{
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) return false;
    uint32_t magic; f.read((char*)&magic, 4);
    if (magic != 0x4D425343) return false;
    unsigned long storedHash; f.read((char*)&storedHash, sizeof(storedHash));
    if (storedHash != paramHash) {
        std::cerr << "Checkpoint param mismatch, ignoring" << std::endl;
        return false;
    }
    int storedCompleted, storedTotal;
    f.read((char*)&storedCompleted, sizeof(storedCompleted));
    f.read((char*)&storedTotal, sizeof(storedTotal));
    if (storedTotal != totalOrient) return false;
    f.read((char*)&energy, sizeof(energy));
    int nAz, nZen; f.read((char*)&nAz, 4); f.read((char*)&nZen, 4);
    // Load M
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t < nZen; ++t)
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c) {
                    double v; f.read((char*)&v, 8);
                    M(p, t)[r][c] = v;
                }
    // Load M_noshadow
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t < nZen; ++t)
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c) {
                    double v; f.read((char*)&v, 8);
                    M_ns(p, t)[r][c] = v;
                }
    completedOrient = storedCompleted;
    f.close();
    return true;
}

static unsigned long HashParams(double wave, double ri_re, int nActs, int nOrient,
                                 double L, double D, int nAz, int nZen)
{
    unsigned long h = 0;
    auto mix = [&](double v) { h ^= std::hash<double>{}(v) + 0x9e3779b9 + (h<<6) + (h>>2); };
    auto mixi = [&](int v) { h ^= std::hash<int>{}(v) + 0x9e3779b9 + (h<<6) + (h>>2); };
    mix(wave); mix(ri_re); mixi(nActs); mixi(nOrient); mix(L); mix(D); mixi(nAz); mixi(nZen);
    return h;
}
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace std;

// Helper: pack Arr2D (N x M grid of n x m matrices) into flat double array
static void Arr2DToFlat(const Arr2D &arr, int N, int M, double *buf)
{
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j) {
            matrix m = arr(i, j);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    buf[((i*M + j)*4 + a)*4 + b] = m[a][b];
        }
}

// Helper: unpack flat double array into Arr2D (replace mode)
static void FlatToArr2D(const double *buf, int N, int M, Arr2D &arr)
{
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j) {
            matrix mt(4, 4);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    mt[a][b] = buf[((i*M + j)*4 + a)*4 + b];
            arr.replace(i, j, mt);
        }
}

// Helper: MPI reduce Arr2D + incomingEnergy, rank 0 gets result
static void MPI_ReduceMueller(HandlerPO *hp, int nAz, int nZen,
                               double &incomingEnergy, int mpi_rank)
{
#ifdef USE_MPI
    int totalDoubles = (nAz+1) * (nZen+1) * 16;
    std::vector<double> sendbuf(totalDoubles), recvbuf(totalDoubles);

    // Reduce M
    Arr2DToFlat(hp->M, nAz+1, nZen+1, sendbuf.data());
    MPI_Reduce(sendbuf.data(), recvbuf.data(), totalDoubles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) FlatToArr2D(recvbuf.data(), nAz+1, nZen+1, hp->M);

    // Reduce M_noshadow
    Arr2DToFlat(hp->M_noshadow, nAz+1, nZen+1, sendbuf.data());
    MPI_Reduce(sendbuf.data(), recvbuf.data(), totalDoubles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) FlatToArr2D(recvbuf.data(), nAz+1, nZen+1, hp->M_noshadow);

    // Reduce incomingEnergy
    double totalEnergy = 0;
    MPI_Reduce(&incomingEnergy, &totalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) incomingEnergy = totalEnergy;
#endif
}

TracerPOTotal::TracerPOTotal(Particle *particle, int nActs,
                             const string &resultFileName)
    : TracerPO(particle, nActs, resultFileName)
{
}

void TracerPOTotal::TraceRandom(const AngleRange &betaRange,
                                const AngleRange &gammaRange)
{
    int betaNorm = (m_symmetry.beta < M_PI_2+FLT_EPSILON && m_symmetry.beta > M_PI_2-FLT_EPSILON) ? 1 : 2;
    double normGamma = gammaRange.number * betaNorm;

    // --coh_orient: legacy mode — HandleBeams accumulates Jones coherently
    // across ALL orientations, then single AddToMueller at end.
    // No chunking, no OpenMP (matches old MBS-raw exactly).
    if (m_cohOrient) {
        m_handler->SetNormIndex(normGamma);
        std::string dir = CreateFolder(m_resultDirName);
        m_resultDirName = dir + m_resultDirName;
        vector<Beam> outBeams;
        for (int i = 0; i <= betaRange.number; ++i) {
            double beta = betaRange.min + i * betaRange.step;
            double dcos;
            CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normGamma, dcos);
            m_handler->SetSinZenith(dcos);
            for (int j = 0; j < gammaRange.number; ++j) {
                double gamma = gammaRange.min + j * gammaRange.step;
                m_particle->Rotate(beta, gamma, 0);
                if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
                m_scattering->ScatterLight(0, 0, outBeams);
                m_handler->HandleBeams(outBeams, dcos);
                m_incomingEnergy += m_scattering->GetIncedentEnergy() * dcos;
                outBeams.clear();
            }
        }
        m_handler->WriteTotalMatricesToFile(m_resultDirName);
        m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
        CalcTimer timer; timer.Start();
        OutputStatisticsPO(timer, (betaRange.number+1)*gammaRange.number, m_resultDirName);
        return;
    }

    // Default: incoherent per-orientation (physically correct)
    // Loop over beta explicitly — enables per-beta Mueller saving

    int nGamma = gammaRange.number;
    int nBeta = betaRange.number + 1;
    int nOrientations = nBeta * nGamma;

    if (m_mpiRank == 0)
        std::cout << "Random grid: " << nBeta << " x " << nGamma
                  << " = " << nOrientations << " orientations" << std::endl;

    CalcTimer timer;
    timer.Start();
    if (m_mpiRank == 0) OutputStartTime(timer);

    m_handler->SetNormIndex(normGamma);

    std::string dir;
    if (m_mpiRank == 0) dir = CreateFolder(m_resultDirName);
#ifdef USE_MPI
    if (m_mpiSize > 1) {
        int dirLen = dir.size();
        MPI_Bcast(&dirLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
        dir.resize(dirLen);
        MPI_Bcast(&dir[0], dirLen, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
#ifdef _WIN32
    m_resultDirName += '\\' + m_resultDirName;
#else
    m_resultDirName = dir + m_resultDirName;
#endif

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) {
        std::cerr << "Error: handler is not HandlerPO in TraceRandom" << std::endl;
        return;
    }

    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    // --save_betas: create output directory
    std::string betaDir;
    if (m_saveBetas && m_mpiRank == 0) {
        betaDir = m_resultDirName + "_betas";
        mkdir(betaDir.c_str(), 0755);
    }

    double phase1_total = 0, phase2_total = 0;
    std::vector<Beam> outBeams;
    long long count = 0;

    // Process beta-by-beta: each beta = one chunk of nGamma orientations
    for (int ib = 0; ib < nBeta; ++ib)
    {
        double beta = betaRange.min + ib * betaRange.step;
        double dcos;
        CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normGamma, dcos);

        // Phase 1: trace all gamma for this beta
        auto tp1 = std::chrono::high_resolution_clock::now();
        std::vector<PreparedOrientation> chunkPrepared(nGamma);

        for (int j = 0; j < nGamma; ++j)
        {
            double gamma = gammaRange.min + j * gammaRange.step;
            m_particle->Rotate(beta, gamma, 0);
            if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
            bool ok = m_scattering->ScatterLight(0, 0, outBeams);
            if (ok) handlerPO->PrepareBeams(outBeams, dcos, chunkPrepared[j]);
            else    chunkPrepared[j].sinZenith = dcos;
            m_incomingEnergy += m_scattering->GetIncedentEnergy() * dcos;
            if (m_mpiRank == 0) OutputProgress(nOrientations, count, ib*nGamma+j, 0, timer, outBeams.size());
            outBeams.clear();
            ++count;
        }
        phase1_total += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tp1).count();

        // Phase 2: diffraction (OpenMP parallel over gamma orientations)
        auto tp2 = std::chrono::high_resolution_clock::now();
        Arr2D betaM(nAz+1, nZen+1, 4, 4); betaM.ClearArr();
        Arr2D betaM_ns(nAz+1, nZen+1, 4, 4); betaM_ns.ClearArr();

        #pragma omp parallel
        {
            Arr2D localM(nAz+1, nZen+1, 4, 4); localM.ClearArr();
            Arr2D localM_ns(nAz+1, nZen+1, 4, 4); localM_ns.ClearArr();
            std::vector<Arr2DC> localJ, localJ_ns;
            if (handlerPO->isCoh) {
                Arr2DC t1(nAz+1,nZen+1,2,2); t1.ClearArr(); localJ.push_back(t1);
                Arr2DC t2(nAz+1,nZen+1,2,2); t2.ClearArr(); localJ_ns.push_back(t2);
            }
            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < nGamma; ++i) {
                if (!chunkPrepared[i].beams.empty())
                    handlerPO->HandleBeamsToLocal(chunkPrepared[i], localM, localJ,
                                                   handlerPO->isCoh ? &localJ_ns : nullptr);
                if (handlerPO->isCoh && !localJ.empty()) {
                    double w = chunkPrepared[i].sinZenith;
                    HandlerPO::AddToMuellerLocal(localJ, w, localM, nAz, nZen);
                    HandlerPO::AddToMuellerLocal(localJ_ns, w, localM_ns, nAz, nZen);
                    localJ[0].ClearArr(); localJ_ns[0].ClearArr();
                }
            }
            #pragma omp critical
            { for (int p=0;p<=nAz;++p) for (int t=0;t<=nZen;++t) {
                betaM.insert(p,t,localM(p,t));
                betaM_ns.insert(p,t,localM_ns(p,t));
            } }
        }
        phase2_total += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tp2).count();
        chunkPrepared.clear(); chunkPrepared.shrink_to_fit();

        // Accumulate into global Mueller
        for (int p=0;p<=nAz;++p) for (int t=0;t<=nZen;++t) {
            handlerPO->M.insert(p,t,betaM(p,t));
            handlerPO->M_noshadow.insert(p,t,betaM_ns(p,t));
        }

        // --save_betas: write per-beta Mueller (phi-averaged)
        if (m_saveBetas && m_mpiRank == 0)
        {
            auto &sphere = handlerPO->m_sphere;
            matrix *Lp = handlerPO->m_Lp;

            // Per-beta contribution
            std::string fname = betaDir + "/beta_" + std::to_string(ib)
                + "_" + std::to_string((int)RadToDeg(beta)) + "deg.dat";
            std::ofstream bf(fname, std::ios::out);
            bf << std::setprecision(10);
            bf << "# Per-beta Mueller: beta=" << RadToDeg(beta) << " deg, dcos=" << dcos
               << ", " << nGamma << " gamma orientations" << std::endl;
            bf << "ScAngle 2pi*dcos M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44";
            for (int iZen = 0; iZen <= nZen; ++iZen) {
                matrix Msum(4,4); Msum.Fill(0.0);
                double radZen = sphere.GetZenith(iZen);
                for (int iAz = 0; iAz <= nAz; ++iAz) {
                    double radAz = -iAz * sphere.azinuthStep;
                    matrix m = betaM(iAz, iZen);
                    (*Lp)[1][1] = cos(2*radAz); (*Lp)[1][2] = sin(2*radAz);
                    (*Lp)[2][1] = -(*Lp)[1][2]; (*Lp)[2][2] = (*Lp)[1][1];
                    Msum += m * (*Lp);
                }
                Msum /= nAz;
                double _2PiDcos = sphere.Compute2PiDcos(iZen);
                bf << std::endl << RadToDeg(radZen) << ' ' << _2PiDcos << ' ';
                bf << Msum;
            }
            bf.close();

            // Also write cumulative
            std::string cfname = betaDir + "/cumul_" + std::to_string(ib)
                + "_" + std::to_string((int)RadToDeg(beta)) + "deg.dat";
            std::ofstream cf(cfname, std::ios::out);
            cf << std::setprecision(10);
            cf << "# Cumulative Mueller up to beta=" << RadToDeg(beta) << " deg"
               << " (" << (ib+1)*nGamma << "/" << nOrientations << " orientations)" << std::endl;
            cf << "ScAngle 2pi*dcos M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44";
            for (int iZen = 0; iZen <= nZen; ++iZen) {
                matrix Msum(4,4); Msum.Fill(0.0);
                double radZen = sphere.GetZenith(iZen);
                for (int iAz = 0; iAz <= nAz; ++iAz) {
                    double radAz = -iAz * sphere.azinuthStep;
                    matrix m = handlerPO->M(iAz, iZen);
                    (*Lp)[1][1] = cos(2*radAz); (*Lp)[1][2] = sin(2*radAz);
                    (*Lp)[2][1] = -(*Lp)[1][2]; (*Lp)[2][2] = (*Lp)[1][1];
                    Msum += m * (*Lp);
                }
                Msum /= nAz;
                double _2PiDcos = sphere.Compute2PiDcos(iZen);
                cf << std::endl << RadToDeg(radZen) << ' ' << _2PiDcos << ' ';
                cf << Msum;
            }
            cf.close();
        }
    }

    MPI_ReduceMueller(handlerPO, nAz, nZen, m_incomingEnergy, m_mpiRank);

    EraseConsoleLine(60);
    if (m_mpiRank == 0) {
        std::cout << "Phase 1 (tracing): " << std::fixed << std::setprecision(2) << phase1_total << " s" << std::endl;
        std::cout << "Phase 2 (diffraction, OpenMP): " << phase2_total << " s" << std::endl;
        std::cout << "Total: " << phase1_total + phase2_total << " s" << std::endl;
        if (m_saveBetas)
            std::cout << "Saved " << nBeta << " beta files to " << betaDir << "/" << std::endl;

        m_handler->WriteTotalMatricesToFile(m_resultDirName);
        m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
        std::swap(handlerPO->M, handlerPO->M_noshadow);
        std::string nsName = m_resultDirName + "_noshadow";
        handlerPO->WriteMatricesToFile(nsName, m_incomingEnergy);
        std::swap(handlerPO->M, handlerPO->M_noshadow);
        OutputStatisticsPO(timer, nOrientations, m_resultDirName);
    }
}


void TracerPOTotal::TraceMonteCarlo(const AngleRange &betaRange,
                                    const AngleRange &gammaRange,
                                    int nOrientations)
{
    // Generate random orientations, then use same chunked+OpenMP pipeline
    std::vector<std::pair<double,double>> orientations;
    std::vector<double> weights;

    unsigned int lo, hi;
    asm volatile("rdtsc" : "=a"(lo), "=d"(hi));
    srand(lo ^ hi);

    for (int i = 0; i < nOrientations; ++i)
    {
        double beta = RandomDouble(0, 1) * betaRange.max;
        double gamma = RandomDouble(0, 1) * gammaRange.max;
        orientations.push_back({beta, gamma});
        weights.push_back(sin(beta));
    }

    std::cout << "Monte Carlo: " << nOrientations << " orientations" << std::endl;

    CalcTimer timer;
    timer.Start();
    if (m_mpiRank == 0) OutputStartTime(timer);

    m_handler->SetNormIndex(1);

    std::string dir = CreateFolder(m_resultDirName);
#ifdef _WIN32
    m_resultDirName += '\\' + m_resultDirName;
#else
    m_resultDirName = dir + m_resultDirName;
#endif

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) { std::cerr << "Error: not HandlerPO" << std::endl; return; }

    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    // MPI split
    int myStart = m_mpiRank * nOrientations / m_mpiSize;
    int myEnd = (m_mpiRank + 1) * nOrientations / m_mpiSize;
    int myCount = myEnd - myStart;

    // Chunked + OpenMP (same as TraceRandom)
    std::vector<Beam> outBeams;
    for (int i = myStart; i < myEnd; ++i)
    {
        m_particle->Rotate(orientations[i].first, orientations[i].second, 0);
        if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
        bool ok = m_scattering->ScatterLight(0, 0, outBeams);

        if (ok)
        {
            PreparedOrientation prepared;
            handlerPO->PrepareBeams(outBeams, weights[i], prepared);

            std::vector<Arr2DC> localJ, localJ_ns;
            if (handlerPO->isCoh) {
                Arr2DC t1(nAz+1,nZen+1,2,2); t1.ClearArr(); localJ.push_back(t1);
                Arr2DC t2(nAz+1,nZen+1,2,2); t2.ClearArr(); localJ_ns.push_back(t2);
            }
            handlerPO->HandleBeamsToLocal(prepared, handlerPO->M, localJ,
                                           handlerPO->isCoh ? &localJ_ns : nullptr);
            if (handlerPO->isCoh && !localJ.empty()) {
                HandlerPO::AddToMuellerLocal(localJ, prepared.sinZenith,
                                              handlerPO->M, nAz, nZen);
                HandlerPO::AddToMuellerLocal(localJ_ns, prepared.sinZenith,
                                              handlerPO->M_noshadow, nAz, nZen);
            }
        }
        m_incomingEnergy += m_scattering->GetIncedentEnergy() * weights[i];
        outBeams.clear();
    }

    MPI_ReduceMueller(handlerPO, nAz, nZen, m_incomingEnergy, m_mpiRank);

    if (m_mpiRank == 0) {
        m_handler->WriteTotalMatricesToFile(m_resultDirName);
        m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
        std::swap(handlerPO->M, handlerPO->M_noshadow);
        std::string nsName = m_resultDirName + "_noshadow";
        handlerPO->WriteMatricesToFile(nsName, m_incomingEnergy);
        std::swap(handlerPO->M, handlerPO->M_noshadow);
        OutputStatisticsPO(timer, nOrientations, m_resultDirName);
    }
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
    int chunkSize = std::max(32, std::min(4096, std::min(nOrientations, (int)(beamBudget / bytesPerOrient))));
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

    // Checkpoint: resume from crash
    std::string ckptPath = m_resultDirName + "_checkpoint.bin";
    unsigned long paramHash = HashParams(m_scattering->m_wave,
        real(m_particle->GetRefractiveIndex()), 0,
        nOrientations, m_particle->MaximalDimention(), 0, nAz, nZen);
    int resumeChunk = 0;
    {
        int completedOrient = 0;
        if (LoadCheckpoint(ckptPath, handlerPO->M, handlerPO->M_noshadow,
                           m_incomingEnergy, completedOrient, nOrientations, paramHash))
        {
            resumeChunk = completedOrient / chunkSize;
            count = completedOrient;
            if (m_mpiRank == 0)
                std::cout << "*** RESUMED from checkpoint: " << completedOrient
                          << "/" << nOrientations << " orientations (chunk " << resumeChunk
                          << "/" << nChunks << ")" << std::endl;
        }
    }

    for (int chunk = resumeChunk; chunk < nChunks; ++chunk)
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
            if (m_mpiRank == 0) OutputProgress(nOrientations, count, iStart + i, 0, timer, outBeams.size());
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
            Arr2D localM_ns(nAz + 1, nZen + 1, 4, 4); // no-shadow
            localM_ns.ClearArr();

            std::vector<Arr2DC> localJ, localJ_ns;
            if (handlerPO->isCoh)
            {
                Arr2DC tmp(nAz + 1, nZen + 1, 2, 2);
                tmp.ClearArr();
                localJ.push_back(tmp);
                Arr2DC tmp2(nAz + 1, nZen + 1, 2, 2);
                tmp2.ClearArr();
                localJ_ns.push_back(tmp2);
            }

            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < thisChunkSize; ++i)
            {
                if (!chunkPrepared[i].beams.empty())
                {
                    handlerPO->HandleBeamsToLocal(chunkPrepared[i], localM, localJ,
                                                   handlerPO->isCoh ? &localJ_ns : nullptr);
                }

                if (handlerPO->isCoh && !localJ.empty())
                {
                    double w = chunkPrepared[i].sinZenith;
                    HandlerPO::AddToMuellerLocal(localJ, w, localM, nAz, nZen);
                    HandlerPO::AddToMuellerLocal(localJ_ns, w, localM_ns, nAz, nZen);
                    localJ[0].ClearArr();
                    localJ_ns[0].ClearArr();
                }
            }

            #pragma omp critical
            {
                for (int p = 0; p <= nAz; ++p)
                    for (int t = 0; t <= nZen; ++t)
                    {
                        handlerPO->M.insert(p, t, localM(p, t));
                        handlerPO->M_noshadow.insert(p, t, localM_ns(p, t));
                    }
            }
        } // end omp parallel

        auto t_p2_end = std::chrono::high_resolution_clock::now();
        phase2_total += std::chrono::duration<double>(t_p2_end - t_p2).count();

        // Free chunk memory immediately
        chunkPrepared.clear();
        chunkPrepared.shrink_to_fit();

        // Save checkpoint after each chunk
        if (m_mpiRank == 0) {
            SaveCheckpoint(ckptPath, handlerPO->M, handlerPO->M_noshadow,
                           m_incomingEnergy, iEnd, nOrientations, paramHash, nAz+1, nZen+1);
        }
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

    // Write no-shadow Mueller matrix (swap M <-> M_noshadow temporarily)
    {
        HandlerPO *hp = dynamic_cast<HandlerPO*>(m_handler);
        if (hp) {
            std::swap(hp->M, hp->M_noshadow);
            std::string nsName = m_resultDirName + "_noshadow";
            hp->WriteMatricesToFile(nsName, m_incomingEnergy);
            std::swap(hp->M, hp->M_noshadow); // swap back
        }
    }

    OutputStatisticsPO(timer, nOrientations, m_resultDirName);

    // Remove checkpoint on successful completion
    std::remove(ckptPath.c_str());
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
            if (m_mpiRank == 0) std::cout << std::endl << "Orientation " << i
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
    if (m_mpiRank == 0) std::cout << "Phase 1 (tracing): 100%" << std::endl;
    if (m_mpiRank == 0) std::cout << "Cached " << cache.totalBeams() << " beams from "
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

void TracerPOTotal::TraceSobolMultiSize(int nOrient, double betaSym, double gammaSym,
                                         const std::vector<double> &x_sizes, double x_ref)
{
    // Generate Sobol orientations (same as TraceFromSobol)
    Sobol2D sobol(42);
    std::vector<double> su, sv;
    sobol.generate(nOrient, su, sv);

    double cosBetaSym = cos(betaSym);
    std::vector<std::pair<double,double>> orientations(nOrient);
    for (int i = 0; i < nOrient; ++i) {
        double beta  = acos(1.0 - (1.0 - cosBetaSym) * su[i]);
        double gamma = gammaSym * sv[i];
        orientations[i] = {beta, gamma};
    }

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) {
        std::cerr << "Error! Handler is not HandlerPO in TraceSobolMultiSize" << std::endl;
        throw std::exception();
    }

    double D_ref = m_particle->MaximalDimention();
    double weight = 1.0 / nOrient;

    std::string dir = CreateFolder(m_resultDirName);
#ifdef _WIN32
    m_resultDirName += '\\' + m_resultDirName;
#else
    m_resultDirName = dir + m_resultDirName;
#endif

    if (m_mpiRank == 0)
        std::cout << "SobolMultiSize: " << nOrient << " orientations, "
                  << x_sizes.size() << " sizes, x_ref=" << x_ref << std::endl;

    CalcTimer timer;
    timer.Start();
    OutputStartTime(timer);

    // Phase 1: Trace and cache beams
    BeamCache cache;
    cache.D_ref = D_ref;
    cache.orientations.resize(nOrient);

    auto t_cache_start = std::chrono::high_resolution_clock::now();

    std::vector<Beam> outBeams;
    long long count = 0;
    for (int i = 0; i < nOrient; ++i) {
        m_particle->Rotate(orientations[i].first, orientations[i].second, 0);
        if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
        bool ok = m_scattering->ScatterLight(0, 0, outBeams);
        double incomingE = m_scattering->GetIncedentEnergy() * weight;
        if (ok)
            handlerPO->CacheBeams(outBeams, weight, D_ref, incomingE, cache.orientations[i]);
        else {
            cache.orientations[i].weight = weight;
            cache.orientations[i].incomingEnergy = incomingE;
        }
        if (m_mpiRank == 0) OutputProgress(nOrient, count, i, 0, timer, outBeams.size());
        outBeams.clear();
        ++count;
    }

    auto t_cache_end = std::chrono::high_resolution_clock::now();
    double cache_sec = std::chrono::duration<double>(t_cache_end - t_cache_start).count();

    EraseConsoleLine(60);
    if (m_mpiRank == 0) std::cout << "Phase 1 (tracing): 100%, "
              << cache.totalBeams() << " beams cached in " << cache_sec << " s" << std::endl;

    // Phase 2: Compute diffraction for all sizes
    auto t_diff_start = std::chrono::high_resolution_clock::now();

    std::vector<Arr2D> results_M;
    std::vector<double> results_energy;
    handlerPO->ComputeFromCache(cache, x_sizes, results_M, results_energy);

    auto t_diff_end = std::chrono::high_resolution_clock::now();
    double diff_sec = std::chrono::duration<double>(t_diff_end - t_diff_start).count();

    if (m_mpiRank == 0)
        std::cout << "Phase 2 (diffraction for " << x_sizes.size() << " sizes): "
                  << diff_sec << " s" << std::endl
                  << "Total: " << cache_sec + diff_sec << " s" << std::endl;

    // Write results for each size
    std::string baseName = m_resultDirName;
    Arr2D &origM = handlerPO->M;

    for (size_t s = 0; s < x_sizes.size(); ++s) {
        Arr2D savedM = origM;
        origM = results_M[s];
        m_incomingEnergy = results_energy[s];

        std::string sizeName = baseName + "_x" + std::to_string((int)x_sizes[s]);
        handlerPO->WriteMatricesToFile(sizeName, m_incomingEnergy);

        // Phi-averaged output
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

            for (int iZen = 0; iZen <= nZen; ++iZen) {
                matrix Msum(4, 4); Msum.Fill(0.0);
                double radZen = sphere.GetZenith(iZen);
                for (int iAz = 0; iAz <= nAz; ++iAz) {
                    double radAz = -iAz * sphere.azinuthStep;
                    matrix m = origM(iAz, iZen);
                    (*Lp)[1][1] = cos(2*radAz); (*Lp)[1][2] = sin(2*radAz);
                    (*Lp)[2][1] = -(*Lp)[1][2]; (*Lp)[2][2] = (*Lp)[1][1];
                    Msum += m * (*Lp);
                }
                Msum /= nAz;
                double _2PiDcos = sphere.Compute2PiDcos(iZen);
                outFile << std::endl << RadToDeg(radZen) << ' ' << _2PiDcos << ' ';
                outFile << Msum;
            }
            outFile.close();
        }
        origM = savedM;
    }

    if (m_mpiRank == 0)
        std::cout << "Results written to " << baseName << "_x*.dat" << std::endl;

    OutputStatisticsPO(timer, nOrient, m_resultDirName);
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

    // MPI: each rank processes a subset of orientations
    int myStart = m_mpiRank * nOrient / m_mpiSize;
    int myEnd = (m_mpiRank + 1) * nOrient / m_mpiSize;
    int myCount = myEnd - myStart;

    if (m_mpiRank == 0)
        std::cout << "Sobol: " << nOrient << " orientations"
                  << (m_mpiSize > 1 ? " (" + std::to_string(myCount) + " per rank, " + std::to_string(m_mpiSize) + " ranks)" : "")
                  << ", beta_sym=" << RadToDeg(betaSym) << " deg, gamma_sym="
                  << RadToDeg(gammaSym) << " deg" << std::endl;

    CalcTimer timer;
    timer.Start();
    if (m_mpiRank == 0) OutputStartTime(timer);

    double weight = 1.0 / nOrient;  // same weight regardless of MPI split
    m_handler->SetNormIndex(1);

    std::string dir;
    if (m_mpiRank == 0) dir = CreateFolder(m_resultDirName);
#ifdef USE_MPI
    if (m_mpiSize > 1) {
        // Broadcast dir name to all ranks
        int dirLen = dir.size();
        MPI_Bcast(&dirLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
        dir.resize(dirLen);
        MPI_Bcast(&dir[0], dirLen, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
#endif
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
    long long beamBudget = std::max(100LL, availMB / 2);
    int chunkSize = std::max(32, std::min(4096, std::min(myCount, (int)(beamBudget * 1024 / 350))));
    int nChunks = (myCount + chunkSize - 1) / chunkSize;

    if (m_mpiRank == 0)
        std::cerr << "Memory: " << availMB << " MB available, chunk="
                  << chunkSize << " orientations (" << nChunks << " chunks)" << std::endl;

    double phase1_total = 0, phase2_total = 0;
    std::vector<Beam> outBeams;
    long long count = 0;

    for (int chunk = 0; chunk < nChunks; ++chunk)
    {
        int iStart = myStart + chunk * chunkSize;
        int iEnd = std::min(iStart + chunkSize, myEnd);
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
            if (m_mpiRank == 0) OutputProgress(nOrient, count, iStart + i, 0, timer, outBeams.size());
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
            Arr2D localM_ns(nAz + 1, nZen + 1, 4, 4);
            localM_ns.ClearArr();
            std::vector<Arr2DC> localJ, localJ_ns;
            if (handlerPO->isCoh) {
                Arr2DC tmp(nAz + 1, nZen + 1, 2, 2);
                tmp.ClearArr();
                localJ.push_back(tmp);
                Arr2DC tmp2(nAz + 1, nZen + 1, 2, 2);
                tmp2.ClearArr();
                localJ_ns.push_back(tmp2);
            }

            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < thisChunk; ++i)
            {
                if (!chunkPrepared[i].beams.empty())
                    handlerPO->HandleBeamsToLocal(chunkPrepared[i], localM, localJ,
                                                   handlerPO->isCoh ? &localJ_ns : nullptr);
                if (handlerPO->isCoh && !localJ.empty()) {
                    double w = chunkPrepared[i].sinZenith;
                    HandlerPO::AddToMuellerLocal(localJ, w, localM, nAz, nZen);
                    HandlerPO::AddToMuellerLocal(localJ_ns, w, localM_ns, nAz, nZen);
                    localJ[0].ClearArr();
                    localJ_ns[0].ClearArr();
                }
            }

            #pragma omp critical
            {
                for (int p = 0; p <= nAz; ++p)
                    for (int t = 0; t <= nZen; ++t) {
                        handlerPO->M.insert(p, t, localM(p, t));
                        handlerPO->M_noshadow.insert(p, t, localM_ns(p, t));
                    }
            }
        }
        phase2_total += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tp2).count();

        chunkPrepared.clear();
        chunkPrepared.shrink_to_fit();
    }

    // MPI: reduce Mueller matrices from all ranks to rank 0
    MPI_ReduceMueller(handlerPO, nAz, nZen, m_incomingEnergy, m_mpiRank);

    EraseConsoleLine(60);
    if (m_mpiRank == 0) {
        std::cout << "Phase 1 (tracing): " << std::fixed
                  << std::setprecision(2) << phase1_total << " s" << std::endl;
        std::cout << "Phase 2 (diffraction, OpenMP): " << phase2_total << " s" << std::endl;
        std::cout << "Total: " << phase1_total + phase2_total << " s" << std::endl;

        m_handler->WriteTotalMatricesToFile(m_resultDirName);
        m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);

        // Write no-shadow Mueller
        std::swap(handlerPO->M, handlerPO->M_noshadow);
        std::string nsName = m_resultDirName + "_noshadow";
        handlerPO->WriteMatricesToFile(nsName, m_incomingEnergy);
        std::swap(handlerPO->M, handlerPO->M_noshadow);

        OutputStatisticsPO(timer, nOrient, m_resultDirName);
    }
}

void TracerPOTotal::TraceAdaptive(double eps, double betaSym, double gammaSym, int maxOrientOverride)
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

    // =========================================================================
    // Incremental adaptive: each iteration adds NEW orientations only.
    // Sobol property: first 2^m points ⊂ first 2^(m+1) points.
    // So doubling N means adding N new points (indices N..2N-1).
    //
    // Mueller is accumulated with weight=1 per batch, then after each
    // doubling: M_total = (M_old + M_new) / 2 = average of equal-size batches.
    // Since Mueller is additive, we store M_accumulated (sum of all batches)
    // and divide by number of batches when extracting results.
    //
    // Cost: N + N + N + ... = 2N (geometric series) instead of N+2N+4N = 7N.
    // Speedup vs restart: ~3.5× for convergence at the same N.
    // =========================================================================
    // Physics-based starting N: oldauto div16 formula
    // Δθ = 0.69 * λ / Dmax * (180/π), orient_step = Δθ / 3
    // N_beta = betaSym / orient_step / 16, N_gamma = gammaSym / orient_step / 16
    // N_start = N_beta × N_gamma rounded down to power of 2
    double Dmax = m_particle->MaximalDimention();
    double delta_deg = 0.69 * m_scattering->m_wave / Dmax * (180.0 / M_PI);
    double orient_step = delta_deg / 3.0;
    int nb16 = std::max(1, (int)(RadToDeg(betaSym) / orient_step / 16));
    int ng16 = std::max(1, (int)(RadToDeg(gammaSym) / orient_step / 16));
    int nStartRaw = nb16 * ng16;
    int nOrient = 1;
    while (nOrient * 2 <= nStartRaw) nOrient *= 2;
    if (nOrient < 64) nOrient = 64;
    std::cout << "  Start estimate (div16): " << nb16 << " x " << ng16
              << " = " << nStartRaw << " -> " << nOrient << " (power of 2)" << std::endl;

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
    // maxOrient from div1 estimate (full diffraction-limited grid)
    int nb1 = std::max(1, (int)(RadToDeg(betaSym) / orient_step));
    int ng1 = std::max(1, (int)(RadToDeg(gammaSym) / orient_step));
    int maxFromPhysics = nb1 * ng1;
    // Round down to power of 2
    int maxP2 = 1;
    while (maxP2 * 2 <= maxFromPhysics) maxP2 *= 2;

    int maxOrient;
    if (maxOrientOverride > 0) {
        maxOrient = 1;
        while (maxOrient < maxOrientOverride) maxOrient *= 2;
        if (maxOrient > maxOrientOverride) maxOrient /= 2;
        if (maxOrient < 64) maxOrient = 64;
    } else {
        // Use physics estimate, but cap by available memory
        int maxFromMem = std::max(1024, (int)(availMB_ad / 2 * 1024 / 350));
        int memP2 = 1; while (memP2 * 2 <= maxFromMem) memP2 *= 2;
        maxOrient = std::min(maxP2, memP2);
        if (maxOrient < 64) maxOrient = 64;
    }
    std::cout << "  Max estimate (div1): " << nb1 << " x " << ng1
              << " = " << maxFromPhysics << " -> " << maxP2 << " (power of 2)" << std::endl;
    std::cerr << "Adaptive: max orientations = " << maxOrient
              << " (" << availMB_ad << " MB available)" << std::endl;

    auto t_start = std::chrono::high_resolution_clock::now();

    // Generate ALL Sobol points up to maxOrient at once (deterministic)
    Sobol2D sobol_gen(42);
    std::vector<double> su_all, sv_all;
    sobol_gen.generate(maxOrient, su_all, sv_all);

    // Accumulator: M stores sum of ALL batches (not averaged yet)
    hp->M.ClearArr();
    hp->M_noshadow.ClearArr();
    hp->CleanJ();
    m_incomingEnergy = 0;
    hp->m_outputEnergy = 0;

    double prevM11_180 = 0, prevM11_90 = 0, prevM11_46 = 0, prevM11_22 = 0;
    double prevCsca = 0;
    int totalOrient = 0;
    int nBatches = 0;
    int convergedCount = 0;
    double ctrlAccum[4] = {0,0,0,0}; // accumulated M11 at 4 control angles

    for (int iter = 0; iter < 15; ++iter)
    {
        int batchStart = totalOrient;
        int batchEnd = nOrient;  // nOrient = target total after this iteration
        int batchSize = batchEnd - batchStart;

        if (batchSize <= 0) { nOrient *= 2; continue; }

        // Trace and diffract only the NEW orientations [batchStart..batchEnd)
        // Each batch uses weight = 1.0/batchSize (self-contained average)
        // After all batches: M_total = sum(M_batch_i) / nBatches

        // MPI: each rank processes a subset of this batch
        int myBatchStart = m_mpiRank * batchSize / m_mpiSize;
        int myBatchEnd = (m_mpiRank + 1) * batchSize / m_mpiSize;
        int myBatchSize = myBatchEnd - myBatchStart;

        double batchEnergy = 0;
        double weight = 1.0 / batchSize;  // global weight (not per-rank)
        std::vector<Beam> outBeams;

        // Phase 1 (sequential): trace and preprocess this rank's portion
        std::vector<PreparedOrientation> batchPrepared(myBatchSize);
        for (int i = 0; i < myBatchSize; ++i)
        {
            int idx = batchStart + myBatchStart + i;
            double beta  = acos(1.0 - (1.0 - cosBetaSym) * su_all[idx]);
            double gamma = gammaSym * sv_all[idx];
            m_particle->Rotate(beta, gamma, 0);
            if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
            bool ok = m_scattering->ScatterLight(0, 0, outBeams);

            if (ok)
                hp->PrepareBeams(outBeams, weight, batchPrepared[i]);
            else
                batchPrepared[i].sinZenith = weight;

            batchEnergy += m_scattering->GetIncedentEnergy() * weight;
            outBeams.clear();
        }

        // Phase 2: diffract ONLY at 4 control angles (not full grid!)
        // Full grid computed once after convergence.
        auto findThetaIdx = [&](double deg) -> int {
            double rad = DegToRad(deg);
            int best = 0; double bestD = 1e30;
            for (int jj = 0; jj <= nZen; ++jj) {
                double d = fabs(hp->m_sphere.GetZenith(jj) - rad);
                if (d < bestD) { bestD = d; best = jj; }
            }
            return best;
        };
        int ctrlIdx[4] = { findThetaIdx(22.0), findThetaIdx(46.0),
                           findThetaIdx(90.0), nZen };

        double batchCtrl[4] = {0,0,0,0};
        for (int i = 0; i < myBatchSize; ++i) {
            if (batchPrepared[i].beams.empty()) continue;
            double m11[4];
            hp->DiffractControlPoints(batchPrepared[i], ctrlIdx, 4, m11);
            for (int k = 0; k < 4; ++k) batchCtrl[k] += m11[k];
        }

        for (int k = 0; k < 4; ++k) ctrlAccum[k] += batchCtrl[k];

        // Free batch memory — NOT stored (final phase re-traces)
        batchPrepared.clear();
        batchPrepared.shrink_to_fit();

        // OLD full Phase 2 removed — replaced by fast control points above
        // Full diffraction done once after convergence via TraceFromSobol
        if (false) // dead code
        #pragma omp parallel
        {
            Arr2D localM(nAz + 1, nZen + 1, 4, 4);
            localM.ClearArr();
            Arr2D localM_ns(nAz + 1, nZen + 1, 4, 4);
            localM_ns.ClearArr();
            std::vector<Arr2DC> localJ, localJ_ns;
            if (hp->isCoh) {
                Arr2DC tmp(nAz+1, nZen+1, 2, 2); tmp.ClearArr();
                localJ.push_back(tmp);
                Arr2DC tmp2(nAz+1, nZen+1, 2, 2); tmp2.ClearArr();
                localJ_ns.push_back(tmp2);
            }

            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < myBatchSize; ++i)
            {
                if (!batchPrepared[i].beams.empty())
                    hp->HandleBeamsToLocal(batchPrepared[i], localM, localJ,
                                            hp->isCoh ? &localJ_ns : nullptr);
                if (hp->isCoh && !localJ.empty()) {
                    double w = batchPrepared[i].sinZenith;
                    HandlerPO::AddToMuellerLocal(localJ, w, localM, nAz, nZen);
                    HandlerPO::AddToMuellerLocal(localJ_ns, w, localM_ns, nAz, nZen);
                    localJ[0].ClearArr();
                    localJ_ns[0].ClearArr();
                }
            }

            #pragma omp critical
            {
                for (int p = 0; p <= nAz; ++p)
                    for (int t = 0; t <= nZen; ++t) {
                        hp->M.insert(p, t, localM(p, t));
                        hp->M_noshadow.insert(p, t, localM_ns(p, t));
                    }
            }
        }

        m_incomingEnergy += batchEnergy;
        totalOrient = batchEnd;
        nBatches++;

        // MPI: reduce only 4 control points + energy (not full Mueller!)
#ifdef USE_MPI
        if (m_mpiSize > 1)
        {
            double sbuf5[5] = {ctrlAccum[0], ctrlAccum[1], ctrlAccum[2], ctrlAccum[3], m_incomingEnergy};
            double rbuf5[5];
            MPI_Allreduce(sbuf5, rbuf5, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            ctrlAccum[0]=rbuf5[0]; ctrlAccum[1]=rbuf5[1]; ctrlAccum[2]=rbuf5[2]; ctrlAccum[3]=rbuf5[3];
            m_incomingEnergy = rbuf5[4];
        }
#endif

        // Control points from fast DiffractControlPoints (not full grid)
        double Csca = (nBatches > 0) ? m_incomingEnergy / nBatches : 0;
        double M11_22  = ((nBatches > 0) ? ctrlAccum[0] / nBatches : 0);
        double M11_46  = ((nBatches > 0) ? ctrlAccum[1] / nBatches : 0);
        double M11_90  = ((nBatches > 0) ? ctrlAccum[2] / nBatches : 0);
        double M11_180 = ((nBatches > 0) ? ctrlAccum[3] / nBatches : 0);

        double dCsca  = (prevCsca > 0)    ? fabs(Csca - prevCsca) / prevCsca : 1.0;
        double dM22   = (prevM11_22 > 0)  ? fabs(M11_22 - prevM11_22) / prevM11_22 : 1.0;
        double dM46   = (prevM11_46 > 0)  ? fabs(M11_46 - prevM11_46) / prevM11_46 : 1.0;
        double dM90   = (prevM11_90 > 0)  ? fabs(M11_90 - prevM11_90) / prevM11_90 : 1.0;
        double dM180  = (prevM11_180 > 0) ? fabs(M11_180 - prevM11_180) / fabs(prevM11_180) : 1.0;
        double dMax = std::max({dCsca, dM22, dM46, dM90, dM180});

        std::cout << std::fixed << std::setprecision(2);
        if (m_mpiRank == 0) std::cout << "  N=" << totalOrient << " (+" << batchSize << ")"
                  << "  dQ=" << dCsca*100 << "%"
                  << "  d22=" << dM22*100 << "%"
                  << "  d46=" << dM46*100 << "%"
                  << "  d90=" << dM90*100 << "%"
                  << "  d180=" << dM180*100 << "%"
                  << "  max=" << dMax*100 << "%"
                  << std::endl;

        bool all_ok = (dMax < eps && iter > 0);

        if (all_ok)
            convergedCount++;
        else
            convergedCount = 0;

        if (convergedCount >= 2)
        {
            if (m_mpiRank == 0) std::cout << "Converged at N=" << totalOrient
                      << " (all 5 controls within " << eps*100
                      << "% for 2 consecutive steps)" << std::endl;
            break;
        }

        prevCsca = Csca;
        prevM11_22 = M11_22;
        prevM11_46 = M11_46;
        prevM11_90 = M11_90;
        prevM11_180 = M11_180;
        nOrient *= 2;
        if (nOrient > maxOrient)
        {
            if (m_mpiRank == 0) std::cout << "WARNING: Max orientations reached (N=" << totalOrient
                      << ", limit=" << maxOrient << "). Target accuracy "
                      << eps*100 << "% may not be achieved." << std::endl;
            std::cout << "  To improve: use --maxorient " << maxOrient*2 << std::endl;
            break;
        }
    }

    // =================================================================
    // FINAL PHASE: full diffraction via TraceFromSobol (re-traces).
    // Sobol is deterministic → same orientations. Phase 1 = ~5% overhead.
    // Memory = O(chunkSize), not O(totalOrient). Chunked + OpenMP + MPI.
    // =================================================================
    if (m_mpiRank == 0)
        std::cout << "\nFinal: full diffraction with N=" << totalOrient
                  << " (re-tracing, chunked)" << std::endl;

    hp->M.ClearArr();
    hp->M_noshadow.ClearArr();
    hp->CleanJ();
    m_incomingEnergy = 0;
    hp->m_outputEnergy = 0;

    // TraceFromSobol handles file output, MPI reduce, and cleanup.
    // It re-traces all orientations (Sobol deterministic → same results)
    // using chunked streaming (O(chunkSize) memory).
    TraceFromSobol(totalOrient, betaSym, gammaSym);

    if (m_mpiRank == 0) {
        auto t_end = std::chrono::high_resolution_clock::now();
        double total_sec = std::chrono::duration<double>(t_end - t_start).count();
        std::cout << "Adaptive total time: " << std::fixed
                  << std::setprecision(1) << total_sec << " s" << std::endl;
    }
}

// =============================================================================
// TraceAutoFull: 3D sequential optimization n → N_phi → N_orient
// Step 1: Increase n until Q_sca converges (cheap: tracing only)
// Step 2: Increase N_phi until Q_sca converges (medium: new grid)
// Step 3: Increase N_orient until all control points converge (expensive)
// =============================================================================
void TracerPOTotal::TraceAutoFull(double eps, double betaSym, double gammaSym,
                                   int maxOrientOverride,
                                   Particle *particle, double wave,
                                   ScatteringRange &conus,
                                   HandlerPOTotal *handler)
{
    auto t_total_start = std::chrono::high_resolution_clock::now();
    HandlerPO *hp = dynamic_cast<HandlerPO*>(m_handler);
    if (!hp) { std::cerr << "Error: not HandlerPO" << std::endl; return; }

    int nAz = hp->m_sphere.nAzimuth;
    int nZen = hp->m_sphere.nZenith;

    std::cout << "===== AUTOFULL: 3D sequential optimization =====" << std::endl;
    std::cout << "Target accuracy: " << eps*100 << "%" << std::endl;
    std::cout << std::endl;

    // =========================================================================
    // STEP 1: Find optimal n (reflections)
    // Run with small N_orient=128 and current N_phi, vary n
    // =========================================================================
    std::cout << "--- Step 1: n convergence (N_orient=128, current N_phi) ---" << std::endl;

    int n_opt = 4;
    double prev_qsca_n = 0;
    int n_converged_count = 0;

    for (int n_test = 4; n_test <= 30; n_test += 2)
    {
        // Set reflection count
        m_scattering->SetMaxReflections(n_test); // hacky but works

        // Run quick computation
        hp->M.ClearArr(); hp->M_noshadow.ClearArr(); hp->CleanJ();
        m_incomingEnergy = 0; hp->m_outputEnergy = 0;

        Sobol2D sobol(42);
        std::vector<double> su, sv;
        sobol.generate(128, su, sv);
        double cosBetaSym = cos(betaSym);
        double weight = 1.0 / 128;
        std::vector<Beam> outBeams;

        for (int i = 0; i < 128; ++i)
        {
            double beta = acos(1.0 - (1.0 - cosBetaSym) * su[i]);
            double gamma = gammaSym * sv[i];
            m_particle->Rotate(beta, gamma, 0);
            if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
            m_scattering->ScatterLight(0, 0, outBeams);
            m_handler->HandleBeams(outBeams, weight);
            m_incomingEnergy += m_scattering->GetIncedentEnergy() * weight;
            outBeams.clear();
        }

        // Extract Q_sca
        double csca = 0;
        for (int j = 0; j <= nZen; ++j) {
            double m11_avg = 0;
            for (int p = 0; p <= nAz; ++p) m11_avg += hp->M(p, j)[0][0];
            m11_avg /= (nAz + 1);
            csca += m11_avg * hp->m_sphere.Compute2PiDcos(j);
        }
        double qsca = (m_incomingEnergy > 0) ? csca / m_incomingEnergy : 0;
        double dq = (prev_qsca_n > 0) ? fabs(qsca - prev_qsca_n) / prev_qsca_n : 1.0;

        std::cout << "  n=" << n_test << " Q_sca=" << std::fixed << std::setprecision(4)
                  << qsca << " dQ=" << std::setprecision(2) << dq*100 << "%" << std::endl;

        if (dq < eps && n_test > 4) n_converged_count++;
        else n_converged_count = 0;

        if (n_converged_count >= 2) {
            n_opt = n_test - 2; // use the value that first converged
            std::cout << "  → n converged at " << n_opt << std::endl;
            break;
        }
        prev_qsca_n = qsca;
        n_opt = n_test;
    }
    m_scattering->SetMaxReflections(n_opt);
    std::cout << std::endl;

    // =========================================================================
    // STEP 2: Find optimal N_phi
    // Run with N_orient=128, n=n_opt, vary N_phi
    // =========================================================================
    std::cout << "--- Step 2: N_phi convergence (n=" << n_opt << ", N_orient=128) ---" << std::endl;

    int phi_opt = 48;
    double prev_qsca_phi = 0;
    int phi_converged_count = 0;

    for (int phi_test = 36; phi_test <= 300; phi_test = (int)(phi_test * 1.5 / 6) * 6)
    {
        if (phi_test < 36) phi_test = 36;

        // Rebuild scattering sphere with new N_phi
        conus.nAzimuth = phi_test;
        conus.azinuthStep = 2.0 * M_PI / phi_test;
        hp->SetScatteringSphere(conus);
        hp->SetTracks(hp->GetTracks());

        int nAz2 = hp->m_sphere.nAzimuth;
        int nZen2 = hp->m_sphere.nZenith;

        hp->M.ClearArr(); hp->M_noshadow.ClearArr(); hp->CleanJ();
        m_incomingEnergy = 0; hp->m_outputEnergy = 0;

        Sobol2D sobol(42);
        std::vector<double> su, sv;
        sobol.generate(128, su, sv);
        double cosBetaSym = cos(betaSym);
        double weight = 1.0 / 128;
        std::vector<Beam> outBeams;

        for (int i = 0; i < 128; ++i)
        {
            double beta = acos(1.0 - (1.0 - cosBetaSym) * su[i]);
            double gamma = gammaSym * sv[i];
            m_particle->Rotate(beta, gamma, 0);
            if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
            m_scattering->ScatterLight(0, 0, outBeams);
            m_handler->HandleBeams(outBeams, weight);
            m_incomingEnergy += m_scattering->GetIncedentEnergy() * weight;
            outBeams.clear();
        }

        double csca = 0;
        for (int j = 0; j <= nZen2; ++j) {
            double m11_avg = 0;
            for (int p = 0; p <= nAz2; ++p) m11_avg += hp->M(p, j)[0][0];
            m11_avg /= (nAz2 + 1);
            csca += m11_avg * hp->m_sphere.Compute2PiDcos(j);
        }
        double qsca = (m_incomingEnergy > 0) ? csca / m_incomingEnergy : 0;
        double dq = (prev_qsca_phi > 0) ? fabs(qsca - prev_qsca_phi) / prev_qsca_phi : 1.0;

        std::cout << "  N_phi=" << phi_test << " Q_sca=" << std::fixed << std::setprecision(4)
                  << qsca << " dQ=" << std::setprecision(2) << dq*100 << "%" << std::endl;

        if (dq < eps && phi_test > 36) phi_converged_count++;
        else phi_converged_count = 0;

        if (phi_converged_count >= 2) {
            phi_opt = phi_test - (int)((phi_test - phi_test/1.5) + 0.5); // approximate previous
            phi_opt = ((phi_opt + 5) / 6) * 6;
            std::cout << "  → N_phi converged at " << phi_opt << std::endl;
            break;
        }
        prev_qsca_phi = qsca;
        phi_opt = phi_test;
    }

    // Set final N_phi
    conus.nAzimuth = phi_opt;
    conus.azinuthStep = 2.0 * M_PI / phi_opt;
    hp->SetScatteringSphere(conus);
    hp->SetTracks(hp->GetTracks());
    std::cout << std::endl;

    // =========================================================================
    // STEP 3: Adaptive N_orient (reuse TraceAdaptive with found n, N_phi)
    // =========================================================================
    std::cout << "--- Step 3: N_orient convergence (n=" << n_opt
              << ", N_phi=" << phi_opt << ") ---" << std::endl;

    // Clean and run adaptive
    hp->M.ClearArr(); hp->M_noshadow.ClearArr(); hp->CleanJ();
    m_incomingEnergy = 0; hp->m_outputEnergy = 0;

    TraceAdaptive(eps, betaSym, gammaSym, maxOrientOverride);

    auto t_total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double>(t_total_end - t_total_start).count();
    std::cout << std::endl << "===== AUTOFULL TOTAL: " << std::fixed
              << std::setprecision(1) << total_time << " s =====" << std::endl;
    std::cout << "Final: n=" << n_opt << ", N_phi=" << phi_opt << std::endl;
}
