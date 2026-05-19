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
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <future>
#include <sys/stat.h>

// ---- Checkpoint save/load for resume after crash ----
static void SaveCheckpoint(const std::string &path, const Arr2D &M, const Arr2D &M_ns,
                            double energy, double outputEnergy,
                            int completedOrient, int totalOrient,
                            unsigned long paramHash, int nAz, int nZen)
{
    std::ofstream f(path, std::ios::binary);
    if (!f.is_open()) return;
    uint32_t magic = 0x4D425344; // "MBSD"
    f.write((char*)&magic, 4);
    f.write((char*)&paramHash, sizeof(paramHash));
    f.write((char*)&completedOrient, sizeof(completedOrient));
    f.write((char*)&totalOrient, sizeof(totalOrient));
    f.write((char*)&energy, sizeof(energy));
    f.write((char*)&outputEnergy, sizeof(outputEnergy));
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
                            double &energy, double &outputEnergy,
                            int &completedOrient, int totalOrient,
                            unsigned long paramHash)
{
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) return false;
    uint32_t magic; f.read((char*)&magic, 4);
    if (magic != 0x4D425344) return false;
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
    f.read((char*)&outputEnergy, sizeof(outputEnergy));
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

static unsigned long HashParams(double wave, double ri_re, double ri_im, int nActs, int nOrient,
                                 double L, double D, int nAz, int nZen)
{
    unsigned long h = 0;
    auto mix = [&](double v) { h ^= std::hash<double>{}(v) + 0x9e3779b9 + (h<<6) + (h>>2); };
    auto mixi = [&](int v) { h ^= std::hash<int>{}(v) + 0x9e3779b9 + (h<<6) + (h>>2); };
    mix(wave); mix(ri_re); mix(ri_im); mixi(nActs); mixi(nOrient); mix(L); mix(D); mixi(nAz); mixi(nZen);
    return h;
}
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace std;
using ::complex;

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

// Helper: MPI reduce Arr2D + energy counters, rank 0 gets result
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

    // Reduce output energy used for absorption/extinction accounting.
    double totalOutputEnergy = 0;
    MPI_Reduce(&hp->m_outputEnergy, &totalOutputEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) hp->m_outputEnergy = totalOutputEnergy;
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
    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) {
        std::cerr << "Error: handler is not HandlerPO in TraceRandom" << std::endl;
        return;
    }

    m_incomingEnergy = 0;
    handlerPO->m_outputEnergy = 0;
    handlerPO->M.ClearArr();
    handlerPO->M_noshadow.ClearArr();
    handlerPO->CleanJ();

    int betaNorm = (m_symmetry.beta < M_PI_2+FLT_EPSILON && m_symmetry.beta > M_PI_2-FLT_EPSILON) ? 1 : 2;
    double normGamma = gammaRange.number * betaNorm; // MBS-raw: normalize by gammaRange.number

    // --coh_orient: legacy mode — HandleBeams accumulates Jones coherently
    // across ALL orientations, then single AddToMueller at end.
    // No chunking, no OpenMP (matches old MBS-raw exactly).
    if (m_cohOrient) {
        m_handler->SetNormIndex(normGamma);
        vector<Beam> outBeams;
        for (int i = 0; i <= betaRange.number; ++i) {
            double beta = betaRange.min + i * betaRange.step;
            double dcos;
            CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normGamma, dcos);
            const bool pole = (fabs(beta) <= FLT_EPSILON
                               || fabs(beta - M_PI) <= FLT_EPSILON);
            const bool fastPole = pole && m_fastPoleGamma;
            const int gammaCount = fastPole ? 1 : gammaRange.number;
            const double gammaWeight = fastPole ? dcos * gammaRange.number : dcos;
            m_handler->SetSinZenith(gammaWeight);
            for (int j = 0; j < gammaCount; ++j) {
                double gamma = gammaRange.min + (j + 0.5) * gammaRange.step;
                m_particle->Rotate(beta, gamma, 0);
                if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
                m_scattering->ScatterLight(0, 0, outBeams);
                m_handler->HandleBeams(outBeams, gammaWeight);
                m_incomingEnergy += m_scattering->GetIncedentEnergy() * gammaWeight;
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

    int nGamma = gammaRange.number; // MBS-raw: same as gammaRange.number
    int nBeta = betaRange.number + 1;
    int nOrientations = nBeta * nGamma;

    if (m_mpiRank == 0)
        std::cout << "Random grid: " << nBeta << " x " << nGamma
                  << " = " << nOrientations << " orientations" << std::endl;

    CalcTimer timer;
    timer.Start();
    if (m_mpiRank == 0) OutputStartTime(timer);

    m_handler->SetNormIndex(normGamma);

    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;
    int nThreads = 1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        nThreads = omp_get_num_threads();
    }
#endif
    int parallelTraceMinGamma = std::max(128, 8 * nThreads);
    const char *traceMinEnv = std::getenv("MBS_TRACE_MIN_GAMMA");
    if (traceMinEnv && *traceMinEnv)
    {
        char *end = nullptr;
        long value = std::strtol(traceMinEnv, &end, 10);
        if (end && *end == '\0' && value > 0)
            parallelTraceMinGamma = (int)value;
    }
    if (m_mpiRank == 0)
    {
        std::ostringstream line;
        line << "Oldauto/random Phase 1: ";
        if (nThreads > 1)
        {
            line << "gamma blocks < " << parallelTraceMinGamma
                 << " trace sequentially; larger blocks trace in parallel ("
                 << nThreads << " threads)";
        }
        else
        {
            line << "sequential tracing (1 thread)";
        }
        std::cout << line.str() << std::endl;
        AppendTextLog(line.str() + "\n");
    }

    // --save_betas: create output directory
    std::string betaDir;
    if (m_saveBetas && m_mpiRank == 0) {
        betaDir = m_resultDirName + "_betas";
        mkdir(betaDir.c_str(), 0755);
    }

    double phase1_total = 0, phase2_total = 0;
    long long count = 0;

    struct BetaBlock
    {
        int ib = 0;
        int gammaCount = 0;
        double beta = 0.0;
        double dcos = 0.0;
        double phase1 = 0.0;
        std::vector<PreparedOrientation> prepared;
        std::vector<double> energies;
        std::vector<double> outputEnergies;
        std::vector<int> beamCounts;
    };

    auto prepareBetaBlock = [&](int ib) {
        BetaBlock block;
        block.ib = ib;
        block.beta = betaRange.min + ib * betaRange.step;
        CalcCsBeta(betaNorm, block.beta, betaRange, gammaRange, normGamma, block.dcos);
        const bool pole = (fabs(block.beta) <= FLT_EPSILON
                           || fabs(block.beta - M_PI) <= FLT_EPSILON);
        const bool fastPole = pole && m_fastPoleGamma;
        block.gammaCount = fastPole ? 1 : nGamma;
        const double gammaWeight = fastPole ? block.dcos * nGamma : block.dcos;

        auto tp1 = std::chrono::high_resolution_clock::now();
        block.prepared.assign(block.gammaCount, PreparedOrientation());
        block.energies.assign(block.gammaCount, 0.0);
        block.outputEnergies.assign(block.gammaCount, 0.0);
        block.beamCounts.assign(block.gammaCount, 0);

        if (nThreads > 1 && block.gammaCount >= parallelTraceMinGamma)
        {
            #pragma omp parallel
            {
                Particle localParticle = *m_particle;
                Scattering *localScatter = m_scattering->CloneFor(&localParticle, &m_incidentLight);
                localScatter->m_wave = m_scattering->m_wave;
                localScatter->restriction = m_scattering->restriction;

                HandlerPO localHandler(&localParticle, &m_incidentLight,
                                       handlerPO->nTheta, m_scattering->m_wave);
                localHandler.ConfigureForThreadLocalPrepare(*handlerPO, localScatter);

                std::vector<Beam> localBeams;

                #pragma omp for schedule(dynamic, 1)
                for (int jj = 0; jj < block.gammaCount; ++jj)
                {
                    double gamma = gammaRange.min + (jj + 0.5) * gammaRange.step;
                    localParticle.Rotate(block.beta, gamma, 0);
                    if (!shadowOff) localScatter->FormShadowBeam(localBeams);
                    bool ok = localScatter->ScatterLight(0, 0, localBeams);
                    block.beamCounts[jj] = (int)localBeams.size();
                    if (ok)
                    {
                        double beforeOutput = localHandler.m_outputEnergy;
                        localHandler.PrepareBeams(localBeams, gammaWeight, block.prepared[jj]);
                        block.outputEnergies[jj] = localHandler.m_outputEnergy - beforeOutput;
                    }
                    else
                    {
                        block.prepared[jj].sinZenith = gammaWeight;
                    }
                    block.energies[jj] = localScatter->GetIncedentEnergy() * gammaWeight;
                    localBeams.clear();
                }

                delete localScatter;
            }
        }
        else
        {
            Particle localParticle = *m_particle;
            Scattering *localScatter = m_scattering->CloneFor(&localParticle, &m_incidentLight);
            localScatter->m_wave = m_scattering->m_wave;
            localScatter->restriction = m_scattering->restriction;

            HandlerPO localHandler(&localParticle, &m_incidentLight,
                                   handlerPO->nTheta, m_scattering->m_wave);
            localHandler.ConfigureForThreadLocalPrepare(*handlerPO, localScatter);

            std::vector<Beam> localBeams;
            for (int jj = 0; jj < block.gammaCount; ++jj)
            {
                double gamma = gammaRange.min + (jj + 0.5) * gammaRange.step;
                localParticle.Rotate(block.beta, gamma, 0);
                if (!shadowOff) localScatter->FormShadowBeam(localBeams);
                bool ok = localScatter->ScatterLight(0, 0, localBeams);
                block.beamCounts[jj] = (int)localBeams.size();
                if (ok)
                {
                    double beforeOutput = localHandler.m_outputEnergy;
                    localHandler.PrepareBeams(localBeams, gammaWeight, block.prepared[jj]);
                    block.outputEnergies[jj] = localHandler.m_outputEnergy - beforeOutput;
                }
                else
                {
                    block.prepared[jj].sinZenith = gammaWeight;
                }
                block.energies[jj] = localScatter->GetIncedentEnergy() * gammaWeight;
                localBeams.clear();
            }

            delete localScatter;
        }

        block.phase1 = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - tp1).count();
        return block;
    };

    const bool overlapCpuGpu = handlerPO->IsGpuEnabled();
    std::future<BetaBlock> nextBeta;
    if (overlapCpuGpu && nBeta > 0)
        nextBeta = std::async(std::launch::async, prepareBetaBlock, 0);

    // Process beta-by-beta: each beta = one chunk of nGamma orientations
    for (int ib = 0; ib < nBeta; ++ib)
    {
        BetaBlock block;
        if (overlapCpuGpu)
        {
            block = nextBeta.get();
            if (ib + 1 < nBeta)
                nextBeta = std::async(std::launch::async, prepareBetaBlock, ib + 1);
        }
        else
        {
            block = prepareBetaBlock(ib);
        }

        const int gammaCount = block.gammaCount;
        std::vector<PreparedOrientation> &chunkPrepared = block.prepared;
        const double beta = block.beta;
        const double dcos = block.dcos;

        for (int jj = 0; jj < gammaCount; ++jj)
        {
            m_incomingEnergy += block.energies[jj];
            handlerPO->m_outputEnergy += block.outputEnergies[jj];
            ++count;
            if (m_mpiRank == 0)
                OutputProgress(nOrientations, count, ib*nGamma+jj, 0, timer, block.beamCounts[jj]);
        }
        phase1_total += block.phase1;

        // Phase 2: diffraction (OpenMP parallel over gamma orientations)
        auto tp2 = std::chrono::high_resolution_clock::now();
        Arr2D betaM(nAz+1, nZen+1, 4, 4); betaM.ClearArr();
        Arr2D betaM_ns(nAz+1, nZen+1, 4, 4); betaM_ns.ClearArr();

        if (handlerPO->IsGpuEnabled())
        {
            Arr2D localM(nAz+1, nZen+1, 4, 4); localM.ClearArr();
            Arr2D localM_ns(nAz+1, nZen+1, 4, 4); localM_ns.ClearArr();
            for (int gpuStart = 0; gpuStart < gammaCount; )
            {
                int gpuBatchSize = handlerPO->SelectGpuOrientationBatchSize(
                    chunkPrepared, gpuStart, gammaCount - gpuStart);
                int gpuEnd = std::min(gpuStart + gpuBatchSize, gammaCount);
                bool ok = handlerPO->IsFftEnabled()
                    ? handlerPO->HandleOrientationsToLocalGpuFftPhi(
                        chunkPrepared, gpuStart, gpuEnd - gpuStart,
                        localM, localM_ns)
                    : handlerPO->HandleOrientationsToLocalGpu(
                        chunkPrepared, gpuStart, gpuEnd - gpuStart,
                        localM, localM_ns);
                if (!ok)
                {
                    std::cerr << "ERROR: --gpu requested but GPU diffraction backend "
                              << "could not process oldauto/random beta block." << std::endl;
                    throw std::runtime_error("GPU diffraction backend failed");
                }
                gpuStart = gpuEnd;
            }
            for (int p=0;p<nAz;++p) for (int t=0;t<=nZen;++t) {
                betaM.insert(p,t,localM(p,t));
                betaM_ns.insert(p,t,localM_ns(p,t));
            }
        }
        else
        {
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
                for (int i = 0; i < gammaCount; ++i) {
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
                { for (int p=0;p<nAz;++p) for (int t=0;t<=nZen;++t) {
                    betaM.insert(p,t,localM(p,t));
                    betaM_ns.insert(p,t,localM_ns(p,t));
                } }
            }
        }
        phase2_total += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tp2).count();
        chunkPrepared.clear(); chunkPrepared.shrink_to_fit();

        // Accumulate into global Mueller
        for (int p=0;p<nAz;++p) for (int t=0;t<=nZen;++t) {
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
                for (int iAz = 0; iAz < nAz; ++iAz) {
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
                for (int iAz = 0; iAz < nAz; ++iAz) {
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
        std::ostringstream line;
        line << "Sequential tracing completed: " << count << "/" << nOrientations
             << " trace calls, phase1=" << std::fixed << std::setprecision(2)
             << phase1_total << " s";
        std::cout << line.str() << std::endl;
        AppendTextLog(line.str() + "\n");
        std::cout << "Phase 1 (tracing): " << std::fixed << std::setprecision(2) << phase1_total << " s" << std::endl;
        std::cout << "Phase 2 (diffraction, "
                  << (handlerPO->IsGpuEnabled() ? "CUDA" : "OpenMP")
                  << "): " << phase2_total << " s" << std::endl;
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

    const double betaSpan = betaRange.max - betaRange.min;
    const double betaIntegral = cos(betaRange.min) - cos(betaRange.max);
    const double betaWeightNorm = (fabs(betaIntegral) > DBL_EPSILON)
        ? betaSpan / (betaIntegral * nOrientations)
        : 0.0;

    for (int i = 0; i < nOrientations; ++i)
    {
        double beta = betaRange.min + RandomDouble(0, 1) * betaSpan;
        double gamma = gammaRange.min + RandomDouble(0, 1) * (gammaRange.max - gammaRange.min);
        orientations.push_back({beta, gamma});
        weights.push_back(sin(beta) * betaWeightNorm);
    }

    std::cout << "Monte Carlo: " << nOrientations << " orientations" << std::endl;

    CalcTimer timer;
    timer.Start();
    if (m_mpiRank == 0) OutputStartTime(timer);

    m_handler->SetNormIndex(1);

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) { std::cerr << "Error: not HandlerPO" << std::endl; return; }

    m_incomingEnergy = 0;
    handlerPO->m_outputEnergy = 0;
    handlerPO->M.ClearArr();
    handlerPO->M_noshadow.ClearArr();
    handlerPO->CleanJ();

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

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO)
    {
        std::cerr << "Error! Handler is not HandlerPO in TraceFromFile" << std::endl;
        throw std::exception();
    }

    m_incomingEnergy = 0;
    handlerPO->m_outputEnergy = 0;
    handlerPO->M.ClearArr();
    handlerPO->M_noshadow.ClearArr();
    handlerPO->CleanJ();

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

    // Checkpoint is opt-in because frequent binary dumps add avoidable I/O
    // overhead to ordinary orientation-file runs.
    std::string ckptPath = m_resultDirName + "_checkpoint.bin";
    const ::complex ri = m_particle->GetRefractiveIndex();
    unsigned long paramHash = HashParams(m_scattering->m_wave,
        real(ri), imag(ri), m_scattering->GetMaxReflections(),
        nOrientations, m_particle->MaximalDimention(), 0, nAz, nZen);
    for (const auto &bg : orientations)
    {
        paramHash ^= std::hash<double>{}(bg.first)
            + 0x9e3779b9 + (paramHash << 6) + (paramHash >> 2);
        paramHash ^= std::hash<double>{}(bg.second)
            + 0x9e3779b9 + (paramHash << 6) + (paramHash >> 2);
    }
    int resumeChunk = 0;
    if (m_enableCheckpoint)
    {
        int completedOrient = 0;
        if (LoadCheckpoint(ckptPath, handlerPO->M, handlerPO->M_noshadow,
                           m_incomingEnergy, handlerPO->m_outputEnergy,
                           completedOrient, nOrientations, paramHash))
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
            if (m_mpiRank == 0) OutputProgress(nOrientations, count + 1, iStart + i, 0, timer, outBeams.size());
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
                for (int p = 0; p < nAz; ++p)
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

        // Save checkpoint after each chunk.
        if (m_enableCheckpoint && m_mpiRank == 0) {
            SaveCheckpoint(ckptPath, handlerPO->M, handlerPO->M_noshadow,
                           m_incomingEnergy, handlerPO->m_outputEnergy,
                           iEnd, nOrientations, paramHash, nAz+1, nZen+1);
        }
    }

    EraseConsoleLine(60);
    if (m_mpiRank == 0)
    {
        std::ostringstream line;
        line << "Sequential tracing completed: " << count << "/" << nOrientations
             << " trace calls, phase1=" << std::fixed << std::setprecision(2)
             << phase1_total << " s";
        std::cout << line.str() << std::endl;
        AppendTextLog(line.str() + "\n");
    }
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

    // Remove checkpoint on successful completion.
    if (m_enableCheckpoint)
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

        OutputProgress(nOrientations, count + 1, i, 0, timer, outBeams.size());
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
    std::string baseName = m_resultDirName;

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

                for (int iAz = 0; iAz < nAz; ++iAz)
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
    // Simple approach: for each size, resize particle and run TraceFromSobol.
    // Uses the fully optimized HandleBeamsToLocal path (batched sincos, OpenMP).
    // Tracing is fast (<1s), dominated by diffraction which is already optimal.

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) {
        std::cerr << "Error! Handler is not HandlerPO in TraceSobolMultiSize" << std::endl;
        throw std::exception();
    }

    double D_ref = m_particle->MaximalDimention();
    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    std::string baseName = m_resultDirName;

    if (m_mpiRank == 0)
        std::cout << "MultiSize: " << x_sizes.size() << " sizes, "
                  << nOrient << " Sobol orientations each" << std::endl;

    auto t_total_start = std::chrono::high_resolution_clock::now();

    for (size_t s = 0; s < x_sizes.size(); ++s)
    {
        // Scale particle to this size
        double D_target = x_sizes[s] / x_ref * D_ref;
        m_particle->Resize(D_target);

        if (m_mpiRank == 0)
            std::cout << "  Size " << (s+1) << "/" << x_sizes.size()
                      << ": x=" << (int)x_sizes[s]
                      << " (D=" << m_particle->MaximalDimention() << ")" << std::endl;

        // Reset Mueller arrays
        handlerPO->M = Arr2D(nAz+1, nZen+1, 4, 4);
        handlerPO->M_noshadow = Arr2D(nAz+1, nZen+1, 4, 4);
        m_incomingEnergy = 0;

        // Use fully optimized TraceFromSobol (batched sincos, parallel Phase 1+2)
        // Save/restore m_resultDirName since TraceFromSobol writes output
        std::string savedName = m_resultDirName;
        std::string sizeName = baseName + "_x" + std::to_string((int)x_sizes[s]);
        m_resultDirName = sizeName;
        TraceFromSobol(nOrient, betaSym, gammaSym);
        m_resultDirName = savedName;
    }

    // Restore particle to original size
    m_particle->Resize(D_ref);

    auto t_total_end = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(t_total_end - t_total_start).count();

    if (m_mpiRank == 0)
        std::cout << "MultiSize total: " << total_sec << " s for "
                  << x_sizes.size() << " sizes" << std::endl;
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

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO)
    {
        std::cerr << "Error! Handler is not HandlerPO in TraceFromSobol" << std::endl;
        throw std::exception();
    }

    m_incomingEnergy = 0;
    handlerPO->m_outputEnergy = 0;
    handlerPO->M.ClearArr();
    handlerPO->M_noshadow.ClearArr();
    handlerPO->CleanJ();

    // Phase 1 traces and preprocesses orientations in parallel using
    // thread-local Particle/Scattering/HandlerPO state. Phase 2 then performs
    // the diffraction loops in parallel with local Mueller accumulators.
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
    if (m_sobolChunkSize > 0)
        chunkSize = std::max(1, std::min(chunkSize, m_sobolChunkSize));
    int nChunks = (myCount + chunkSize - 1) / chunkSize;
    int nThreads = 1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        nThreads = omp_get_num_threads();
    }
#endif

    if (m_mpiRank == 0)
    {
        std::ostringstream log;
        log << "Memory: " << availMB << " MB available, chunk="
            << chunkSize << " orientations (" << nChunks << " chunks), "
            << nThreads << " threads";
        std::cerr << log.str() << std::endl;
        AppendTextLog(log.str() + "\n");
    }

    double phase1_total = 0, phase2_total = 0;
    std::vector<Beam> outBeams;
    long long count = 0;

    for (int chunk = 0; chunk < nChunks; ++chunk)
    {
        int iStart = myStart + chunk * chunkSize;
        int iEnd = std::min(iStart + chunkSize, myEnd);
        int thisChunk = iEnd - iStart;

        // Phase 1: trace and preprocess this chunk in parallel. Each thread
        // has its own Particle, Scattering and HandlerPO scratch state.
        auto tp1 = std::chrono::high_resolution_clock::now();
        std::vector<PreparedOrientation> chunkPrepared(thisChunk);
        std::vector<double> chunkEnergies(thisChunk, 0);
        std::vector<double> chunkOutputEnergies(thisChunk, 0);

        #pragma omp parallel
        {
            // Thread-local copies of all mutable tracing/preprocessing state.
            Particle localParticle = *m_particle;
            Scattering *localScatter = m_scattering->CloneFor(&localParticle, &m_incidentLight);
            localScatter->m_wave = m_scattering->m_wave;
            localScatter->restriction = m_scattering->restriction;

            HandlerPO localHandler(&localParticle, &m_incidentLight,
                                   handlerPO->nTheta, m_scattering->m_wave);
            localHandler.ConfigureForThreadLocalPrepare(*handlerPO, localScatter);

            std::vector<Beam> localBeams;

            #pragma omp for schedule(dynamic, 4)
            for (int i = 0; i < thisChunk; ++i)
            {
                int idx = iStart + i;
                localParticle.Rotate(orientations[idx].first, orientations[idx].second, 0);
                if (!shadowOff) localScatter->FormShadowBeam(localBeams);
                bool ok = localScatter->ScatterLight(0, 0, localBeams);
                if (ok)
                {
                    double beforeOutput = localHandler.m_outputEnergy;
                    localHandler.PrepareBeams(localBeams, weight, chunkPrepared[i]);
                    chunkOutputEnergies[i] = localHandler.m_outputEnergy - beforeOutput;
                }
                else
                {
                    chunkPrepared[i].sinZenith = weight;
                }
                chunkEnergies[i] = localScatter->GetIncedentEnergy() * weight;
                localBeams.clear();
            }

            delete localScatter;
        }

        // Accumulate scalar counters sequentially after the parallel region.
        for (int i = 0; i < thisChunk; ++i) {
            m_incomingEnergy += chunkEnergies[i];
            handlerPO->m_outputEnergy += chunkOutputEnergies[i];
            count++;
        }
        phase1_total += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tp1).count();

        // Phase 2: parallel diffraction for this chunk
        auto tp2 = std::chrono::high_resolution_clock::now();

        if (handlerPO->IsGpuEnabled())
        {
            Arr2D localM(nAz + 1, nZen + 1, 4, 4);
            localM.ClearArr();
            Arr2D localM_ns(nAz + 1, nZen + 1, 4, 4);
            localM_ns.ClearArr();

            for (int gpuStart = 0; gpuStart < thisChunk; )
            {
                int gpuBatchSize = handlerPO->SelectGpuOrientationBatchSize(
                    chunkPrepared, gpuStart, thisChunk - gpuStart);
                int gpuEnd = std::min(gpuStart + gpuBatchSize, thisChunk);
                bool ok = handlerPO->IsFftEnabled()
                    ? handlerPO->HandleOrientationsToLocalGpuFftPhi(
                        chunkPrepared, gpuStart, gpuEnd - gpuStart,
                        localM, localM_ns)
                    : handlerPO->HandleOrientationsToLocalGpu(
                        chunkPrepared, gpuStart, gpuEnd - gpuStart,
                        localM, localM_ns);
                if (!ok)
                {
                    std::cerr << "ERROR: --gpu requested but GPU diffraction backend "
                              << "could not process this chunk." << std::endl;
                    throw std::runtime_error("GPU diffraction backend failed");
                }
                gpuStart = gpuEnd;
            }

            for (int p = 0; p < nAz; ++p)
                for (int t = 0; t <= nZen; ++t)
                {
                    handlerPO->M.insert(p, t, localM(p, t));
                    handlerPO->M_noshadow.insert(p, t, localM_ns(p, t));
                }
        }
        else
        {
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
                    for (int p = 0; p < nAz; ++p)
                        for (int t = 0; t <= nZen; ++t) {
                            handlerPO->M.insert(p, t, localM(p, t));
                            handlerPO->M_noshadow.insert(p, t, localM_ns(p, t));
                        }
                }
            }
        }
        phase2_total += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tp2).count();

        if (m_mpiRank == 0)
            OutputProgress(nOrient, count, iEnd - 1, chunk + 1, timer, -1);

        chunkPrepared.clear();
        chunkPrepared.shrink_to_fit();
    }

    // MPI: reduce Mueller matrices from all ranks to rank 0
    MPI_ReduceMueller(handlerPO, nAz, nZen, m_incomingEnergy, m_mpiRank);

    EraseConsoleLine(60);
    if (m_mpiRank == 0) {
        std::cout << "Phase 1 (tracing): " << std::fixed
                  << std::setprecision(2) << phase1_total << " s" << std::endl;
        std::cout << "Phase 2 (diffraction, "
                  << (handlerPO->IsGpuEnabled() ? "CUDA" : "OpenMP")
                  << "): " << phase2_total << " s" << std::endl;
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

void TracerPOTotal::TraceAdaptiveTheta(int nOrient, double betaSym, double gammaSym,
                                        double eps, int maxDepth, bool gridOnly)
{
    HandlerPO *hp = dynamic_cast<HandlerPO*>(m_handler);
    if (!hp) { TraceFromSobol(nOrient, betaSym, gammaSym); return; }

    // Step 1: Trace and cache beams
    Sobol2D sobol(42);
    std::vector<double> su, sv;
    sobol.generate(nOrient, su, sv);
    double cosBetaSym = cos(betaSym);

    double weight = 1.0 / nOrient;
    std::vector<Beam> outBeams;

    // Cache all PreparedOrientations (for DiffractAtThetas)
    std::vector<PreparedOrientation> allPrepared(nOrient);
    CalcTimer timer; timer.Start();
    if (m_mpiRank == 0) {
        std::cout << "AdaptiveTheta: tracing " << nOrient << " orientations..." << std::endl;
        OutputStartTime(timer);
    }

    for (int i = 0; i < nOrient; ++i) {
        double beta = acos(1.0 - (1.0 - cosBetaSym) * su[i]);
        double gamma = gammaSym * sv[i];
        m_particle->Rotate(beta, gamma, 0);
        if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
        bool ok = m_scattering->ScatterLight(0, 0, outBeams);
        if (ok) hp->PrepareBeams(outBeams, weight, allPrepared[i]);
        else    allPrepared[i].sinZenith = weight;
        m_incomingEnergy += m_scattering->GetIncedentEnergy() * weight;
        outBeams.clear();
    }

    auto t_trace_end = std::chrono::high_resolution_clock::now();
    if (m_mpiRank == 0) std::cout << "Tracing done." << std::endl;

    // Step 2: Compute M11 at arbitrary theta values (summed over all
    // orientations, phi-averaged). Batch theta candidates so each
    // orientation/beam/phi coefficient setup is reused across all candidates.
    double theta_eval_seconds = 0.0;
    bool thetaEvalUsedGpu = false;
    auto computeM11Batch = [&](const std::vector<double> &theta_rads) {
        auto t0 = std::chrono::high_resolution_clock::now();
        std::vector<double> totals(theta_rads.size(), 0.0);
        if (theta_rads.empty())
            return totals;

        const char *gpuAdaptiveThetaOff = std::getenv("MBS_GPU_ADAPTIVE_THETA_OFF");
        if (hp->IsGpuEnabled()
            && !(gpuAdaptiveThetaOff && gpuAdaptiveThetaOff[0] == '1' && gpuAdaptiveThetaOff[1] == '\0')
            && hp->DiffractThetasGpu(allPrepared, theta_rads.data(),
                                     (int)theta_rads.size(), totals))
        {
            thetaEvalUsedGpu = true;
            theta_eval_seconds += std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t0).count();
            return totals;
        }

#ifdef _OPENMP
        #pragma omp parallel
        {
            std::vector<double> localTotals(theta_rads.size(), 0.0);
            std::vector<double> m11_one(theta_rads.size(), 0.0);
            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < nOrient; ++i) {
                if (allPrepared[i].beams.empty()) continue;
                std::fill(m11_one.begin(), m11_one.end(), 0.0);
                hp->DiffractAtThetas(allPrepared[i], theta_rads.data(),
                                     (int)theta_rads.size(), m11_one.data());
                for (size_t k = 0; k < localTotals.size(); ++k)
                    localTotals[k] += m11_one[k];
            }
            #pragma omp critical
            {
                for (size_t k = 0; k < totals.size(); ++k)
                    totals[k] += localTotals[k];
            }
        }
#else
        std::vector<double> m11_one(theta_rads.size(), 0.0);
        for (int i = 0; i < nOrient; ++i) {
            if (allPrepared[i].beams.empty()) continue;
            std::fill(m11_one.begin(), m11_one.end(), 0.0);
            hp->DiffractAtThetas(allPrepared[i], theta_rads.data(),
                                 (int)theta_rads.size(), m11_one.data());
            for (size_t k = 0; k < totals.size(); ++k)
                totals[k] += m11_one[k];
        }
#endif
        theta_eval_seconds += std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - t0).count();
        return totals;
    };

    // Step 3: Adaptive bisection
    // Start with coarse grid
    double minTheta = 0.0, maxTheta = M_PI;
    double minStep = m_scattering->m_wave / (8.0 * m_particle->MaximalDimention()); // λ/(8D) = quarter-Nyquist

    // Initial points (30 uniform)
    int nInitial = 30;
    std::map<double, double> thetaM11; // theta_rad -> M11
    std::vector<double> initialThetas;
    initialThetas.reserve(nInitial + 1);
    for (int i = 0; i <= nInitial; ++i) {
        double t = minTheta + i * (maxTheta - minTheta) / nInitial;
        initialThetas.push_back(t);
    }
    std::vector<double> initialM11 = computeM11Batch(initialThetas);
    for (size_t i = 0; i < initialThetas.size(); ++i)
        thetaM11[initialThetas[i]] = initialM11[i];
    if (m_mpiRank == 0)
        std::cout << "Initial grid: " << thetaM11.size() << " points" << std::endl;

    // Recursive bisection
    int totalAdded = 0;
    for (int depth = 0; depth < maxDepth; ++depth) {
        struct ThetaCandidate
        {
            double theta;
            double interp;
        };
        std::vector<ThetaCandidate> candidates;
        auto it = thetaM11.begin();
        auto prev = it; ++it;
        while (it != thetaM11.end()) {
            double t1 = prev->first, t2 = it->first;
            double m1 = prev->second, m2 = it->second;
            double dt = t2 - t1;
            if (dt < minStep) { prev = it; ++it; continue; }

            double tmid = (t1 + t2) * 0.5;
            double m_interp = (m1 + m2) * 0.5;
            candidates.push_back({tmid, m_interp});
            prev = it; ++it;
        }

        std::vector<double> toEval;
        toEval.reserve(candidates.size());
        for (const ThetaCandidate &candidate : candidates)
            toEval.push_back(candidate.theta);
        std::vector<double> actual = computeM11Batch(toEval);

        size_t addedThisDepth = 0;
        for (size_t ci = 0; ci < candidates.size(); ++ci) {
            double tmid = candidates[ci].theta;
            double m_interp = candidates[ci].interp;
            double m_actual = actual[ci];
            double denom = std::max(fabs(m_actual), fabs(m_interp));
            double relErr = (denom > 1e-30) ? fabs(m_actual - m_interp) / denom : 0;
            if (relErr > eps) {
                thetaM11[tmid] = m_actual;
                ++addedThisDepth;
            }
        }
        totalAdded += (int)addedThisDepth;
        if (m_mpiRank == 0)
            std::cout << "  Depth " << depth << ": +" << addedThisDepth
                      << " points (total " << thetaM11.size() << ")" << std::endl;
        if (addedThisDepth == 0) break;
    }

    if (m_mpiRank == 0)
        std::cout << "Final adaptive grid: " << thetaM11.size()
                  << " theta points (theta eval "
                  << std::fixed << std::setprecision(2)
                  << theta_eval_seconds << " s"
                  << (thetaEvalUsedGpu ? ", GPU" : "")
                  << ")" << std::endl;

    // Step 4: Set the adaptive theta grid and re-run full TraceFromSobol
    // Free cached beams (we'll re-trace with new grid)
    allPrepared.clear();
    allPrepared.shrink_to_fit();

    std::vector<double> finalThetas;
    for (auto &kv : thetaM11)
        finalThetas.push_back(kv.first);

    // Set the new grid
    hp->m_sphere.thetaValues = finalThetas;
    hp->m_sphere.isNonUniform = true;
    hp->m_sphere.nZenith = finalThetas.size() - 1;
    hp->m_sphere.zenithStart = finalThetas.front();
    hp->m_sphere.zenithEnd = finalThetas.back();
    hp->m_sphere.zenithStep = (hp->m_sphere.zenithEnd - hp->m_sphere.zenithStart) / hp->m_sphere.nZenith;
    hp->m_sphere.ComputeSphereDirections(m_incidentLight);

    // Reallocate Mueller arrays
    int nAz = hp->m_sphere.nAzimuth;
    int nZen = hp->m_sphere.nZenith;
    hp->M = Arr2D(nAz + 1, nZen + 1, 4, 4);
    hp->M_noshadow = Arr2D(nAz + 1, nZen + 1, 4, 4);
    hp->SetScatteringSphere(hp->m_sphere);

    if (gridOnly) {
        if (m_mpiRank == 0)
            std::cout << "Adaptive theta grid set (" << finalThetas.size()
                      << " points). Full computation deferred." << std::endl;
        return;
    }

    if (m_mpiRank == 0)
        std::cout << "Full diffraction on " << finalThetas.size()
                  << " theta × " << (nAz+1) << " phi (re-tracing)..." << std::endl;

    // Re-run full computation with new grid
    m_incomingEnergy = 0;
    TraceFromSobol(nOrient, betaSym, gammaSym);
}

void TracerPOTotal::TraceAdaptive(double eps, double betaSym, double gammaSym, int maxOrientOverride)
{
    std::cout << "Adaptive mode: target relative change = " << eps << std::endl;
    std::cout << "  Convergence criteria: incoming energy, M11(22/46/90/180°)" << std::endl;
    std::cout << "  beta_sym=" << RadToDeg(betaSym) << " deg, gamma_sym="
              << RadToDeg(gammaSym) << " deg" << std::endl;
    std::vector<std::string> adaptiveLogLines;
    if (m_mpiRank == 0)
    {
        std::ostringstream log;
        log << "===== ADAPTIVE SOBOL =====\n";
        log << "target relative change = " << eps << "\n";
        log << "controls: incoming energy, M11(22), M11(46), M11(90), M11(180)\n";
        log << "beta_sym=" << RadToDeg(betaSym) << " deg, gamma_sym="
            << RadToDeg(gammaSym) << " deg\n";
        adaptiveLogLines.push_back(log.str());
    }

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
    // Δθ = 0.69 * λ / Dmax * (180/π), orient_step = Δθ / m_ringPoints
    // N_beta = betaSym / orient_step / 16, N_gamma = gammaSym / orient_step / 16
    // N_start = N_beta × N_gamma rounded down to power of 2
    double Dmax = m_particle->MaximalDimention();
    double delta_deg = 0.69 * m_scattering->m_wave / Dmax * (180.0 / M_PI);
    double orient_step = delta_deg / std::max(1, m_ringPoints);
    int nb16 = std::max(1, (int)(RadToDeg(betaSym) / orient_step / 16));
    int ng16 = std::max(1, (int)(RadToDeg(gammaSym) / orient_step / 16));
    int nStartRaw = nb16 * ng16;
    int nOrient = 1;
    while (nOrient * 2 <= nStartRaw) nOrient *= 2;
    if (nOrient < 64) nOrient = 64;
    std::cout << "  Start estimate (div16): " << nb16 << " x " << ng16
              << " = " << nStartRaw << " -> " << nOrient << " (power of 2)" << std::endl;
    if (m_mpiRank == 0)
    {
        std::ostringstream line;
        line << "Start estimate (div16): " << nb16 << " x " << ng16
             << " = " << nStartRaw << " -> " << nOrient << " (power of 2)\n";
        adaptiveLogLines.push_back(line.str());
    }

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
        // Physics cap: div2 of full grid (div1 is overkill, hours of compute)
        // User can override with --maxorient for div1 if needed
        int maxDiv2 = std::max(1024, maxP2 / 2);
        int p2 = 1; while (p2 * 2 <= maxDiv2) p2 *= 2;
        maxOrient = p2;
    }
    std::cout << "  Max estimate (div1): " << nb1 << " x " << ng1
              << " = " << maxFromPhysics << ", div2 cap -> " << maxOrient << std::endl;
    std::cerr << "Adaptive: max orientations = " << maxOrient
              << " (" << availMB_ad << " MB available)" << std::endl;
    if (m_mpiRank == 0)
    {
        std::ostringstream line;
        line << "Max estimate (div1): " << nb1 << " x " << ng1
             << " = " << maxFromPhysics << ", div2 cap -> " << maxOrient << "\n";
        line << "Adaptive: max orientations = " << maxOrient
             << " (" << availMB_ad << " MB available)\n";
        adaptiveLogLines.push_back(line.str());
    }

    auto t_start = std::chrono::high_resolution_clock::now();

    // Generate ALL Sobol points up to maxOrient at once (deterministic)
    Sobol2D sobol_gen(42);
    std::vector<double> su_all, sv_all;
    sobol_gen.generate(maxOrient, su_all, sv_all);

    // Accumulator: each batch is averaged internally and then combined by
    // its orientation count. This keeps incremental Sobol identical to a
    // single run over the same prefix of Sobol points.
    hp->M.ClearArr();
    hp->M_noshadow.ClearArr();
    hp->CleanJ();
    m_incomingEnergy = 0;
    hp->m_outputEnergy = 0;

    double prevM11_180 = 0, prevM11_90 = 0, prevM11_46 = 0, prevM11_22 = 0;
    double prevEnergy = 0;
    int totalOrient = 0;
    int convergedCount = 0;
    double ctrlWeightedSum[4] = {0,0,0,0}; // sum(batch_average * batch_size)
    double energyWeightedSum = 0.0;         // sum(batch_average * batch_size)

    for (int iter = 0; iter < 15; ++iter)
    {
        int batchStart = totalOrient;
        int batchEnd = nOrient;  // nOrient = target total after this iteration
        int batchSize = batchEnd - batchStart;

        if (batchSize <= 0) { nOrient *= 2; continue; }

        // Trace and diffract only the NEW orientations [batchStart..batchEnd).
        // Each batch is averaged with weight=1/batchSize, then accumulated
        // with weight=batchSize/totalOrient. This is required because after
        // the first doubling the added Sobol batches have different sizes.

        // MPI: each rank processes a subset of this batch
        int myBatchStart = m_mpiRank * batchSize / m_mpiSize;
        int myBatchEnd = (m_mpiRank + 1) * batchSize / m_mpiSize;
        int myBatchSize = myBatchEnd - myBatchStart;

        double batchEnergy = 0;
        double weight = 1.0 / batchSize;  // global weight (not per-rank)
        std::vector<Beam> outBeams;

        // Control angle indices (computed once)
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

        // Process in sub-chunks to limit memory (max 4096 per chunk)
        int subChunkMax = std::min(4096, myBatchSize);
        for (int sc = 0; sc < myBatchSize; sc += subChunkMax)
        {
            int scSize = std::min(subChunkMax, myBatchSize - sc);
            std::vector<PreparedOrientation> chunkPrep(scSize);

            for (int i = 0; i < scSize; ++i)
            {
                int idx = batchStart + myBatchStart + sc + i;
                double beta  = acos(1.0 - (1.0 - cosBetaSym) * su_all[idx]);
                double gamma = gammaSym * sv_all[idx];
                m_particle->Rotate(beta, gamma, 0);
                if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
                bool ok = m_scattering->ScatterLight(0, 0, outBeams);
                if (ok)
                    hp->PrepareBeams(outBeams, weight, chunkPrep[i]);
                else
                    chunkPrep[i].sinZenith = weight;
                batchEnergy += m_scattering->GetIncedentEnergy() * weight;
                outBeams.clear();
            }

            // Diffract at 4 control angles only
            for (int i = 0; i < scSize; ++i) {
                if (chunkPrep[i].beams.empty()) continue;
                double m11[4];
                hp->DiffractControlPoints(chunkPrep[i], ctrlIdx, 4, m11);
                for (int k = 0; k < 4; ++k) batchCtrl[k] += m11[k];
            }
        } // sub-chunks

        // MPI: reduce only THIS BATCH's control points + energy.
#ifdef USE_MPI
        if (m_mpiSize > 1)
        {
            double sbuf5[5] = {batchCtrl[0], batchCtrl[1], batchCtrl[2], batchCtrl[3], batchEnergy};
            double rbuf5[5];
            MPI_Allreduce(sbuf5, rbuf5, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            for (int k = 0; k < 4; ++k)
                batchCtrl[k] = rbuf5[k];
            batchEnergy = rbuf5[4];
        }
#endif

        for (int k = 0; k < 4; ++k)
            ctrlWeightedSum[k] += batchCtrl[k] * batchSize;
        energyWeightedSum += batchEnergy * batchSize;

        totalOrient = batchEnd;

        // Control points from fast DiffractControlPoints (not full grid).
        double energyAvg = (totalOrient > 0) ? energyWeightedSum / totalOrient : 0;
        double M11_22  = (totalOrient > 0) ? ctrlWeightedSum[0] / totalOrient : 0;
        double M11_46  = (totalOrient > 0) ? ctrlWeightedSum[1] / totalOrient : 0;
        double M11_90  = (totalOrient > 0) ? ctrlWeightedSum[2] / totalOrient : 0;
        double M11_180 = (totalOrient > 0) ? ctrlWeightedSum[3] / totalOrient : 0;

        auto relChange = [](double current, double previous) {
            double denom = std::max(std::fabs(previous), 1e-30);
            return std::fabs(current - previous) / denom;
        };

        double dEnergy = (iter > 0) ? relChange(energyAvg, prevEnergy) : 1.0;
        double dM22    = (iter > 0) ? relChange(M11_22, prevM11_22) : 1.0;
        double dM46    = (iter > 0) ? relChange(M11_46, prevM11_46) : 1.0;
        double dM90    = (iter > 0) ? relChange(M11_90, prevM11_90) : 1.0;
        double dM180   = (iter > 0) ? relChange(M11_180, prevM11_180) : 1.0;
        double dMax = std::max({dEnergy, dM22, dM46, dM90, dM180});

        std::cout << std::fixed << std::setprecision(2);
        if (m_mpiRank == 0)
        {
            std::ostringstream line;
            line << std::fixed << std::setprecision(2);
            line << "  N=" << totalOrient << " (+" << batchSize << ")"
                 << "  dE=" << dEnergy*100 << "%"
                 << "  d22=" << dM22*100 << "%"
                 << "  d46=" << dM46*100 << "%"
                 << "  d90=" << dM90*100 << "%"
                 << "  d180=" << dM180*100 << "%"
                 << "  max=" << dMax*100 << "%";
            std::cout << line.str() << std::endl;
            adaptiveLogLines.push_back(line.str() + "\n");
        }

        bool all_ok = (dMax < eps && iter > 0);

        if (all_ok)
            convergedCount++;
        else
            convergedCount = 0;

        if (convergedCount >= 2)
        {
            if (m_mpiRank == 0)
            {
                std::ostringstream line;
                line << "Converged at N=" << totalOrient
                     << " (all 5 controls within " << eps*100
                     << "% for 2 consecutive steps)";
                std::cout << line.str() << std::endl;
                adaptiveLogLines.push_back(line.str() + "\n");
            }
            break;
        }

        prevEnergy = energyAvg;
        prevM11_22 = M11_22;
        prevM11_46 = M11_46;
        prevM11_90 = M11_90;
        prevM11_180 = M11_180;
        nOrient *= 2;
        if (nOrient > maxOrient)
        {
            if (m_mpiRank == 0)
            {
                std::ostringstream line;
                line << "WARNING: Max orientations reached (N=" << totalOrient
                     << ", limit=" << maxOrient << "). Target accuracy "
                     << eps*100 << "% may not be achieved.";
                std::cout << line.str() << std::endl;
                std::cout << "  To improve: use --maxorient " << maxOrient*2 << std::endl;
                adaptiveLogLines.push_back(line.str() + "\n");
                adaptiveLogLines.push_back("  To improve: use --maxorient " + std::to_string(maxOrient*2) + "\n");
            }
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
        std::ostringstream log;
        log << "\n";
        for (const std::string &line : adaptiveLogLines)
            log << line;
        log << "Adaptive total time: " << std::fixed << std::setprecision(1)
            << total_sec << " s\n";
        AppendTextLog(log.str());
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
            for (int p = 0; p < nAz; ++p) m11_avg += hp->M(p, j)[0][0];
            m11_avg /= nAz;
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
            for (int p = 0; p < nAz2; ++p) m11_avg += hp->M(p, j)[0][0];
            m11_avg /= nAz2;
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
