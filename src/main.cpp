#include <iostream>
#include <fstream>
#include <assert.h>
#include <float.h>
#include <string>
#include <sys/stat.h>
#include <cstdlib>
#include <thread>
#include <set>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>

#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include "HexagonalAggregate.h"
#include "ConcaveHexagonal.h"
#include "CertainAggregate.h"
#include "Bullet.h"
#include "BulletRosette.h"
#include "global.h"
#include "TracerGO.h"
#include "TracerPOTotal.h"
#include "ArgPP.h"
#include "Tracks.h"
#include "HandlerPOTotal.h"
#include "HandlerTotalGO.h"
#include "HandlerTracksGO.h"
#include "Droxtal.h"
#include "GpuSupport.h"

#ifdef _OUTPUT_NRG_CONV
ofstream energyFile("energy.dat", ios::out);
double SSconfined=0;
int bcount=0;
#endif

using namespace std;
using namespace chrono;
using ::complex;

int DefaultPhysicalCoreCount()
{
#ifdef __linux__
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    std::set<std::string> cores;
    std::string physicalId = "0";
    std::string coreId;
    while (std::getline(cpuinfo, line))
    {
        if (line.compare(0, 12, "physical id") == 0)
        {
            size_t pos = line.find(':');
            if (pos != std::string::npos)
                physicalId = line.substr(pos + 1);
        }
        else if (line.compare(0, 7, "core id") == 0)
        {
            size_t pos = line.find(':');
            if (pos != std::string::npos)
            {
                coreId = line.substr(pos + 1);
                cores.insert(physicalId + ":" + coreId);
            }
        }
    }
    if (!cores.empty())
        return (int)cores.size();
#endif
    unsigned int n = std::thread::hardware_concurrency();
    return n > 1 ? (int)((n + 1) / 2) : 1;
}

void ConfigureDefaultOpenMP()
{
#ifdef _OPENMP
    if (std::getenv("OMP_NUM_THREADS") != nullptr)
        return;
    int cores = std::max(1, DefaultPhysicalCoreCount());
    omp_set_num_threads(cores);
    setenv("OMP_NUM_THREADS", std::to_string(cores).c_str(), 0);
#endif
}

void ConfigureOpenMPThreads(int threads)
{
#ifdef _OPENMP
    omp_set_num_threads(threads);
    setenv("OMP_NUM_THREADS", std::to_string(threads).c_str(), 1);
#else
    (void)threads;
    std::cerr << "WARNING: --threads ignored: OpenMP is not enabled in this build." << std::endl;
#endif
}

void ApplyBeamCutoffOptions(const ArgPP &args, Handler *handler)
{
    if (args.IsCatched("beam_cutoff"))
        handler->m_targetEps = args.GetDoubleValue("beam_cutoff", 0);
    if (args.IsCatched("beam_cutoff_j"))
        handler->m_beamCutoffJRel = args.GetDoubleValue("beam_cutoff_j", 0);
    if (args.IsCatched("beam_cutoff_area"))
        handler->m_beamCutoffAreaRel = args.GetDoubleValue("beam_cutoff_area", 0);
    if (args.IsCatched("beam_cutoff_importance"))
        handler->m_beamCutoffImportanceRel = args.GetDoubleValue("beam_cutoff_importance", 0);
}

void ApplyPOOutputOptions(const ArgPP &args, HandlerPO *handler, double wave)
{
    if (!handler)
        return;
    handler->SetFullOnly(!args.IsCatched("noshadow_output"));
    handler->m_otPhaseAverage = args.IsCatched("ot_phase_avg");
    if (args.IsCatched("ot_phase_shift"))
        handler->m_otFarReferencePath =
            args.GetDoubleValue("ot_phase_shift", 0) * wave;
}

void ApplyTraceCutoffOptions(const ArgPP &args, Scattering *scattering)
{
    if (args.IsCatched("trace_cutoff"))
    {
        double eps = args.GetDoubleValue("trace_cutoff", 0);
        scattering->m_traceCutoffJRel = eps;
        scattering->m_traceCutoffAreaRel = eps;
    }
    if (args.IsCatched("trace_cutoff_j"))
        scattering->m_traceCutoffJRel = args.GetDoubleValue("trace_cutoff_j", 0);
    if (args.IsCatched("trace_cutoff_area"))
        scattering->m_traceCutoffAreaRel = args.GetDoubleValue("trace_cutoff_area", 0);
    if (args.IsCatched("trace_cutoff_importance"))
        scattering->m_traceCutoffImportanceRel = args.GetDoubleValue("trace_cutoff_importance", 0);
    if (args.IsCatched("trace_max_beams"))
        scattering->m_traceMaxBeams = args.GetIntValue("trace_max_beams", 0);
    if (args.IsCatched("gpu_trace"))
        scattering->m_gpuTracePrefilter = true;
}

enum class ParticleType : int
{
    Hexagonal = 1,
    Bullet = 2,
    BulletRosette = 3,
    Droxtal = 4,
    ConcaveHexagonal = 10,
    TiltedHexagonal = 11,
    HexagonalAggregate = 12,
    CertainAggregate = 999
};

Tracks trackGroups;

void SetArgRules(ArgPP &parser)
{
    int zero = 0;
    parser.AddRule("p", '+', true); // particle (type, size, ...)
    parser.AddRule("ri", 2); // refractive index (Re and Im parts)
    parser.AddRule("n", 1, true); // number of internal reflection (optional for --autofull)
    parser.AddRule("pf", 1, true); // particle (filename)
    parser.AddRule("rs", 1, true, "pf"); // resize particle (new size)
    parser.AddRule("k_eq", 1, true); // resize particle to equivalent-volume size parameter
    parser.AddRule("fixed", 2, true); // fixed orientarion (beta, gamma)
    parser.AddRule("random", 2, true); // random orientarion (beta number, gamma number)
    parser.AddRule("montecarlo", 1, true); // random orientarion (beta number, gamma number)
    parser.AddRule("orientfile", 1, true); // orientations from file
    parser.AddRule("karczewski", 0, true); // use Karczewski polarization matrix
    parser.AddRule("go", 0, true); // geometrical optics method
    parser.AddRule("po", 0, true); // phisical optics method
    parser.AddRule("w", 1, true); // wavelength
    parser.AddRule("b", 2, true); // beta range (begin, end)
    parser.AddRule("g", 2, true); // gamma range (begin, end)
    parser.AddRule("grid", '+', true); /* backscattering grid:
 (radius, Nphi, Ntheta) when 3 parameters
 (theta1, theta2, Nphi, Ntheta) when 4 parameters*/
    parser.AddRule("point", zero, true, "po"); // calculate only backscatter point
    parser.AddRule("tr", 1, true); // file with trajectories
    parser.AddRule("all", 0, true); // calculate all trajectories
    parser.AddRule("abs", zero, true, "w"); // accounting of absorbtion
    parser.AddRule("abs_points", 1, true); // absorption samples: 1=center, all=all polygon vertices
    parser.AddRule("close", 0, true); // closing of program after calculation
    parser.AddRule("o", 1, true); // output folder name
    parser.AddRule("gr", zero, true);
    parser.AddRule("filter", 1, true); // scattering angle filter
    parser.AddRule("shadow", zero, true);
    parser.AddRule("incoh", zero, true);
    parser.AddRule("jones", zero, true);
    parser.AddRule("shadow_off", zero, true);
    parser.AddRule("full_only", zero, true);
    parser.AddRule("noshadow_output", zero, true);
    parser.AddRule("forced_nonconvex", zero, true);
    parser.AddRule("forced_convex", zero, true);
    parser.AddRule("r", 1, true); // restriction ratio for small beams when intersection (100 by default)
    parser.AddRule("log", 1, true); // time of writing progress (in seconds)
    parser.AddRule("multigrid", 3, true); // multi-size: Dmin Dmax Nsizes (log scale)
    parser.AddRule("multikeq", 3, true); // multi-size by k_eq: Kmin Kmax Nsizes (log scale)
    parser.AddRule("multikeq_list", 1, true); // multi-size by exact k_eq values from file
    parser.AddRule("multigrid_parallel", 1, true); // run multigrid sizes as child processes
    parser.AddRule("multigrid_threads", 1, true); // per-child OpenMP threads for multigrid_parallel
    parser.AddRule("save_betas", 0, true); // save intermediate Mueller for each beta to betas/ subfolder
    parser.AddRule("checkpoint", 0, true); // enable checkpoint save/resume for long orientfile runs
    parser.AddRule("tgrid", 1, true); // non-uniform theta grid file
    parser.AddRule("beam_cutoff", 1, true); // common relative beam cutoff
    parser.AddRule("beam_cutoff_j", 1, true); // relative |J|^2 beam cutoff
    parser.AddRule("beam_cutoff_area", 1, true); // relative area beam cutoff
    parser.AddRule("beam_cutoff_importance", 1, true); // relative |J|^2*area beam cutoff
    parser.AddRule("ot_phase_avg", 0, true); // average optical-theorem extinction over far-reference phase
    parser.AddRule("ot_phase_shift", 1, true); // diagnostic OT far-reference phase shift in wavelengths
    parser.AddRule("trace_cutoff", 1, true); // common relative tracing prune cutoff
    parser.AddRule("trace_cutoff_j", 1, true); // relative |J|^2 tracing prune cutoff
    parser.AddRule("trace_cutoff_area", 1, true); // relative area tracing prune cutoff
    parser.AddRule("trace_cutoff_importance", 1, true); // relative |J|^2*area tracing prune cutoff
    parser.AddRule("trace_max_beams", 1, true); // max traced beam nodes per orientation
    parser.AddRule("gpu_trace", 0, true); // experimental CUDA tracing prefilter
    parser.AddRule("sobol", 1, true); // Sobol quasi-random orientations (number, power of 2)
    parser.AddRule("sobol_seed", 2, true); // Sobol orientations with nested Owen scramble seed (N seed)
    parser.AddRule("sobol_ring", 2, true); // Sobol beta x uniform gamma ring (Nbeta Ngamma)
    parser.AddRule("hammersley", 1, true); // Hammersley orientation set (N)
    parser.AddRule("lattice", 1, true); // rank-1 lattice orientation set (N)
    parser.AddRule("lattice_z", 2, true); // rank-1 lattice with explicit generator (N Z)
    parser.AddRule("euler_quad", 2, true); // high-order Euler quadrature (Nbeta Ngamma)
    parser.AddRule("auto_tgrid", 1, true); // adaptive theta grid (arg: tolerance, e.g. 0.05)
    parser.AddRule("auto_phi", 0, true); // auto-select N_phi based on size parameter
    parser.AddRule("nphi", 1, true); // override N_phi (takes priority over --grid and --auto_phi)
    parser.AddRule("adaptive", 1, true); // adaptive convergence (target relative accuracy)
    parser.AddRule("autofull", 1, true); // full 3D sequential: n → N_phi → N_orient
    parser.AddRule("oldautofull", 1, true); // autofull search + oldauto regular final grid
    parser.AddRule("owen_avg", 1, true); // --autofull final: average K nested Owen Sobol seeds
    parser.AddRule("owen_seeds", '+', true); // explicit seeds for --autofull final averaging
    parser.AddRule("auto", 1, true); // full auto: auto_tgrid + auto_phi + adaptive (one arg: eps)
    parser.AddRule("maxorient", 1, true); // max orientations for adaptive (power of 2)
    parser.AddRule("chunk", 1, true); // max Sobol orientations per memory chunk
    parser.AddRule("oldauto", 1, true); // physics-based: div2/div4/div8 of diffraction-limited grid
    parser.AddRule("ring_points", 1, true); // points per diffraction ring for orientation estimates
    parser.AddRule("mirror_gamma", 0, true); // use mirror symmetry: gamma fundamental range is halved
    parser.AddRule("threads", 1, true); // OpenMP worker threads
    parser.AddRule("gpu", 0, true); // enable CUDA GPU backend
    parser.AddRule("cpu", 0, true); // force CPU backend in GPU-default builds
    parser.AddRule("fft", 0, true); // enable experimental FFT angular interpolation backend
    parser.AddRule("coh_orient", 0, true); // coherent across orientations (legacy mode)
    parser.AddRule("pole", 0, true); // fast pole shortcut: one gamma value at beta poles
    parser.AddRule("legacy_sign", 0, true); // use old (+) Fresnel sign for forward direction
    parser.AddRule("sym", 2, true); // symmetry override: beta_factor gamma_factor (e.g. --sym 2 6)
    parser.AddRule("help", 0, true); // print help
}

void PrintHelp()
{
    using namespace std;
    cout << "MBS-fast: Physical Optics for Ice Crystals\n"
         << "Usage: mbs_po --po [orientation] [options] -p TYPE L D [-w LAMBDA] [--ri Re Im] [-n N]\n\n"

         << "=== Particle ===\n"
         << "  -p TYPE L D [extra]    Particle: 1=hex, 2=bullet, 3=rosette, 4=droxtal,\n"
         << "                         10=concave hex, 12=hex aggregate\n"
         << "  --pf FILE              Particle from .obj file\n"
         << "  --rs SIZE              Resize particle to Dmax=SIZE (with --pf)\n"
         << "  --k_eq X               Resize particle so 2*pi*r_eq/lambda = X\n"
         << "  --ri Re Im             Refractive index (default 1.31 0)\n"
         << "  -w LAMBDA              Wavelength in um (default 0.532)\n"
         << "  -n N                   Max internal reflections (default 6)\n\n"

         << "=== Method ===\n"
         << "  --po                   Physical optics (default)\n"
         << "  --go                   Geometric optics\n\n"

         << "=== Orientation ===\n"
         << "  --fixed BETA GAMMA     Single orientation (degrees)\n"
         << "  --random Nb Ng         Regular beta x gamma grid\n"
         << "  --sobol N              Sobol quasi-random (N orientations)\n"
         << "  --sobol_seed N S       Sobol with nested Owen scramble seed S\n"
         << "  --sobol_ring Nb Ng     Sobol beta x shifted uniform gamma ring\n"
         << "  --hammersley N         Hammersley low-discrepancy orientations\n"
         << "  --lattice N            Rank-1 lattice orientations\n"
         << "  --lattice_z N Z        Rank-1 lattice with explicit generator Z\n"
         << "  --euler_quad Nb Ng     High-order: Gauss in cos(beta) x periodic gamma\n"
         << "  --montecarlo N         Monte Carlo random (N orientations)\n"
         << "  --adaptive EPS         Adaptive Sobol (converge to EPS relative accuracy)\n"
         << "  --auto EPS             Full auto: adaptive theta + phi + orientations\n"
         << "  --autofull EPS         Full auto including n search\n"
         << "  --oldautofull EPS      Autofull search, oldauto regular final orientation grid\n"
         << "  --owen_avg K           With --autofull, average K nested Owen final seeds (default 5)\n"
         << "  --owen_seeds S...      Explicit final Owen seeds for --autofull averaging\n"
         << "  --oldauto DIV          Physics-based grid (div2/div4/div8 of diffraction limit)\n"
         << "  --ring_points N        Points per diffraction ring for orientation estimates (default 3)\n"
         << "  --mirror_gamma         Use mirror symmetry: gamma range is halved (60 -> 30 deg)\n"
         << "  --orientfile FILE      Orientations from file (beta gamma per line)\n"
         << "  --b B1 B2              Beta range in degrees for --random\n"
         << "  --g G1 G2              Gamma range in degrees for --random\n"
         << "  --maxorient N          Max orientations for adaptive (power of 2)\n"
         << "  --chunk N              Max orientations/gamma values per memory chunk\n"
         << "  --coh_orient           Coherent across orientations (legacy)\n"
         << "  --pole                 Fast pole gamma: use one gamma value at beta poles\n\n"

         << "=== Scattering grid ===\n"
         << "  --grid T1 T2 Nphi Nth  Theta range [T1,T2] deg, Nphi azimuth, Nth zenith\n"
         << "  --grid R Nphi Nth      Backscatter cone with angular radius R deg\n"
         << "  --tgrid FILE           Non-uniform theta grid from file (degrees, one per line)\n"
         << "  --auto_tgrid EPS       Adaptive theta grid via bisection (tolerance EPS)\n"
         << "  --auto_phi             Auto N_phi = x/5 + 48\n"
         << "  --nphi N               Override N_phi (highest priority)\n\n"

         << "=== Grid priority ===\n"
         << "  Theta: --tgrid > --grid > --auto_tgrid > --auto > default\n"
         << "  Phi:   --nphi  > --grid > --auto_phi   > --auto > default\n\n"

         << "=== Optimization ===\n"
         << "  --threads N            OpenMP worker threads (default: physical cores)\n"
         << "  --gpu                  Use CUDA GPU backend for diffraction (default in gpu/ build)\n"
         << "  --cpu                  Force CPU backend in gpu/ build\n"
         << "  --fft                  With --gpu, use experimental cuFFT phi interpolation backend\n"
         << "                         Auto modes select phi factor from EPS/Nphi; env MBS_FFT_PHI_FACTOR=N overrides\n"
         << "  --beam_cutoff EPS      Set both relative beam cutoffs below; skip if either matches\n"
         << "  --beam_cutoff_j EPS    Skip beams with |J|^2/max < EPS (0 disables J test)\n"
         << "  --beam_cutoff_area EPS Skip beams with area/max < EPS (0 disables area test)\n"
         << "  --beam_cutoff_importance EPS Skip beams with |J|^2*area/max < EPS\n"
         << "  --ot_phase_avg         Average OT extinction over one far-reference phase period\n"
         << "  --ot_phase_shift F     Diagnostic OT phase shift in wavelengths (default 0)\n"
         << "  --trace_cutoff EPS     Stop tracing internal beams if either relative test matches\n"
         << "  --trace_cutoff_j EPS   Trace prune by |J|^2/initial-max < EPS\n"
         << "  --trace_cutoff_area EPS Trace prune by area/initial-max < EPS\n"
         << "  --trace_cutoff_importance EPS Trace prune by |J|^2*area/max < EPS\n"
         << "  --trace_max_beams N    Abort one orientation after N traced beam nodes (0 disables)\n"
         << "  --gpu_trace            Experimental CUDA prefilter for nonconvex tracing candidates\n"
         << "                         Env: MBS_GPU_TRACE_BATCH_BEAMS=1024, MBS_GPU_TRACE_MIN_CANDIDATES=8192\n"
         << "  -r RATIO               Beam area restriction ratio (default 100)\n"
         << "  --sym Sb Sg            Override symmetry: beta/Sb, 360/Sg degrees\n"
         << "  --filter DEG           Restrict output to backscattering cone\n"
         << "  --point                Backscatter point mode (legacy; not for optimized --random)\n"
         << "  --shadow               Legacy flag, currently no effect\n"
         << "  --shadow_off           Disable shadow beam\n"
         << "  --full_only            Do not compute/write no-shadow Mueller output (default)\n"
         << "  --noshadow_output      Also compute/write _noshadow Mueller output\n"
         << "  --forced_convex        Force convex-particle processing\n"
         << "  --forced_nonconvex     Force non-convex-particle processing\n\n"

         << "=== Trajectories ===\n"
         << "  --tr FILE              Load trajectory file\n"
         << "  --all                  Calculate all loaded trajectories\n"
         << "  --gr                   Output trajectory groups\n\n"

         << "=== Multi-size ===\n"
         << "  --multigrid Dmin Dmax N  N sizes from Dmin to Dmax (log scale)\n"
         << "                           -p must specify largest particle\n\n"
         << "  --multikeq Kmin Kmax N   N equivalent-size parameters k_eq from Kmin to Kmax (log scale)\n"
         << "  --multikeq_list FILE     Exact k_eq values, one per line; --oldauto traces max once\n"
         << "  --multigrid_parallel N   Run multigrid/multikeq as N child processes\n"
         << "                           Useful with --gpu; default child threads is 1 unless overridden\n"
         << "  --multigrid_threads N    OpenMP threads per child in --multigrid_parallel\n\n"


         << "=== Output ===\n"
         << "  -o NAME                Output path/name\n"
         << "  --close                Exit after computation\n"
         << "  --save_betas           Save per-beta Mueller to _betas/ folder\n"
         << "  --checkpoint           Enable checkpoint save/resume for --orientfile\n"
         << "  --incoh                Incoherent per-beam Mueller (no Jones sum)\n"
         << "  --jones                Output Jones matrices\n"
         << "  --abs                  Enable absorption (requires Im(ri) > 0)\n"
         << "  --abs_points N|all     Absorption samples: 1=center (default), all=all polygon vertices\n"
         << "  --karczewski           Use Karczewski polarization matrix\n"
         << "  --legacy_sign          Use old (+) Fresnel sign\n"
         << "  --log SEC              Progress output interval (seconds)\n\n"
         << "  --help, -h             Print this help\n\n"

         << "=== Examples ===\n"
         << "  # Full auto, hex column 100x70 um\n"
         << "  mbs_po --po --auto 0.05 -p 1 100 70 -w 0.532 --ri 1.31 0 -n 8 --close\n\n"
         << "  # Manual Sobol with adaptive theta\n"
         << "  mbs_po --po --sobol 1024 --auto_tgrid 0.05 --auto_phi -p 1 100 70 ...\n\n"
         << "  # Multi-size scan\n"
         << "  mbs_po --po --sobol 1024 --multigrid 50 500 20 -p 1 500 62.5 ...\n\n"
         << "  # Physics-based grid\n"
         << "  mbs_po --po --oldauto 8 -p 1 200 25 -w 0.532 --ri 1.31 0 -n 4 ...\n"
         << endl;
}

static int RoundUpToMultiple(int value, int multiple)
{
    if (value <= 0 || multiple <= 1)
        return value;
    return ((value + multiple - 1) / multiple) * multiple;
}

static int MirrorSafePhiCount(ArgPP &parser, int nphi)
{
    if (nphi <= 0)
        return nphi;
    if (parser.IsCatched("mirror_gamma") || parser.IsCatched("fft"))
    {
        int rounded = RoundUpToMultiple(nphi, 6);
        if (rounded != nphi)
        {
            std::cerr << "WARNING: "
                      << (parser.IsCatched("mirror_gamma") ? "--mirror_gamma" : "--fft")
                      << " requires N_phi divisible by 6 for exact phi-sector "
                      << "mapping; using N_phi=" << rounded
                      << " instead of " << nphi << std::endl;
        }
        return rounded;
    }
    if (nphi % 6 != 0)
    {
        std::cerr << "WARNING: N_phi=" << nphi
                  << " is not divisible by 6. This is allowed for direct "
                  << "calculation, but exact hexagonal phi-sector reuse "
                  << "requires N_phi % 6 == 0." << std::endl;
    }
    return nphi;
}

static void SetRangeNphi(ArgPP &parser, ScatteringRange &range, int nphi)
{
    nphi = MirrorSafePhiCount(parser, nphi);
    if (nphi > 0)
    {
        range.nAzimuth = nphi;
        range.azinuthStep = 2.0 * M_PI / nphi;
    }
}

static void MakeMirrorGammaRangeConsistent(ArgPP &parser, Particle *particle,
                                           AngleRange &gamma)
{
    if (!parser.IsCatched("mirror_gamma"))
        return;

    if (!parser.IsCatched("g"))
    {
        gamma.min = 0.0;
        gamma.max = 0.5 * particle->GetSymmetry().gamma;
        gamma.norm = gamma.max - gamma.min;
        gamma.step = gamma.norm / gamma.number;
        std::cout << "Gamma mirror symmetry: random gamma range reduced to 0.."
                  << RadToDeg(gamma.max) << " deg" << std::endl;
    }

    if (gamma.number % 2 != 0)
    {
        int old = gamma.number;
        gamma.number += 1;
        gamma.step = gamma.norm / gamma.number;
        std::cerr << "WARNING: --mirror_gamma works best with an even number "
                  << "of gamma points in the half-domain; using "
                  << gamma.number << " instead of " << old << std::endl;
    }
}

/// Apply --nphi override if present (highest priority for N_phi)
void ApplyNphiOverride(ArgPP &parser, ScatteringRange &range)
{
    if (parser.IsCatched("nphi"))
    {
        int nphi = parser.GetIntValue("nphi", 0);
        SetRangeNphi(parser, range, nphi);
    }
}

void ApplyAbsorptionPointOption(ArgPP &parser, Handler *handler)
{
    if (!parser.IsCatched("abs_points"))
        return;

    std::string value = parser.GetStringValue("abs_points", 0);
    int n = 1;

    if (value == "all")
    {
        n = -1;
    }
    else
    {
        try
        {
            n = std::stoi(value);
        }
        catch (...)
        {
            n = 0;
        }
    }

    if (n == 0 || n < -1)
    {
        std::cerr << "WARNING: --abs_points must be positive or 'all'; using 1 (center)." << std::endl;
        n = 1;
    }

    handler->SetAbsorptionPointCount(n);
}

ScatteringRange SetConus(ArgPP &parser)
{
    ScatteringRange range(0, M_PI, 1, 1); // placeholder

    if (parser.GetArgNumber("grid") == 3)
    {
        double radius = parser.GetDoubleValue("grid", 0);
        int nAz = MirrorSafePhiCount(parser, parser.GetDoubleValue("grid", 1));
        int nZen = parser.GetDoubleValue("grid", 2);
        range = ScatteringRange(M_PI - DegToRad(radius), M_PI, nAz, nZen);
    }
    else if (parser.GetArgNumber("grid") == 4)
    {
        double zenStart = parser.GetDoubleValue("grid", 0);
        double zenEnd = parser.GetDoubleValue("grid", 1);

        if (zenStart < zenEnd)
        {
            int nAz = MirrorSafePhiCount(parser, parser.GetDoubleValue("grid", 2));
            int nZen = parser.GetDoubleValue("grid", 3);
            range = ScatteringRange(DegToRad(zenStart), DegToRad(zenEnd), nAz, nZen);
        }
        else
        {
            std::cerr << "ERROR!!!!!!! In \"grid\" arg 1 must be less than arg 2";
            exit(1);
        }
    }
    else if (parser.IsCatched("grid"))
    {
        std::cerr << "ERROR!!!!!!! Wrong \"grid\" argument number";
        exit(1);
    }
    // else: no --grid → use placeholder (0, pi, 1, 1), will be overwritten by --auto_tgrid or --tgrid

    // Apply non-uniform theta grid if --tgrid is specified
    if (parser.IsCatched("tgrid"))
    {
        std::string tgridFile = parser.GetStringValue("tgrid", 0);
        if (!range.LoadThetaGrid(tgridFile))
        {
            std::cerr << "ERROR: Cannot load theta grid from file: " << tgridFile << std::endl;
            exit(1);
        }
        // If --tgrid without --grid, set default N_phi=48
        if (!parser.IsCatched("grid"))
        {
            SetRangeNphi(parser, range, 48);
        }
        std::cout << "Non-uniform theta grid: " << (range.nZenith + 1)
                  << " points from " << RadToDeg(range.zenithStart)
                  << " to " << RadToDeg(range.zenithEnd) << " deg" << std::endl;
    }

    return range;
}

// (ParseMultigridFile removed — multigrid now takes Dmin Dmax N directly)

/// Generate auto theta grid based on size parameter x = pi * D / lambda
void ApplyAutoThetaGrid(ScatteringRange &range, double D, double wave)
{
    if (wave <= 0 || D <= 0) return;

    double x = M_PI * D / wave;
    double peak_width_deg = 180.0 / x; // diffraction peak half-width in degrees

    std::vector<double> thetas; // in degrees

    // Adaptive grid with 4 zones:
    // 1. Fine: resolve forward diffraction peak (Δθ = peak_width/10)
    // 2. Transition: gradual coarsening (geometric spacing)
    // 3. Medium: side scattering (Δθ = 1°)
    // 4. Coarse: backscattering where Mueller oscillates slowly (Δθ = 2°)

    // Zone 1: Forward peak [0, fine_end]
    // Fine step: resolve peak with at least 20 points across half-width
    double fine_step = std::max(0.002, std::min(peak_width_deg / 20.0, 0.5));
    // Extend fine zone to capture side lobes (important for Q_sca integration)
    double fine_end = std::min(10.0 * peak_width_deg, 20.0);
    // For small x (<30): peak is wide, fine zone covers more
    if (x < 30) fine_end = std::min(15.0 * peak_width_deg, 40.0);

    for (double t = 0; t <= fine_end + 1e-9; t += fine_step)
        thetas.push_back(t);

    // Zone 2: Transition [fine_end, transition_end] — geometric spacing
    double transition_end = std::max(fine_end + 1.0, 20.0);
    if (x < 30) transition_end = std::max(fine_end + 1.0, 40.0);

    if (fine_end < transition_end)
    {
        // Geometric spacing with max step ≤ 1°
        double t = fine_end;
        double step = fine_step;
        while (t < transition_end - 1e-9)
        {
            step = std::min(step * 1.3, 1.0); // grow by 30% per step, max 1°
            t += step;
            if (t > transition_end) t = transition_end;
            thetas.push_back(t);
        }
    }

    // Zone 3: Medium [transition_end, 120°] — 0.5° for small x, 1° otherwise
    // BUT: add fine zones around halo angles (22° and 46° for hex prisms)
    double medium_step = (x < 100) ? 0.5 : 1.0;

    // Halo fine zones: 0.1° step around 22° and 46° (±4°)
    // Only for particles large enough to produce halos (x > 50)
    double halo_step = std::max(0.05, std::min(0.2, peak_width_deg));
    double halo1_center = 22.0;  // 22° halo (two prism faces)
    double halo2_center = 46.0;  // 46° halo (basal + prism face)
    double halo_half = 4.0;      // ±4° around halo center

    for (double t = transition_end + medium_step; t <= 120.0 + 1e-9; )
    {
        thetas.push_back(t);

        // Check if next step enters a halo zone
        bool in_halo1 = (x > 50) && (t >= halo1_center - halo_half) && (t <= halo1_center + halo_half);
        bool in_halo2 = (x > 50) && (t >= halo2_center - halo_half) && (t <= halo2_center + halo_half);

        if (in_halo1 || in_halo2)
            t += halo_step;
        else
            t += medium_step;
    }

    // Zone 4: Side/back [120°, 175°] — 1° step
    for (double t = 121.0; t <= 175.0 + 1e-9; t += 1.0)
        thetas.push_back(t);

    // Zone 5: Near-backscattering [175°, 180°] — 0.25° step (LDR, depolarization)
    for (double t = 175.25; t <= 180.0 + 1e-9; t += 0.25)
        thetas.push_back(t);

    // Remove duplicates and sort
    std::sort(thetas.begin(), thetas.end());
    std::vector<double> unique_thetas;
    unique_thetas.push_back(thetas[0]);
    for (size_t i = 1; i < thetas.size(); ++i)
    {
        if (thetas[i] - unique_thetas.back() > 0.001)
            unique_thetas.push_back(thetas[i]);
    }

    // Convert to radians and load into range
    range.thetaValues.clear();
    for (double t : unique_thetas)
        range.thetaValues.push_back(DegToRad(t));

    range.isNonUniform = true;
    range.nZenith = (int)range.thetaValues.size() - 1;
    range.zenithStart = range.thetaValues.front();
    range.zenithEnd = range.thetaValues.back();
    range.zenithStep = (range.zenithEnd - range.zenithStart) / range.nZenith;

    std::cout << "Auto theta grid: x=" << x << ", peak_width=" << peak_width_deg
              << " deg, " << range.thetaValues.size() << " points" << std::endl;
}

std::vector<double> GenerateLogSizes(double minValue, double maxValue, int count)
{
    if (count < 1)
        count = 1;
    std::vector<double> values;
    values.reserve(count);
    if (count == 1)
    {
        values.push_back(maxValue);
        return values;
    }
    double logMin = std::log(minValue);
    double logMax = std::log(maxValue);
    for (int i = 0; i < count; ++i)
    {
        double v = std::exp(logMin + (logMax - logMin) * i / (count - 1));
        values.push_back(v);
    }
    values.back() = maxValue;
    return values;
}

std::vector<unsigned int> BuildOwenAverageSeeds(const ArgPP &args)
{
    static const unsigned int defaults[] = {
        7u, 123u, 777u, 2026u, 9001u,
        42u, 314159u, 271828u, 1618033u, 112358u
    };

    std::vector<unsigned int> seeds;
    if (args.IsCatched("owen_seeds"))
    {
        unsigned n = args.GetArgNumber("owen_seeds");
        for (unsigned i = 0; i < n; ++i)
            seeds.push_back((unsigned int)std::max(0, args.GetIntValue("owen_seeds", i)));
        return seeds;
    }

    if (!args.IsCatched("owen_avg"))
    {
        int defaultCount = 5;
        if (const char *env = std::getenv("MBS_AUTOFULL_DEFAULT_OWEN_AVG"))
        {
            char *end = nullptr;
            long parsed = std::strtol(env, &end, 10);
            if (end && *end == '\0')
                defaultCount = (int)parsed;
        }
        if (defaultCount <= 0)
            return seeds;
        int count = std::max(1, defaultCount);
        seeds.reserve(count);
        for (int i = 0; i < count; ++i)
        {
            if (i < (int)(sizeof(defaults) / sizeof(defaults[0])))
                seeds.push_back(defaults[i]);
            else
                seeds.push_back((unsigned int)(2654435761u * (unsigned int)(i + 1) + 1013904223u));
        }
        return seeds;
    }

    int count = std::max(1, args.GetIntValue("owen_avg", 0));
    seeds.reserve(count);
    for (int i = 0; i < count; ++i)
    {
        if (i < (int)(sizeof(defaults) / sizeof(defaults[0])))
            seeds.push_back(defaults[i]);
        else
            seeds.push_back((unsigned int)(2654435761u * (unsigned int)(i + 1) + 1013904223u));
    }
    return seeds;
}

std::vector<double> ReadSizeList(const std::string &fileName)
{
    std::ifstream in(fileName.c_str());
    if (!in)
        throw std::runtime_error("cannot open size list: " + fileName);

    std::vector<double> values;
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line))
    {
        ++lineNo;
        std::istringstream ss(line);
        std::string first;
        if (!(ss >> first))
            continue;
        if (!first.empty() && first[0] == '#')
            continue;

        char *end = nullptr;
        errno = 0;
        double value = std::strtod(first.c_str(), &end);
        if (errno != 0 || end == first.c_str() || (end && *end != '\0') || value <= 0.0)
        {
            std::ostringstream msg;
            msg << "bad positive size in " << fileName << ":" << lineNo;
            throw std::runtime_error(msg.str());
        }
        values.push_back(value);
    }

    if (values.empty())
        throw std::runtime_error("empty size list: " + fileName);
    return values;
}

std::string SizeLabel(double value)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6) << value;
    std::string s = ss.str();
    while (!s.empty() && s.back() == '0')
        s.pop_back();
    if (!s.empty() && s.back() == '.')
        s.pop_back();
    for (char &c : s)
    {
        if (c == '.')
            c = 'p';
        else if (c == '-')
            c = 'm';
    }
    return s.empty() ? "0" : s;
}

std::string ShellQuote(const std::string &s)
{
    std::string out = "'";
    for (char c : s)
    {
        if (c == '\'')
            out += "'\\''";
        else
            out += c;
    }
    out += "'";
    return out;
}

bool FlagMatches(const std::string &arg, const std::string &name)
{
    return arg == "--" + name || (name.size() == 1 && arg == "-" + name);
}

int ParallelMultigridRemoveCount(const std::string &arg)
{
    struct Item { const char *name; int values; };
    static const Item items[] = {
        {"multigrid", 3}, {"multikeq", 3}, {"multikeq_list", 1}, {"multigrid_parallel", 1},
        {"multigrid_threads", 1}, {"rs", 1}, {"k_eq", 1},
        {"threads", 1}, {"o", 1}
    };
    for (const Item &item : items)
    {
        if (FlagMatches(arg, item.name))
            return item.values;
    }
    return -1;
}

int RunParallelMultigrid(int argc, const char *argv[], const ArgPP &args, bool useGpu)
{
    bool byD = args.IsCatched("multigrid");
    bool byKEq = args.IsCatched("multikeq");
    if (byD == byKEq)
    {
        std::cerr << "ERROR: --multigrid_parallel requires exactly one of --multigrid or --multikeq." << std::endl;
        return 1;
    }
    if (byD && !args.IsCatched("pf"))
    {
        std::cerr << "ERROR: parallel --multigrid currently requires --pf, because it resizes file particles with --rs." << std::endl;
        return 1;
    }

    int jobs = args.GetIntValue("multigrid_parallel", 0);
    if (jobs < 1)
    {
        std::cerr << "ERROR: --multigrid_parallel must be >= 1." << std::endl;
        return 1;
    }
    int childThreads = args.IsCatched("multigrid_threads")
        ? args.GetIntValue("multigrid_threads", 0)
        : (useGpu ? 1 : std::max(1, DefaultPhysicalCoreCount() / jobs));
    if (childThreads < 1)
    {
        std::cerr << "ERROR: --multigrid_threads must be >= 1." << std::endl;
        return 1;
    }

    std::string key = byKEq ? "multikeq" : "multigrid";
    double minValue = args.GetDoubleValue(key, 0);
    double maxValue = args.GetDoubleValue(key, 1);
    int count = args.GetIntValue(key, 2);
    if (minValue <= 0 || maxValue <= 0)
    {
        std::cerr << "ERROR: multigrid sizes must be positive." << std::endl;
        return 1;
    }
    std::vector<double> sizes = GenerateLogSizes(minValue, maxValue, count);

    std::string baseOut = args.IsCatched("o") ? args.GetStringValue("o", 0) : "multigrid_parallel";
    mkdir(baseOut.c_str(), 0755);

    std::vector<std::string> baseArgs;
    baseArgs.reserve(argc + 8);
    baseArgs.push_back(argv[0]);
    for (int i = 1; i < argc; ++i)
    {
        std::string a(argv[i]);
        int skip = ParallelMultigridRemoveCount(a);
        if (skip >= 0)
        {
            i += skip;
            continue;
        }
        baseArgs.push_back(a);
    }

    std::cout << "Parallel multigrid: " << sizes.size() << " "
              << (byKEq ? "k_eq" : "Dmax") << " sizes, jobs=" << jobs
              << ", child_threads=" << childThreads
              << (useGpu ? " (--gpu)" : "") << std::endl;
    std::cout << "Output root: " << baseOut << std::endl;

    struct RunningChild { pid_t pid; std::string label; };
    std::vector<RunningChild> running;
    int failures = 0;

    auto waitOne = [&]() {
        int status = 0;
        pid_t pid = wait(&status);
        if (pid <= 0)
            return;
        auto it = std::find_if(running.begin(), running.end(),
            [pid](const RunningChild &c) { return c.pid == pid; });
        std::string label = (it == running.end()) ? std::to_string(pid) : it->label;
        if (it != running.end())
            running.erase(it);
        if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
        {
            ++failures;
            std::cerr << "Parallel multigrid child failed: " << label
                      << " status=" << status << std::endl;
        }
        else
        {
            std::cout << "Parallel multigrid child done: " << label << std::endl;
        }
    };

    for (double value : sizes)
    {
        while ((int)running.size() >= jobs)
            waitOne();

        std::string label = std::string(byKEq ? "keq" : "D") + SizeLabel(value);
        std::vector<std::string> child = baseArgs;
        child.push_back("--threads");
        child.push_back(std::to_string(childThreads));
        if (byKEq)
        {
            child.push_back("--k_eq");
            child.push_back(std::to_string(value));
        }
        else
        {
            child.push_back("--rs");
            child.push_back(std::to_string(value));
        }
        child.push_back("-o");
        child.push_back(baseOut + "/" + label);

        std::cout << "Starting " << label << ":";
        for (const std::string &s : child)
            std::cout << " " << ShellQuote(s);
        std::cout << std::endl;

        pid_t pid = fork();
        if (pid < 0)
        {
            std::cerr << "ERROR: fork failed: " << std::strerror(errno) << std::endl;
            ++failures;
            continue;
        }
        if (pid == 0)
        {
            std::string logPath = baseOut + "/" + label + ".run.log";
            int fd = open(logPath.c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0644);
            if (fd >= 0)
            {
                dup2(fd, STDOUT_FILENO);
                dup2(fd, STDERR_FILENO);
                close(fd);
            }
            std::vector<char*> execArgs;
            execArgs.reserve(child.size() + 1);
            for (std::string &s : child)
                execArgs.push_back(const_cast<char*>(s.c_str()));
            execArgs.push_back(nullptr);
            execvp(execArgs[0], execArgs.data());
            std::cerr << "ERROR: exec failed: " << std::strerror(errno) << std::endl;
            _exit(127);
        }
        running.push_back({pid, label});
    }

    while (!running.empty())
        waitOne();

    if (failures)
        std::cerr << "Parallel multigrid finished with " << failures << " failed child process(es)." << std::endl;
    else
        std::cout << "Parallel multigrid finished successfully." << std::endl;
    return failures ? 1 : 0;
}

AngleRange GetRange(const ArgPP &parser, const std::string &key,
                    Particle *particle)
{
    int number;
    double min, max;

    if (key == "b")
    {
        number = parser.GetIntValue("random", 0);

        if (parser.IsCatched("b"))
        {
            min = DegToRad(parser.GetDoubleValue(key, 0));
            max = DegToRad(parser.GetDoubleValue(key, 1));
        }
        else
        {
            min = 0;
            max = particle->GetSymmetry().beta;
        }
    }
    else if (key == "g")
    {
        number = parser.GetIntValue("random", 1);

        if (parser.IsCatched("g"))
        {
            min = DegToRad(parser.GetDoubleValue(key, 0));
            max = DegToRad(parser.GetDoubleValue(key, 1));
        }
        else
        {
            min = 0;
            max = particle->GetSymmetry().gamma;
        }
    }
    else
    {
        cerr << "Error! " << __FUNCTION__;
        throw std::exception();
    }

    return AngleRange(min, max, number);
}

int main(int argc, const char* argv[])
{
    ConfigureDefaultOpenMP();

    int mpi_rank = 0, mpi_size = 1;
#ifdef USE_MPI
    MPI_Init(&argc, const_cast<char***>(&argv));
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_rank == 0 && mpi_size > 1)
        std::cout << "MPI: " << mpi_size << " processes" << std::endl;
#endif

    RenameConsole("MBS");

    std::string additionalSummary;

    // Build command line string for logging
    std::string cmdLine;
    for (int i = 0; i < argc; ++i) {
        if (i > 0) cmdLine += " ";
        cmdLine += argv[i];
    }
    additionalSummary = "Command: " + cmdLine + "\n";

    if (argc <= 1)
    {
        PrintHelp();
        return 0;
    }

    // Check --help before full parse (avoids required arg errors)
    bool rawGpu = false;
    bool rawCpu = false;
    for (int i = 1; i < argc; ++i) {
        std::string a(argv[i]);
        if (a == "--help" || a == "-h") { PrintHelp(); return 0; }
        if (a == "--gpu") rawGpu = true;
        if (a == "--cpu") rawCpu = true;
    }
    if (rawGpu && rawCpu)
    {
        std::cerr << "ERROR: --gpu and --cpu are mutually exclusive." << std::endl;
        return 1;
    }

    ArgPP args;
    SetArgRules(args);
    args.Parse(argc, argv);

    if (args.IsCatched("threads"))
    {
        int threads = args.GetIntValue("threads", 0);
        if (threads < 1)
        {
            std::cerr << "ERROR: --threads must be >= 1." << std::endl;
            return 1;
        }
        ConfigureOpenMPThreads(threads);
    }

    if (args.IsCatched("gpu") && args.IsCatched("cpu"))
    {
        std::cerr << "ERROR: --gpu and --cpu are mutually exclusive." << std::endl;
        return 1;
    }

#if defined(USE_CUDA) && defined(MBS_GPU_DEFAULT_ON)
    const bool defaultGpu = true;
#else
    const bool defaultGpu = false;
#endif
    const bool useGpu = args.IsCatched("gpu") || (defaultGpu && !args.IsCatched("cpu"));

    if (useGpu)
    {
        GpuDeviceInfo gpuInfo;
        std::string gpuError;
        if (!CheckGpuRuntime(gpuInfo, gpuError))
        {
            std::cerr << "ERROR: CUDA GPU backend is unavailable: "
                      << gpuError << std::endl;
            return 1;
        }

        std::cout << "GPU backend: " << FormatGpuInfo(gpuInfo)
                  << (defaultGpu && !args.IsCatched("gpu") ? " (default)" : "")
                  << std::endl;
    }

    const bool useFft = args.IsCatched("fft");
    if (useFft && !useGpu)
    {
        std::cerr << "ERROR: --fft currently requires --gpu." << std::endl;
        return 1;
    }

    if (args.IsCatched("p") == args.IsCatched("pf"))
    {
        std::cerr << "ERROR: specify exactly one particle source: -p ... or --pf FILE." << std::endl;
        return 1;
    }
    if (args.IsCatched("rs") && args.IsCatched("k_eq"))
    {
        std::cerr << "ERROR: --rs and --k_eq are both size controls; use only one." << std::endl;
        return 1;
    }
    if (args.IsCatched("multigrid_parallel"))
    {
        return RunParallelMultigrid(argc, argv, args, useGpu);
    }

    double re = args.GetDoubleValue("ri", 0);
    double im = args.GetDoubleValue("ri", 1);
    double wave = args.IsCatched("w") ? args.GetDoubleValue("w") : 0;
    int ringPoints = args.IsCatched("ring_points") ? args.GetIntValue("ring_points", 0) : 3;
    if (ringPoints < 1)
    {
        std::cerr << "WARNING: --ring_points must be >= 1; using 3." << std::endl;
        ringPoints = 3;
    }

    if (wave <= 0)
    {
        std::cerr << "WARNING: wavelength is not positive; r_eq is reported but k_eq is set to 0." << std::endl;
    }

    // Enable absorption automatically when Im(ri) != 0
    bool isAbs = args.IsCatched("abs") || (im != 0);
    ::complex refrIndex = ::complex(re, isAbs ? im : 0);

    // TODO: AggregateBuilder

    Particle *particle = nullptr;


    additionalSummary += "Particle: ";

    if (args.IsCatched("pf"))
    {
        std::string filename = args.GetStringValue("pf", 0);
        particle = new Particle();
        particle->SetFromFile(filename);
        particle->SetRefractiveIndex(refrIndex);

        double origDMax = particle->MaximalDimention();
        additionalSummary += "from file: " + filename + "\n";
        additionalSummary += "\tOriginal Dmax: " + std::to_string(origDMax);

        if (args.IsCatched("rs"))
        {
//            std::string dimension = args.GetStringValue("resize", 0);
            double newSize = args.GetDoubleValue("rs", 0);

//            if (dimension == "d")
//            {
                particle->Resize(newSize);
//            }
//            else if (dimension == "v")
//            {
//                double oldV = particle->Volume();
//                double ratio = newSize/oldV;
//                ratio = pow(ratio, 1.0/3);
//                particle->Scale(ratio);
//            }
        }

        double newDMax = particle->MaximalDimention();
        additionalSummary += ", new Dmax: " + std::to_string(newDMax)
                + ", resize factor: " + std::to_string(newDMax/origDMax) + '\n';
    }
    else
    {
        ParticleType type = (ParticleType)args.GetIntValue("p", 0);
        double height;
        double diameter;

        double sup;
        int num;

        switch (type)
        {
        case ParticleType::Hexagonal:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            particle = new Hexagonal(refrIndex, diameter, height);
            break;
        case ParticleType::Bullet:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            sup = (diameter*sqrt(3)*tan(DegToRad(62)))/4;
            particle = new Bullet(refrIndex, diameter, height, sup);
            break;
        case ParticleType::BulletRosette:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            sup = (args.GetArgNumber("p") == 4)
                    ? args.GetDoubleValue("p", 3)
                    : diameter*sqrt(3)*tan(DegToRad(62))/4;
            particle = new BulletRosette(refrIndex, diameter, height, sup);
            break;
        case ParticleType::Droxtal:
            sup = args.GetDoubleValue("p", 3);
            particle = new Droxtal(refrIndex, DegToRad(32.35), DegToRad(71.81), sup);
            break;
//		case ParticleType::TiltedHexagonal:
//			sup = parser.argToValue<double>(vec[3]);
//			particle = new TiltedHexagonal(r, hh, ri, sup);
//			break;
        case ParticleType::ConcaveHexagonal:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            sup = args.GetDoubleValue("p", 3);
            particle = new ConcaveHexagonal(refrIndex, diameter, height, sup);
            break;
        case ParticleType::HexagonalAggregate:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            num = args.GetIntValue("p", 3);
            particle = new HexagonalAggregate(refrIndex, diameter, height, num);
            break;
        case ParticleType::CertainAggregate:
            sup = args.GetDoubleValue("p", 1);
            particle = new CertainAggregate(refrIndex, sup);
            break;
        default:
            assert(false && "ERROR! Incorrect type of particle.");
            break;
        }
    }

    if (args.IsCatched("k_eq"))
    {
        double targetKEq = args.GetDoubleValue("k_eq", 0);
        if (targetKEq <= 0)
        {
            std::cerr << "ERROR: --k_eq must be positive." << std::endl;
            return 1;
        }
        if (wave <= 0)
        {
            std::cerr << "ERROR: --k_eq requires positive wavelength (-w)." << std::endl;
            return 1;
        }

        double currentVolume = particle->Volume();
        double currentReq = (currentVolume > 0)
            ? pow(3.0 * currentVolume / (4.0 * M_PI), 1.0 / 3.0)
            : 0.0;
        if (currentReq <= DBL_EPSILON)
        {
            std::cerr << "ERROR: cannot apply --k_eq because current equivalent radius is zero." << std::endl;
            return 1;
        }

        double targetReq = targetKEq * wave / (2.0 * M_PI);
        double ratio = targetReq / currentReq;
        double oldDMax = particle->MaximalDimention();
        particle->Resize(oldDMax * ratio);
        additionalSummary += "\tResized by k_eq: target " + std::to_string(targetKEq)
                + ", old r_eq: " + std::to_string(currentReq)
                + ", new Dmax: " + std::to_string(particle->MaximalDimention())
                + ", resize factor: " + std::to_string(ratio) + "\n";
    }

    additionalSummary += "\tRefractive index: " + to_string(re);

    if (fabs(im) > FLT_EPSILON)
    {
        additionalSummary +=  " + i\n";
    }

    additionalSummary += "\n";

    particle->Output("particle_for_check.dat");
    double particleArea = particle->Area();
    double particleVolume = particle->Volume();
    double rEq = (particleVolume > 0)
        ? pow(3.0 * particleVolume / (4.0 * M_PI), 1.0 / 3.0)
        : 0.0;
    double kEq = (wave > 0) ? (2.0 * M_PI * rEq / wave) : 0.0;
    additionalSummary += "\tArea:" + to_string(particleArea) + "\n";
    additionalSummary += "\tVolume:" + to_string(particleVolume) + "\n";
    additionalSummary += "\tr_eq:" + to_string(rEq) + "\n";
    additionalSummary += "\tk_eq:" + to_string(kEq)
            + " (2*pi*r_eq/lambda)\n\n";

    int reflNum = args.IsCatched("n") ? (int)args.GetDoubleValue("n") : 6; // default n=6
    additionalSummary += "Number of secondary reflections: " + to_string(reflNum) + "\n";

    string dirName;
    if (args.IsCatched("o"))
    {
        dirName = args.GetStringValue("o");
    }
    else
    {
        // Auto-generate unique output name: results/run_YYYYMMDD_HHMMSS
        time_t now = time(nullptr);
        struct tm *t = localtime(&now);
        char ts[20];
        strftime(ts, sizeof(ts), "%Y%m%d_%H%M%S", t);
        string outDir = "results";
        mkdir(outDir.c_str(), 0755);
        dirName = outDir + "/run_" + ts;
        cout << "Output: " << dirName << endl;
    }
    size_t pos = 0;

    while ((pos = dirName.find('%', pos)) != string::npos)
    {
        size_t start = pos;
        ++pos;
        string key;

        while (dirName[pos] != '_' && pos < dirName.size())
        {
            key += dirName[pos++];
        }

        if (key.size())
        {
            string val;
            int i = 0;

            if (isdigit(key[0]))
            {
                i = (int)key[0]-48;
                key = key.substr(1);
            }

            val += args.GetStringValue(key, i);

            dirName.replace(start, pos-start, val);
            pos = start + val.size();
        }
    }

    // Create output folder once, before any tracer runs
    {
        string dir;
        if (mpi_rank == 0) dir = CreateFolder(dirName);
#ifdef USE_MPI
        if (mpi_size > 1) {
            int dirLen = dir.size();
            MPI_Bcast(&dirLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
            dir.resize(dirLen);
            MPI_Bcast(&dir[0], dirLen, MPI_CHAR, 0, MPI_COMM_WORLD);
            int nameLen = dirName.size();
            MPI_Bcast(&nameLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
            dirName.resize(nameLen);
            MPI_Bcast(&dirName[0], nameLen, MPI_CHAR, 0, MPI_COMM_WORLD);
        }
#endif
        dirName = dir + dirName;
    }

    bool isOutputGroups = args.IsCatched("gr");
    additionalSummary += "Wavelength (um): " + to_string(wave) + "\n";

    if (args.IsCatched("forced_nonconvex"))
    {
        particle->isConcave = true;
    }

    if (args.IsCatched("forced_convex"))
    {
        particle->isConcave = false;
    }

    if (args.IsCatched("tr"))
    {
        string trackFileName = args.GetStringValue("tr");
        trackGroups.ImportTracks(particle->nFacets, trackFileName);
        trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
        trackGroups.shouldOutputGroups = args.IsCatched("gr");
    }
    else
    {

    }

    int nTheta = 1;
    if (args.IsCatched("grid"))
    {
        int gridArgs = args.GetArgNumber("grid");
        nTheta = (gridArgs == 3)
            ? (int)args.GetDoubleValue("grid", 2)
            : (int)args.GetDoubleValue("grid", 3);
    }

    additionalSummary += "Method: ";

    if (args.IsCatched("po"))
    {
        additionalSummary += "Physical optics";

        if (args.IsCatched("karczewski"))
        {
            cout << "Note: --karczewski does not affect M11 (same Frobenius norm as RotateJones)." << endl;
            cout << "      Experimental flag for M33/M34/M44 research only." << endl;
        }

        if (args.IsCatched("fixed"))
        {
            additionalSummary += ", fixed orientation\n\n";

            ScatteringRange bsCone = SetConus(args);

            double beta  = args.GetDoubleValue("fixed", 0);
            double gamma = args.GetDoubleValue("fixed", 1);

            TracerPO tracer(particle, reflNum, dirName);
            tracer.m_summary = additionalSummary;
            ApplyTraceCutoffOptions(args, tracer.m_scattering);

            HandlerPO *handler = new HandlerPO(particle, &tracer.m_incidentLight,
                                               nTheta, wave);
            if (args.IsCatched("r"))
            {
                tracer.m_scattering->restriction = args.GetDoubleValue("r", 0);
            }

            handler->isCoh = !args.IsCatched("incoh");
            handler->SetGpuEnabled(useGpu);
            handler->SetFftEnabled(useFft);
            ApplyPOOutputOptions(args, handler, wave);
            handler->m_legacySign = args.IsCatched("legacy_sign");
            handler->useKarczewski = args.IsCatched("karczewski");
            handler->outputJones = args.IsCatched("jones");
            ApplyNphiOverride(args, bsCone);
            handler->SetScatteringSphere(bsCone);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);
            ApplyAbsorptionPointOption(args, handler);
            ApplyBeamCutoffOptions(args, handler);

            if (args.IsCatched("log"))
            {
                tracer.m_logTime = args.GetIntValue("log");
            }

            tracer.shadowOff = args.IsCatched("shadow_off");
            tracer.SetIsOutputGroups(isOutputGroups);
            tracer.SetHandler(handler);
            tracer.TraceFixed(beta, gamma);
        }
        else if (args.IsCatched("oldauto"))
        {
            if (args.IsCatched("multikeq") || args.IsCatched("multikeq_list"))
            {
                double KmaxRef = 0.0;
                if (args.IsCatched("multikeq_list"))
                {
                    std::vector<double> listed = ReadSizeList(args.GetStringValue("multikeq_list", 0));
                    KmaxRef = *std::max_element(listed.begin(), listed.end());
                }
                else
                {
                    KmaxRef = args.GetDoubleValue("multikeq", 1);
                }
                if (KmaxRef <= 0 || wave <= 0)
                {
                    std::cerr << "ERROR: --multikeq requires positive Kmax/list max and wavelength." << std::endl;
                    return 1;
                }
                double v = particle->Volume();
                double req = (v > 0) ? pow(3.0 * v / (4.0 * M_PI), 1.0 / 3.0) : 0.0;
                if (req <= DBL_EPSILON)
                {
                    std::cerr << "ERROR: cannot apply --multikeq because current equivalent radius is zero." << std::endl;
                    return 1;
                }
                double targetReq = KmaxRef * wave / (2.0 * M_PI);
                double ratio = targetReq / req;
                particle->Resize(particle->MaximalDimention() * ratio);
                additionalSummary += "\tShared multikeq reference resize: k_eq="
                    + std::to_string(KmaxRef)
                    + ", Dmax=" + std::to_string(particle->MaximalDimention())
                    + ", resize factor=" + std::to_string(ratio) + "\n";
            }
            else if (args.IsCatched("multigrid"))
            {
                double DmaxRef = args.GetDoubleValue("multigrid", 1);
                if (DmaxRef <= 0)
                {
                    std::cerr << "ERROR: --multigrid requires positive Dmax." << std::endl;
                    return 1;
                }
                double oldD = particle->MaximalDimention();
                particle->Resize(DmaxRef);
                additionalSummary += "\tShared multigrid reference resize: Dmax="
                    + std::to_string(DmaxRef)
                    + ", resize factor=" + std::to_string(particle->MaximalDimention() / oldD)
                    + "\n";
            }

            // Physics-based orientation grid from diffraction angular step
            // --oldauto DIV: DIV = 2,4,8 (divisor for full diffraction-limited grid)
            int div = args.GetIntValue("oldauto", 0);
            if (div < 1) div = 8;

            double L = particle->MaximalDimention();
            double D_particle = L; // will be overridden below
            // D = 6.96 * sqrt(L) from Excel formula
            // But actual D is from -p argument, use particle geometry
            // L is the height (first -p arg), D is the diameter (second -p arg)
            // particle->MaximalDimention() returns max(L,D)

            // Diffraction angular step: Δθ = 0.69 * λ / L * (180/π)
            double delta_theta_deg = 0.69 * wave / L * (180.0 / M_PI);
            int points_per_ring = ringPoints;

            // Angular step for orientation grid
            double orient_step = delta_theta_deg / points_per_ring;

            // Get symmetry from particle
            double betaSym_deg = RadToDeg(particle->GetSymmetry().beta);
            double gammaSym_deg = RadToDeg(particle->GetSymmetry().gamma);
            bool mirrorGamma = args.IsCatched("mirror_gamma");
            double gammaRange_deg = mirrorGamma ? 0.5 * gammaSym_deg : gammaSym_deg;

            // Full grid: N = sym_range / step
            int N_beta_full = (int)ceil(betaSym_deg / orient_step);
            int N_gamma_full = (int)ceil(gammaRange_deg / orient_step);

            // Divide by div factor, round UP
            int N_beta = (N_beta_full + div - 1) / div;
            int N_gamma = (N_gamma_full + div - 1) / div;
            if (N_beta < 3) N_beta = 3;
            if (N_gamma < 3) N_gamma = 3;
            if (mirrorGamma && (N_gamma % 2 != 0))
            {
                int old = N_gamma;
                ++N_gamma;
                std::cerr << "WARNING: --mirror_gamma rounded oldauto gamma "
                          << "points to an even half-domain count: "
                          << old << " -> " << N_gamma << std::endl;
            }

            // N_phi priority: --nphi > --grid > default 360.
            // --grid has two forms: radius Nphi Ntheta, or theta1 theta2 Nphi Ntheta.
            int N_phi = 360;
            if (args.IsCatched("nphi"))
            {
                N_phi = args.GetIntValue("nphi", 0);
            }
            else if (args.IsCatched("grid"))
            {
                int gridArgs = args.GetArgNumber("grid");
                N_phi = (gridArgs == 3)
                    ? (int)args.GetDoubleValue("grid", 1)
                    : (int)args.GetDoubleValue("grid", 2);
            }
            N_phi = MirrorSafePhiCount(args, N_phi);

            cout << "=== oldauto (physics-based) ===" << endl;
            cout << "  Dmax=" << L << " um, lambda=" << wave << " um" << endl;
            cout << "  Diffraction step: " << delta_theta_deg << " deg" << endl;
            cout << "  Orient step: " << orient_step << " deg (" << points_per_ring
                 << " pts/ring)" << endl;
            cout << "  Full grid: " << N_beta_full << " x " << N_gamma_full
                 << " = " << (long long)N_beta_full * N_gamma_full << endl;
            cout << "  div" << div << ": " << N_beta << " x " << N_gamma
                 << " = " << N_beta * N_gamma << " orientations" << endl;
            if (mirrorGamma)
                cout << "  Gamma mirror symmetry: using 0.." << gammaRange_deg
                     << " deg instead of 0.." << gammaSym_deg << " deg" << endl;
            cout << "  N_phi=" << N_phi << endl;
            cout << "  n=" << reflNum << endl;

            additionalSummary += ", oldauto div" + to_string(div) + "\n\n";
            if (mirrorGamma)
                additionalSummary += "\tGamma mirror symmetry: 0.."
                    + to_string(gammaRange_deg) + " deg\n";

            // Build conus
            ScatteringRange conus = (args.IsCatched("grid") || args.IsCatched("tgrid"))
                ? SetConus(args)
                : ScatteringRange(0, M_PI, N_phi, 1);

            // Apply theta grid (only if no explicit --grid or --tgrid)
            if (args.IsCatched("tgrid") || args.IsCatched("grid")) {
                // tgrid or grid already loaded in SetConus above
            } else {
                ApplyAutoThetaGrid(conus, L, wave);
            }

            // Set N_phi (only if --grid didn't set it explicitly)
            if (!args.IsCatched("grid")) {
                SetRangeNphi(args, conus, N_phi);
            }

            TracerPOTotal *tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { if (args.IsCatched("log")) tpt->m_logTime = args.GetIntValue("log"); tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); tpt->m_enableCheckpoint = args.IsCatched("checkpoint"); tpt->m_fastPoleGamma = args.IsCatched("pole"); tpt->m_mirrorGamma = args.IsCatched("mirror_gamma"); tpt->m_sobolChunkSize = args.IsCatched("chunk") ? std::max(1, args.GetIntValue("chunk", 0)) : 0; tpt->m_ringPoints = ringPoints; } }
            tracer->m_scattering->m_wave = wave;
            ApplyTraceCutoffOptions(args, tracer->m_scattering);
            tracer->shadowOff = args.IsCatched("shadow_off");
            trackGroups.push_back(TrackGroup());
            HandlerPOTotal *handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                         nTheta, wave);

            cout << additionalSummary;
            tracer->m_summary = additionalSummary;

            handler->isCoh = !args.IsCatched("incoh");
            handler->SetGpuEnabled(useGpu);
            handler->SetFftEnabled(useFft);
            ApplyPOOutputOptions(args, handler, wave);
            handler->m_legacySign = args.IsCatched("legacy_sign");
            handler->useKarczewski = args.IsCatched("karczewski");
            ApplyNphiOverride(args, conus);
            handler->SetScatteringSphere(conus);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);
            ApplyAbsorptionPointOption(args, handler);
            ApplyBeamCutoffOptions(args, handler);

            tracer->SetIsOutputGroups(isOutputGroups);
            tracer->SetHandler(handler);

            // Use --random with computed N_beta, N_gamma
            AngleRange betaRange(0, particle->GetSymmetry().beta, N_beta);
            AngleRange gammaRange(0, DegToRad(gammaRange_deg), N_gamma);

            if (args.IsCatched("multikeq") || args.IsCatched("multikeq_list"))
            {
                std::vector<double> kSizes;
                std::string rangeText;
                if (args.IsCatched("multikeq_list"))
                {
                    std::string listFile = args.GetStringValue("multikeq_list", 0);
                    kSizes = ReadSizeList(listFile);
                    double Kmin = *std::min_element(kSizes.begin(), kSizes.end());
                    double Kmax = *std::max_element(kSizes.begin(), kSizes.end());
                    rangeText = std::to_string(Kmin) + ".." + std::to_string(Kmax)
                        + " from " + listFile;
                }
                else
                {
                    double Kmin = args.GetDoubleValue("multikeq", 0);
                    double Kmax = args.GetDoubleValue("multikeq", 1);
                    int nSizes = args.GetIntValue("multikeq", 2);
                    if (Kmin <= 0 || Kmax <= 0)
                    {
                        std::cerr << "ERROR: --multikeq values must be positive." << std::endl;
                        return 1;
                    }
                    kSizes = GenerateLogSizes(Kmin, Kmax, nSizes);
                    rangeText = std::to_string(Kmin) + ".." + std::to_string(Kmax);
                }
                if (kSizes.empty())
                {
                    std::cerr << "ERROR: --multikeq has no sizes." << std::endl;
                    return 1;
                }
                double currentVolumeForMulti = particle->Volume();
                double currentReqForMulti = (currentVolumeForMulti > 0)
                    ? pow(3.0 * currentVolumeForMulti / (4.0 * M_PI), 1.0 / 3.0)
                    : 0.0;
                double kRef = (wave > 0 && currentReqForMulti > 0)
                    ? (2.0 * M_PI * currentReqForMulti / wave) : 0.0;
                double xRef = M_PI * particle->MaximalDimention() / wave;
                if (kRef <= 0 || xRef <= 0)
                {
                    std::cerr << "ERROR: --multikeq requires positive wavelength and particle volume." << std::endl;
                    return 1;
                }
                std::vector<double> xSizes;
                std::vector<std::string> labels;
                xSizes.reserve(kSizes.size());
                labels.reserve(kSizes.size());
                for (double k : kSizes)
                {
                    xSizes.push_back(xRef * (k / kRef));
                    labels.push_back("keq" + SizeLabel(k));
                }
                std::cout << "Shared oldauto multikeq: k_eq " << rangeText
                          << " (" << kSizes.size() << " sizes), reference k_eq="
                          << kRef << std::endl;
                tracer->TraceRandomMultiSize(betaRange, gammaRange, xSizes, labels);
            }
            else if (args.IsCatched("multigrid"))
            {
                double Dmin = args.GetDoubleValue("multigrid", 0);
                double Dmax_mg = args.GetDoubleValue("multigrid", 1);
                int nSizes = args.GetIntValue("multigrid", 2);
                if (Dmin <= 0 || Dmax_mg <= 0)
                {
                    std::cerr << "ERROR: --multigrid values must be positive." << std::endl;
                    return 1;
                }
                std::vector<double> dSizes = GenerateLogSizes(Dmin, Dmax_mg, nSizes);
                double Dref = particle->MaximalDimention();
                double xRef = M_PI * Dref / wave;
                std::vector<double> xSizes;
                std::vector<std::string> labels;
                xSizes.reserve(dSizes.size());
                labels.reserve(dSizes.size());
                for (double d : dSizes)
                {
                    xSizes.push_back(xRef * (d / Dref));
                    labels.push_back("D" + SizeLabel(d));
                }
                std::cout << "Shared oldauto multigrid: D " << Dmin << ".." << Dmax_mg
                          << " (" << dSizes.size() << " log sizes), reference Dmax="
                          << Dref << std::endl;
                tracer->TraceRandomMultiSize(betaRange, gammaRange, xSizes, labels);
            }
            else
            {
                tracer->TraceRandom(betaRange, gammaRange);
            }
            delete handler;
        }
        else if (args.IsCatched("random"))
        {
            // Warn if --auto/--adaptive also specified (random takes priority)
            if (args.IsCatched("auto") || args.IsCatched("autofull")
                || args.IsCatched("oldautofull") || args.IsCatched("adaptive"))
                std::cerr << "WARNING: --random overrides --auto/--adaptive. Use --sobol instead." << std::endl;
            additionalSummary += ", random orientation\n\n";
            AngleRange beta = GetRange(args, "b", particle);
            AngleRange gamma = GetRange(args, "g", particle);
            MakeMirrorGammaRangeConsistent(args, particle, gamma);

            HandlerPO *handler;

            if (args.IsCatched("point"))
            {
                std::cerr << "ERROR: --point is not supported in the optimized PO --random path. "
                          << "Use --grid with a narrow backscatter cone or remove --point." << std::endl;
                exit(1);
//                 TracerBackScatterPoint tracer(particle, reflNum, dirName);
//                 tracer.m_scattering->m_wave = wave;
//                 ApplyTraceCutoffOptions(args, tracer.m_scattering);
//                 if (args.IsCatched("r"))
//                 {
//                     tracer.m_scattering->restriction = args.GetDoubleValue("r", 0);
//                 }
//                 cout << additionalSummary;
//                 tracer.m_summary = additionalSummary;

//                 handler = new HandlerBackScatterPoint(particle, &tracer.m_incidentLight,
//                                                       nTheta, wave);
//                 double normIndex = gamma.step/gamma.norm;
// //                handler->isCoh = !args.IsCatched("incoh");
//              handler->m_legacySign = args.IsCatched("legacy_sign");
//                 handler->SetNormIndex(normIndex);
//                 handler->SetTracks(&trackGroups);
//                 trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
//                 handler->SetAbsorptionAccounting(isAbs);

//                 if (args.IsCatched("filter"))
//                 {
//                     handler->SetBackScatteringConus(DegToRad(args.GetDoubleValue("filter")));
//                 }

//                 if (args.IsCatched("log"))
//                 {
//                     tracer.m_logTime = args.GetIntValue("log");
//                 }
//                 tracer.SetIsOutputGroups(isOutputGroups);
//                 tracer.SetHandler(handler);
//                 tracer.TraceRandom(beta, gamma);
            }
            else
            {
                TracerPO *tracer;
                ScatteringRange conus = SetConus(args);

                if (args.IsCatched("all"))
                {
                    tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { if (args.IsCatched("log")) tpt->m_logTime = args.GetIntValue("log"); tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); tpt->m_enableCheckpoint = args.IsCatched("checkpoint"); tpt->m_fastPoleGamma = args.IsCatched("pole"); tpt->m_mirrorGamma = args.IsCatched("mirror_gamma"); tpt->m_sobolChunkSize = args.IsCatched("chunk") ? std::max(1, args.GetIntValue("chunk", 0)) : 0; tpt->m_ringPoints = ringPoints; } }
                    tracer->m_scattering->m_wave = wave;
            ApplyTraceCutoffOptions(args, tracer->m_scattering);
                    if (args.IsCatched("r"))
                    {
                        tracer->m_scattering->restriction = args.GetDoubleValue("r", 0);
                    }
                    tracer->shadowOff = args.IsCatched("shadow_off");

                    trackGroups.push_back(TrackGroup());
                    handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                                 nTheta, wave);
                    handler->normIndexGamma = gamma.step/gamma.norm;

                    if (args.IsCatched("filter"))
                    {
                        handler->SetBackScatteringConus(DegToRad(args.GetDoubleValue("filter")));
                    }
                }
                else
                {
                    // Use TracerPOTotal for OpenMP + batched sincos acceleration
                    tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { if (args.IsCatched("log")) tpt->m_logTime = args.GetIntValue("log"); tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); tpt->m_enableCheckpoint = args.IsCatched("checkpoint"); tpt->m_fastPoleGamma = args.IsCatched("pole"); tpt->m_mirrorGamma = args.IsCatched("mirror_gamma"); tpt->m_sobolChunkSize = args.IsCatched("chunk") ? std::max(1, args.GetIntValue("chunk", 0)) : 0; tpt->m_ringPoints = ringPoints; } }
                    tracer->m_scattering->m_wave = wave;
            ApplyTraceCutoffOptions(args, tracer->m_scattering);
                    tracer->shadowOff = args.IsCatched("shadow_off");
                    if (args.IsCatched("r"))
                    {
                        tracer->m_scattering->restriction = args.GetDoubleValue("r", 0);
                    }

                    trackGroups.push_back(TrackGroup());
                    handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                                 nTheta, wave);

                    if (args.IsCatched("filter"))
                    {
                        handler->SetBackScatteringConus(DegToRad(args.GetDoubleValue("filter")));
                    }
                }

                cout << additionalSummary;
                tracer->m_summary = additionalSummary;

                handler->isCoh = !args.IsCatched("incoh");
                handler->SetGpuEnabled(useGpu);
                handler->SetFftEnabled(useFft);
                ApplyPOOutputOptions(args, handler, wave);
                handler->m_legacySign = args.IsCatched("legacy_sign");
                handler->useKarczewski = args.IsCatched("karczewski");
                ApplyNphiOverride(args, conus);
                handler->SetScatteringSphere(conus);
                handler->SetTracks(&trackGroups);
                handler->SetAbsorptionAccounting(isAbs);
                ApplyAbsorptionPointOption(args, handler);
                ApplyBeamCutoffOptions(args, handler);

                if (args.IsCatched("log"))
                {
                    tracer->m_logTime = args.GetIntValue("log");
                }
                tracer->SetIsOutputGroups(isOutputGroups);
                tracer->SetHandler(handler);

                // Adaptive theta grid for --random
                if (args.IsCatched("auto_tgrid") && !args.IsCatched("grid") && !args.IsCatched("tgrid"))
                {
                    TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer);
                    if (tpt) {
                        double tgridEps = args.GetDoubleValue("auto_tgrid", 0);
                        if (tgridEps <= 0) tgridEps = 0.05;
                        double betaSym2 = particle->GetSymmetry().beta;
                        double gammaSym2 = particle->GetSymmetry().gamma;
                        // Probe count from div16 physics estimate
                        double Dmax2 = particle->MaximalDimention();
                        double dd2 = 0.69 * wave / Dmax2 * (180.0 / M_PI) / ringPoints;
                        int nb2 = std::max(1, (int)(RadToDeg(betaSym2) / dd2 / 16));
                        int ng2 = std::max(1, (int)(RadToDeg(gammaSym2) / dd2 / 16));
                        int nProbe = nb2 * ng2;
                        int p2 = 1; while (p2 * 2 <= nProbe) p2 *= 2;
                        nProbe = std::max(64, p2);
                        tpt->TraceAdaptiveTheta(nProbe, betaSym2, gammaSym2, tgridEps, 8, true);
                        // Grid is now set in handler->m_sphere. TraceRandom will use it.
                    }
                }

                tracer->TraceRandom(beta, gamma);
            }

            delete handler;
        }
        else if (args.IsCatched("montecarlo"))
        {
            additionalSummary += ", Monte Carlo method\n\n";
            // montecarlo doesn't use --random, get ranges from particle symmetry
            double betaMax = particle->GetSymmetry().beta;
            double gammaMax = particle->GetSymmetry().gamma;
            if (args.IsCatched("mirror_gamma"))
                gammaMax *= 0.5;
            int nOr = args.GetIntValue("montecarlo", 0);
            if (args.IsCatched("mirror_gamma") && (nOr % 2 != 0))
            {
                std::cerr << "WARNING: --mirror_gamma rounded montecarlo "
                          << "gamma grid count to an even half-domain count: "
                          << nOr << " -> " << (nOr + 1) << std::endl;
                ++nOr;
            }
            AngleRange beta(0, betaMax, nOr);
            AngleRange gamma(0, gammaMax, nOr);

            HandlerPO *handler;

            TracerPOTotal *tracer;
            ScatteringRange conus = SetConus(args);

            tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { if (args.IsCatched("log")) tpt->m_logTime = args.GetIntValue("log"); tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); tpt->m_enableCheckpoint = args.IsCatched("checkpoint"); tpt->m_fastPoleGamma = args.IsCatched("pole"); tpt->m_mirrorGamma = args.IsCatched("mirror_gamma"); tpt->m_sobolChunkSize = args.IsCatched("chunk") ? std::max(1, args.GetIntValue("chunk", 0)) : 0; tpt->m_ringPoints = ringPoints; } }
            tracer->m_scattering->m_wave = wave;
            ApplyTraceCutoffOptions(args, tracer->m_scattering);
            tracer->shadowOff = args.IsCatched("shadow_off");
            if (args.IsCatched("r"))
                tracer->m_scattering->restriction = args.GetDoubleValue("r", 0);
            trackGroups.push_back(TrackGroup());
            handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                         nTheta, wave);

            cout << additionalSummary;
            tracer->m_summary = additionalSummary;

            handler->isCoh = !args.IsCatched("incoh");
            handler->SetGpuEnabled(useGpu);
            handler->SetFftEnabled(useFft);
            ApplyPOOutputOptions(args, handler, wave);
            handler->m_legacySign = args.IsCatched("legacy_sign");
            ApplyNphiOverride(args, conus);
            handler->SetScatteringSphere(conus);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);
            ApplyAbsorptionPointOption(args, handler);
            ApplyBeamCutoffOptions(args, handler);

            if (args.GetArgNumber("n") == 3 &&
                args.GetStringValue("n", 2) == "fixed")
            {
                handler->m_fixedItr = reflNum;
            }

            tracer->SetIsOutputGroups(isOutputGroups);
            tracer->SetHandler(handler);

            tracer->TraceMonteCarlo(beta, gamma, nOr);
        }
        else if (args.IsCatched("orientfile"))
        {
            additionalSummary += ", orientations from file\n\n";

            TracerPOTotal *tracer;
            ScatteringRange conus = SetConus(args);

            tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { if (args.IsCatched("log")) tpt->m_logTime = args.GetIntValue("log"); tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); tpt->m_enableCheckpoint = args.IsCatched("checkpoint"); tpt->m_fastPoleGamma = args.IsCatched("pole"); tpt->m_mirrorGamma = args.IsCatched("mirror_gamma"); tpt->m_sobolChunkSize = args.IsCatched("chunk") ? std::max(1, args.GetIntValue("chunk", 0)) : 0; tpt->m_ringPoints = ringPoints; } }
            tracer->m_scattering->m_wave = wave;
            ApplyTraceCutoffOptions(args, tracer->m_scattering);
            tracer->shadowOff = args.IsCatched("shadow_off");
            trackGroups.push_back(TrackGroup());
            HandlerPOTotal *handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                         nTheta, wave);

            cout << additionalSummary;
            tracer->m_summary = additionalSummary;

            handler->isCoh = !args.IsCatched("incoh");
            handler->SetGpuEnabled(useGpu);
            handler->SetFftEnabled(useFft);
            ApplyPOOutputOptions(args, handler, wave);
            handler->m_legacySign = args.IsCatched("legacy_sign");
            handler->useKarczewski = args.IsCatched("karczewski");
            ApplyNphiOverride(args, conus);
            handler->SetScatteringSphere(conus);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);
            ApplyAbsorptionPointOption(args, handler);
            ApplyBeamCutoffOptions(args, handler);

            tracer->SetIsOutputGroups(isOutputGroups);
            tracer->SetHandler(handler);

            std::string orientFileName = args.GetStringValue("orientfile", 0);

            tracer->TraceFromFile(orientFileName);

            delete handler;
        }
        else if (args.IsCatched("sobol") || args.IsCatched("sobol_seed")
              || args.IsCatched("sobol_ring")
              || args.IsCatched("hammersley") || args.IsCatched("lattice")
              || args.IsCatched("lattice_z")
              || args.IsCatched("euler_quad")
              || args.IsCatched("adaptive") || args.IsCatched("auto")
              || args.IsCatched("autofull") || args.IsCatched("oldautofull"))
        {
            // --auto EPS implies --adaptive EPS + --auto_tgrid + --auto_phi
            bool isOldAutoFull = args.IsCatched("oldautofull");
            bool isAuto = args.IsCatched("auto") || args.IsCatched("autofull")
                || isOldAutoFull;
            bool isAutoFull = args.IsCatched("autofull") || isOldAutoFull;
            bool isAdaptive = args.IsCatched("adaptive") || isAuto;
            bool isSobolSeed = args.IsCatched("sobol_seed");
            bool isSobolRing = args.IsCatched("sobol_ring");
            bool isHammersley = args.IsCatched("hammersley");
            bool isLattice = args.IsCatched("lattice") || args.IsCatched("lattice_z");
            bool isEulerQuad = args.IsCatched("euler_quad");
            bool oldAutoFullMultiSize = isOldAutoFull
                && (args.IsCatched("multikeq") || args.IsCatched("multikeq_list")
                    || args.IsCatched("multigrid"));
            double oldAutoFullKRef = 0.0;
            double oldAutoFullDRef = 0.0;
            if (oldAutoFullMultiSize)
            {
                if (args.IsCatched("multikeq") || args.IsCatched("multikeq_list"))
                {
                    if (args.IsCatched("multikeq_list"))
                    {
                        std::vector<double> listed =
                            ReadSizeList(args.GetStringValue("multikeq_list", 0));
                        if (listed.empty())
                        {
                            std::cerr << "ERROR: --multikeq_list has no sizes." << std::endl;
                            return 1;
                        }
                        oldAutoFullKRef =
                            *std::max_element(listed.begin(), listed.end());
                    }
                    else
                    {
                        oldAutoFullKRef = args.GetDoubleValue("multikeq", 1);
                    }
                    if (oldAutoFullKRef <= 0 || wave <= 0)
                    {
                        std::cerr << "ERROR: --oldautofull multikeq requires positive Kmax/list max and wavelength." << std::endl;
                        return 1;
                    }
                    double v = particle->Volume();
                    double req = (v > 0)
                        ? pow(3.0 * v / (4.0 * M_PI), 1.0 / 3.0) : 0.0;
                    if (req <= DBL_EPSILON)
                    {
                        std::cerr << "ERROR: cannot apply --oldautofull multikeq because current equivalent radius is zero." << std::endl;
                        return 1;
                    }
                    double targetReq = oldAutoFullKRef * wave / (2.0 * M_PI);
                    double ratio = targetReq / req;
                    particle->Resize(particle->MaximalDimention() * ratio);
                    additionalSummary += "\tOldautofull multikeq reference resize: k_eq="
                        + std::to_string(oldAutoFullKRef)
                        + ", Dmax=" + std::to_string(particle->MaximalDimention())
                        + ", resize factor=" + std::to_string(ratio) + "\n";
                }
                else
                {
                    oldAutoFullDRef = args.GetDoubleValue("multigrid", 1);
                    if (oldAutoFullDRef <= 0)
                    {
                        std::cerr << "ERROR: --oldautofull multigrid requires positive Dmax." << std::endl;
                        return 1;
                    }
                    double oldD = particle->MaximalDimention();
                    particle->Resize(oldAutoFullDRef);
                    additionalSummary += "\tOldautofull multigrid reference resize: Dmax="
                        + std::to_string(oldAutoFullDRef)
                        + ", resize factor=" + std::to_string(particle->MaximalDimention() / oldD)
                        + "\n";
                }
            }
            if (isAdaptive)
                additionalSummary += ", adaptive Sobol\n\n";
            else if (isSobolSeed)
                additionalSummary += ", Sobol nested Owen seed\n\n";
            else if (isSobolRing)
                additionalSummary += ", Sobol beta x gamma ring\n\n";
            else if (isHammersley)
                additionalSummary += ", Hammersley orientations\n\n";
            else if (isLattice)
                additionalSummary += ", rank-1 lattice orientations\n\n";
            else if (isEulerQuad)
                additionalSummary += ", Euler high-order quadrature\n\n";
            else
                additionalSummary += ", Sobol quasi-random\n\n";

            TracerPOTotal *tracer;
            // --grid optional for --auto/--autofull/--tgrid
            ScatteringRange conus = (args.IsCatched("grid") || args.IsCatched("tgrid"))
                ? SetConus(args)
                : ScatteringRange(0, M_PI, 1, 1);

            tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { if (args.IsCatched("log")) tpt->m_logTime = args.GetIntValue("log"); tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); tpt->m_enableCheckpoint = args.IsCatched("checkpoint"); tpt->m_fastPoleGamma = args.IsCatched("pole"); tpt->m_mirrorGamma = args.IsCatched("mirror_gamma"); tpt->m_oldAutoFullFinal = isOldAutoFull; tpt->m_sobolChunkSize = args.IsCatched("chunk") ? std::max(1, args.GetIntValue("chunk", 0)) : 0; tpt->m_ringPoints = ringPoints; } }
            tracer->m_scattering->m_wave = wave;
            ApplyTraceCutoffOptions(args, tracer->m_scattering);
            tracer->shadowOff = args.IsCatched("shadow_off");
            trackGroups.push_back(TrackGroup());
            HandlerPOTotal *handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                         nTheta, wave);

            // Auto theta grid: use adaptive bisection for both --auto and --auto_tgrid
            // (no static zones anymore)

            // Apply auto_phi if requested (or implied by --auto)
            // --grid phi has priority: don't overwrite explicit N_phi
            if ((args.IsCatched("auto_phi") || isAuto) && !args.IsCatched("grid"))
            {
                double D = particle->MaximalDimention();
                double x = (wave > 0) ? M_PI * D / wave : 100;
                // Phi convergence depends on x (linear, validated by --autofull):
                //   x=100->72, x=300->108, x=600->168, x=1000->252, x=2000->360
                // Validated: autofull found phi=168 for x=618 (matches formula)
                // Formula: N_phi = min(360, x/5 + 48), rounded to multiple of 6
                int nPhi_raw = std::max(48, std::min(360, (int)(x / 5.0 + 48)));
                int nPhi = ((nPhi_raw + 5) / 6) * 6;

                conus.nAzimuth = nPhi;
                conus.azinuthStep = 2.0 * M_PI / nPhi;
                cout << "Auto phi: x=" << x << ", N_phi=" << nPhi << endl;
            }

            cout << additionalSummary;
            tracer->m_summary = additionalSummary;

            handler->isCoh = !args.IsCatched("incoh");
            handler->SetGpuEnabled(useGpu);
            handler->SetFftEnabled(useFft);
            ApplyPOOutputOptions(args, handler, wave);
            handler->m_legacySign = args.IsCatched("legacy_sign");
            handler->useKarczewski = args.IsCatched("karczewski");
            ApplyNphiOverride(args, conus);
            if (useFft && isAuto)
            {
                if (!HandlerPO::HasNumericFftPhiFactorOverride())
                {
                    double fftEps = isAutoFull
                        ? (isOldAutoFull ? args.GetDoubleValue("oldautofull", 0)
                                         : args.GetDoubleValue("autofull", 0))
                        : args.GetDoubleValue("auto", 0);
                    int fftFactor = HandlerPO::SelectAutoFftPhiFactor(conus.nAzimuth, fftEps);
                    handler->SetFftPhiFactor(fftFactor);
                    int directPhi = fftFactor > 1
                        ? std::max(1, conus.nAzimuth / fftFactor)
                        : conus.nAzimuth;
                    cout << (isAutoFull ? "Initial auto FFT phi factor: "
                                         : "Auto FFT phi factor: ")
                         << "N_phi=" << conus.nAzimuth
                         << ", factor=" << fftFactor
                         << ", direct N_phi=" << directPhi
                         << ", eps=" << fftEps << endl;
                }
                else
                {
                    cout << "Auto FFT phi factor: using MBS_FFT_PHI_FACTOR override" << endl;
                }
            }
            handler->SetScatteringSphere(conus);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);
            ApplyAbsorptionPointOption(args, handler);
            ApplyBeamCutoffOptions(args, handler);

            tracer->SetIsOutputGroups(isOutputGroups);
            tracer->SetHandler(handler);

            double betaSym = particle->GetSymmetry().beta;
            double gammaSym = particle->GetSymmetry().gamma;

            // Override symmetry: --sym beta_factor gamma_factor
            // e.g. --sym 2 6 means β∈[0,π/2], γ∈[0,2π/6] = 2× and 6× reduction
            // --sym 1 1 means full sphere (no symmetry reduction)
            if (args.IsCatched("sym"))
            {
                int symBeta = args.GetIntValue("sym", 0);
                int symGamma = args.GetIntValue("sym", 1);
                betaSym = M_PI / std::max(symBeta, 1);
                gammaSym = 2.0 * M_PI / std::max(symGamma, 1);
                cout << "Symmetry override: beta_sym=" << RadToDeg(betaSym)
                     << " deg (/" << symBeta << "), gamma_sym=" << RadToDeg(gammaSym)
                     << " deg (/" << symGamma << ")" << endl;
                if (isAuto)
                    std::cerr << "WARNING: --sym overrides auto-symmetry. "
                              << "Usually NOT needed with --auto (particle symmetry is auto-detected)."
                              << std::endl;
            }

            // Beam cutoff is opt-in only.  Autofull accuracy targets must not
            // silently discard beams; use --beam_cutoff* explicitly.

            if (isAutoFull)
            {
                double epsAdapt = isOldAutoFull
                    ? args.GetDoubleValue("oldautofull", 0)
                    : args.GetDoubleValue("autofull", 0);
                int maxOrientUser = args.IsCatched("maxorient")
                    ? args.GetIntValue("maxorient", 0) : 0;
                tracer->m_owenAverageSeeds = BuildOwenAverageSeeds(args);
                if (!tracer->m_owenAverageSeeds.empty())
                {
                    std::cout << "Autofull final Owen averaging seeds:";
                    for (unsigned int seed : tracer->m_owenAverageSeeds)
                        std::cout << ' ' << seed;
                    std::cout << std::endl;
                }
                if (isOldAutoFull)
                    std::cout << "Oldautofull: final orientation averaging uses"
                              << " regular oldauto beta/gamma grid" << std::endl;
                if (isOldAutoFull && args.IsCatched("mirror_gamma"))
                {
                    gammaSym *= 0.5;
                    std::cout << "Oldautofull mirror gamma: using 0.."
                              << RadToDeg(gammaSym)
                              << " deg half-domain for final regular grid"
                              << std::endl;
                }
                else if (!isOldAutoFull && args.IsCatched("mirror_gamma"))
                {
                    std::cerr << "WARNING: --mirror_gamma with --autofull/Sobol "
                              << "projects an irregular orientation set. Use "
                              << "--oldautofull for regular oldauto-style "
                              << "mirror averaging." << std::endl;
                }
                tracer->TraceAutoFull(epsAdapt, betaSym, gammaSym, maxOrientUser,
                                      particle, wave, conus, handler);
                if (oldAutoFullMultiSize)
                {
                    if (!tracer->m_lastOldAutoFullGridValid)
                    {
                        std::cerr << "ERROR: --oldautofull did not produce a regular final grid for multigrid/multikeq." << std::endl;
                        return 1;
                    }

                    AngleRange betaRange(0.0, tracer->m_lastOldAutoFullBetaSym,
                                         tracer->m_lastOldAutoFullNBeta);
                    AngleRange gammaRange(0.0, tracer->m_lastOldAutoFullGammaSym,
                                          tracer->m_lastOldAutoFullNGamma);

                    std::cout << std::endl
                              << "Oldautofull shared multi-size: reusing tuned"
                              << " reference grid n=" << tracer->m_lastOldAutoFullN
                              << ", N_phi=" << tracer->m_lastOldAutoFullNphi
                              << ", beta/gamma="
                              << (tracer->m_lastOldAutoFullNBeta + 1)
                              << " x " << tracer->m_lastOldAutoFullNGamma
                              << " (div" << tracer->m_lastOldAutoFullDiv
                              << ")" << std::endl;

                    std::vector<double> xSizes;
                    std::vector<std::string> labels;
                    if (args.IsCatched("multikeq") || args.IsCatched("multikeq_list"))
                    {
                        std::vector<double> kSizes;
                        std::string rangeText;
                        if (args.IsCatched("multikeq_list"))
                        {
                            std::string listFile = args.GetStringValue("multikeq_list", 0);
                            kSizes = ReadSizeList(listFile);
                            if (kSizes.empty())
                            {
                                std::cerr << "ERROR: --multikeq_list has no sizes." << std::endl;
                                return 1;
                            }
                            double Kmin = *std::min_element(kSizes.begin(), kSizes.end());
                            double Kmax = *std::max_element(kSizes.begin(), kSizes.end());
                            rangeText = std::to_string(Kmin) + ".."
                                + std::to_string(Kmax) + " from " + listFile;
                        }
                        else
                        {
                            double Kmin = args.GetDoubleValue("multikeq", 0);
                            double Kmax = args.GetDoubleValue("multikeq", 1);
                            int nSizes = args.GetIntValue("multikeq", 2);
                            if (Kmin <= 0 || Kmax <= 0)
                            {
                                std::cerr << "ERROR: --multikeq values must be positive." << std::endl;
                                return 1;
                            }
                            kSizes = GenerateLogSizes(Kmin, Kmax, nSizes);
                            rangeText = std::to_string(Kmin) + ".."
                                + std::to_string(Kmax);
                        }

                        double currentVolumeForMulti = particle->Volume();
                        double currentReqForMulti = (currentVolumeForMulti > 0)
                            ? pow(3.0 * currentVolumeForMulti / (4.0 * M_PI), 1.0 / 3.0)
                            : 0.0;
                        double kRef = (wave > 0 && currentReqForMulti > 0)
                            ? (2.0 * M_PI * currentReqForMulti / wave) : 0.0;
                        double xRef = M_PI * particle->MaximalDimention() / wave;
                        if (kRef <= 0 || xRef <= 0)
                        {
                            std::cerr << "ERROR: --oldautofull multikeq requires positive wavelength and particle volume." << std::endl;
                            return 1;
                        }
                        xSizes.reserve(kSizes.size());
                        labels.reserve(kSizes.size());
                        for (double k : kSizes)
                        {
                            xSizes.push_back(xRef * (k / kRef));
                            labels.push_back("keq" + SizeLabel(k));
                        }
                        std::cout << "Shared oldautofull multikeq: k_eq "
                                  << rangeText << " (" << kSizes.size()
                                  << " sizes), reference k_eq=" << kRef
                                  << std::endl;
                    }
                    else
                    {
                        double Dmin = args.GetDoubleValue("multigrid", 0);
                        double Dmax_mg = args.GetDoubleValue("multigrid", 1);
                        int nSizes = args.GetIntValue("multigrid", 2);
                        if (Dmin <= 0 || Dmax_mg <= 0)
                        {
                            std::cerr << "ERROR: --multigrid values must be positive." << std::endl;
                            return 1;
                        }
                        std::vector<double> dSizes =
                            GenerateLogSizes(Dmin, Dmax_mg, nSizes);
                        double Dref = particle->MaximalDimention();
                        double xRef = M_PI * Dref / wave;
                        xSizes.reserve(dSizes.size());
                        labels.reserve(dSizes.size());
                        for (double d : dSizes)
                        {
                            xSizes.push_back(xRef * (d / Dref));
                            labels.push_back("D" + SizeLabel(d));
                        }
                        std::cout << "Shared oldautofull multigrid: D "
                                  << Dmin << ".." << Dmax_mg << " ("
                                  << dSizes.size()
                                  << " log sizes), reference Dmax=" << Dref
                                  << std::endl;
                    }

                    const bool useSharedOldAutoFullMulti =
                        std::getenv("MBS_OLDAUTOFULL_SHARED_MULTI")
                        && std::atoi(std::getenv(
                            "MBS_OLDAUTOFULL_SHARED_MULTI")) != 0;
                    if (useSharedOldAutoFullMulti)
                    {
                        std::cerr << "WARNING: MBS_OLDAUTOFULL_SHARED_MULTI=1 "
                                  << "uses experimental reference-size beam "
                                  << "cache reuse for --oldautofull "
                                  << "multikeq/multigrid. For non-convex file "
                                  << "particles this can be inaccurate; the "
                                  << "release default is independent retrace."
                                  << std::endl;
                        tracer->TraceRandomMultiSize(betaRange, gammaRange,
                                                     xSizes, labels);
                    }
                    else
                    {
                        const double xRefIndependent =
                            M_PI * particle->MaximalDimention() / wave;
                        tracer->TraceRandomMultiSizeIndependent(
                            betaRange, gammaRange, xSizes, labels,
                            xRefIndependent);
                    }
                }
            }
            else if (isAdaptive)
            {
                if (args.IsCatched("sobol") || args.IsCatched("sobol_seed")
                    || args.IsCatched("sobol_ring")
                    || args.IsCatched("hammersley") || args.IsCatched("lattice")
                    || args.IsCatched("lattice_z")
                    || args.IsCatched("euler_quad"))
                    std::cerr << "WARNING: explicit orientation rule ignored (--auto/--adaptive overrides with adaptive orientations)." << std::endl;
                double epsAdapt = isAuto ? args.GetDoubleValue("auto", 0)
                                         : args.GetDoubleValue("adaptive", 0);
                int maxOrientUser = args.IsCatched("maxorient")
                    ? args.GetIntValue("maxorient", 0) : 0;

                // Adaptive theta grid probe (bisection) before adaptive orient
                if (!args.IsCatched("tgrid") && !args.IsCatched("grid"))
                {
                    double tgridEps = args.IsCatched("auto_tgrid")
                        ? args.GetDoubleValue("auto_tgrid", 0) : 0.05;
                    if (tgridEps <= 0) tgridEps = 0.05;
                    // Probe count from div16 physics estimate
                    double Dmax_p = particle->MaximalDimention();
                    double dd_p = 0.69 * wave / Dmax_p * (180.0 / M_PI) / ringPoints;
                    int nb_p = std::max(1, (int)(RadToDeg(betaSym) / dd_p / 16));
                    int ng_p = std::max(1, (int)(RadToDeg(gammaSym) / dd_p / 16));
                    int nProbeAuto = nb_p * ng_p;
                    int p2a = 1; while (p2a * 2 <= nProbeAuto) p2a *= 2;
                    nProbeAuto = std::max(64, p2a);
                    tracer->TraceAdaptiveTheta(nProbeAuto, betaSym, gammaSym, tgridEps, 8, true);
                    // Grid now set. TraceAdaptive will compute full Mueller.
                }

                tracer->TraceAdaptive(epsAdapt, betaSym, gammaSym, maxOrientUser);
            }
            else
            {
                bool useSobolSeed = args.IsCatched("sobol_seed");
                bool useSobolRing = args.IsCatched("sobol_ring");
                bool useHammersley = args.IsCatched("hammersley");
                bool useLattice = args.IsCatched("lattice") || args.IsCatched("lattice_z");
                bool useEulerQuad = args.IsCatched("euler_quad");
                int nSobol = args.IsCatched("sobol") ? args.GetIntValue("sobol", 0) : 0;

                if (useSobolSeed)
                {
                    if (args.IsCatched("sobol") || args.IsCatched("euler_quad"))
                        std::cerr << "WARNING: --sobol/--euler_quad ignored because --sobol_seed is selected." << std::endl;
                    if (args.IsCatched("multigrid"))
                        std::cerr << "WARNING: --sobol_seed currently runs single-size; --multigrid ignored." << std::endl;
                    int nOrient = args.GetIntValue("sobol_seed", 0);
                    unsigned int seed = (unsigned int)std::max(0, args.GetIntValue("sobol_seed", 1));
                    if (args.IsCatched("auto_tgrid") && !args.IsCatched("grid") && !args.IsCatched("tgrid"))
                    {
                        double tgridEps = args.GetDoubleValue("auto_tgrid", 0);
                        if (tgridEps <= 0) tgridEps = 0.05;
                        tracer->TraceAdaptiveTheta(std::max(64, nOrient), betaSym, gammaSym, tgridEps, 8, true);
                    }
                    tracer->TraceFromSobolSeed(nOrient, seed, betaSym, gammaSym);
                }
                else if (useSobolRing)
                {
                    if (args.IsCatched("sobol") || args.IsCatched("euler_quad"))
                        std::cerr << "WARNING: --sobol/--euler_quad ignored because --sobol_ring is selected." << std::endl;
                    if (args.IsCatched("multigrid"))
                        std::cerr << "WARNING: --sobol_ring currently runs single-size; --multigrid ignored." << std::endl;
                    int nBeta = args.GetIntValue("sobol_ring", 0);
                    int nGamma = args.GetIntValue("sobol_ring", 1);
                    if (args.IsCatched("auto_tgrid") && !args.IsCatched("grid") && !args.IsCatched("tgrid"))
                    {
                        double tgridEps = args.GetDoubleValue("auto_tgrid", 0);
                        if (tgridEps <= 0) tgridEps = 0.05;
                        int nProbe = std::max(64, nBeta * nGamma);
                        tracer->TraceAdaptiveTheta(nProbe, betaSym, gammaSym, tgridEps, 8, true);
                    }
                    tracer->TraceFromSobolRing(nBeta, nGamma, betaSym, gammaSym);
                }
                else if (useHammersley)
                {
                    if (args.IsCatched("sobol") || args.IsCatched("euler_quad")
                        || args.IsCatched("lattice"))
                        std::cerr << "WARNING: other explicit orientation rules ignored because --hammersley is selected." << std::endl;
                    if (args.IsCatched("multigrid"))
                        std::cerr << "WARNING: --hammersley currently runs single-size; --multigrid ignored." << std::endl;
                    int nOrient = args.GetIntValue("hammersley", 0);
                    if (args.IsCatched("auto_tgrid") && !args.IsCatched("grid") && !args.IsCatched("tgrid"))
                    {
                        double tgridEps = args.GetDoubleValue("auto_tgrid", 0);
                        if (tgridEps <= 0) tgridEps = 0.05;
                        tracer->TraceAdaptiveTheta(std::max(64, nOrient), betaSym, gammaSym, tgridEps, 8, true);
                    }
                    tracer->TraceFromHammersley(nOrient, betaSym, gammaSym);
                }
                else if (useLattice)
                {
                    if (args.IsCatched("sobol") || args.IsCatched("euler_quad"))
                        std::cerr << "WARNING: other explicit orientation rules ignored because --lattice is selected." << std::endl;
                    bool explicitGenerator = args.IsCatched("lattice_z");
                    int nOrient = explicitGenerator
                        ? args.GetIntValue("lattice_z", 0)
                        : args.GetIntValue("lattice", 0);
                    int latticeGenerator = explicitGenerator
                        ? args.GetIntValue("lattice_z", 1) : 0;
                    if (args.IsCatched("auto_tgrid") && !args.IsCatched("grid") && !args.IsCatched("tgrid"))
                    {
                        double tgridEps = args.GetDoubleValue("auto_tgrid", 0);
                        if (tgridEps <= 0) tgridEps = 0.05;
                        tracer->TraceAdaptiveTheta(std::max(64, nOrient), betaSym, gammaSym, tgridEps, 8, true);
                    }
                    if (args.IsCatched("multigrid"))
                    {
                        double Dmin = args.GetDoubleValue("multigrid", 0);
                        double Dmax_mg = args.GetDoubleValue("multigrid", 1);
                        int nSizes = args.GetIntValue("multigrid", 2);
                        if (nSizes < 2) nSizes = 2;

                        double D_current = particle->MaximalDimention();
                        double x_ref = M_PI * D_current / wave;

                        std::vector<double> x_sizes;
                        double logMin = log(Dmin), logMax = log(Dmax_mg);
                        for (int i = 0; i < nSizes; ++i) {
                            double D_user = exp(logMin + (logMax - logMin) * i / (nSizes - 1));
                            double x = x_ref * (D_user / Dmax_mg);
                            x_sizes.push_back(x);
                        }
                        x_sizes.back() = x_ref;

                        cout << "Multigrid: " << nSizes << " sizes, D=" << Dmin
                             << ".." << Dmax_mg << " (log), Dmax_particle=" << D_current
                             << ", x_ref=" << x_ref << endl;
                        cout << "  x range: " << x_sizes.front() << " .. " << x_sizes.back() << endl;

                        if (D_current < Dmax_mg * 0.99)
                            std::cerr << "WARNING: -p particle Dmax=" << D_current
                                      << " < multigrid Dmax=" << Dmax_mg
                                      << ". Use -p with largest size." << endl;

                        tracer->TraceLatticeMultiSize(nOrient, betaSym, gammaSym,
                                                      x_sizes, x_ref,
                                                      latticeGenerator);
                    }
                    else
                    {
                        if (explicitGenerator)
                            tracer->TraceFromLatticeGenerator(
                                nOrient, latticeGenerator, betaSym, gammaSym);
                        else
                            tracer->TraceFromLattice(nOrient, betaSym, gammaSym);
                    }
                }
                else if (useEulerQuad)
                {
                    if (args.IsCatched("sobol"))
                        std::cerr << "WARNING: --sobol ignored because --euler_quad is selected." << std::endl;
                    int nBeta = args.GetIntValue("euler_quad", 0);
                    int nGamma = args.GetIntValue("euler_quad", 1);
                    if (args.IsCatched("auto_tgrid") && !args.IsCatched("grid") && !args.IsCatched("tgrid"))
                    {
                        double tgridEps = args.GetDoubleValue("auto_tgrid", 0);
                        if (tgridEps <= 0) tgridEps = 0.05;
                        int nProbe = std::max(64, nBeta * nGamma);
                        tracer->TraceAdaptiveTheta(nProbe, betaSym, gammaSym, tgridEps, 8, true);
                    }
                    if (args.IsCatched("multigrid"))
                    {
                        double Dmin = args.GetDoubleValue("multigrid", 0);
                        double Dmax_mg = args.GetDoubleValue("multigrid", 1);
                        int nSizes = args.GetIntValue("multigrid", 2);
                        if (nSizes < 2) nSizes = 2;

                        double D_current = particle->MaximalDimention();
                        double x_ref = M_PI * D_current / wave;

                        std::vector<double> x_sizes;
                        double logMin = log(Dmin), logMax = log(Dmax_mg);
                        for (int i = 0; i < nSizes; ++i) {
                            double D_user = exp(logMin + (logMax - logMin) * i / (nSizes - 1));
                            double x = x_ref * (D_user / Dmax_mg);
                            x_sizes.push_back(x);
                        }
                        x_sizes.back() = x_ref;

                        cout << "Multigrid: " << nSizes << " sizes, D=" << Dmin
                             << ".." << Dmax_mg << " (log), Dmax_particle=" << D_current
                             << ", x_ref=" << x_ref << endl;
                        cout << "  x range: " << x_sizes.front() << " .. " << x_sizes.back() << endl;

                        if (D_current < Dmax_mg * 0.99)
                            std::cerr << "WARNING: -p particle Dmax=" << D_current
                                      << " < multigrid Dmax=" << Dmax_mg
                                      << ". Use -p with largest size." << endl;

                        tracer->TraceEulerQuadratureMultiSize(nBeta, nGamma,
                                                              betaSym, gammaSym,
                                                              x_sizes, x_ref);
                    }
                    else
                    {
                        tracer->TraceFromEulerQuadrature(nBeta, nGamma, betaSym, gammaSym);
                    }
                }
                else if (args.IsCatched("multigrid"))
                {
                    double Dmin = args.GetDoubleValue("multigrid", 0);
                    double Dmax_mg = args.GetDoubleValue("multigrid", 1);
                    int nSizes = args.GetIntValue("multigrid", 2);
                    if (nSizes < 2) nSizes = 2;

                    // Generate x_sizes in log scale
                    // x_ref = x of the -p particle (reference, largest)
                    // x_sizes are SCALED relative to x_ref:
                    //   x(D) = x_ref * (D / Dmax_particle)
                    // This ensures x_sizes[Dmax] == x_ref exactly.
                    double D_current = particle->MaximalDimention();
                    double x_ref = M_PI * D_current / wave;

                    // Scale Dmin/Dmax to MaximalDimention scale
                    // User gives "characteristic D" but particle Dmax may differ
                    // Last size MUST equal x_ref for correctness
                    std::vector<double> x_sizes;
                    double logMin = log(Dmin), logMax = log(Dmax_mg);
                    for (int i = 0; i < nSizes; ++i) {
                        double D_user = exp(logMin + (logMax - logMin) * i / (nSizes - 1));
                        double x = x_ref * (D_user / Dmax_mg);
                        x_sizes.push_back(x);
                    }
                    // Force last = x_ref
                    x_sizes.back() = x_ref;

                    cout << "Multigrid: " << nSizes << " sizes, D=" << Dmin
                         << ".." << Dmax_mg << " (log), Dmax_particle=" << D_current
                         << ", x_ref=" << x_ref << endl;
                    cout << "  x range: " << x_sizes.front() << " .. " << x_sizes.back() << endl;

                    if (D_current < Dmax_mg * 0.99)
                        std::cerr << "WARNING: -p particle Dmax=" << D_current
                                  << " < multigrid Dmax=" << Dmax_mg
                                  << ". Use -p with largest size." << endl;

                    tracer->TraceSobolMultiSize(nSobol, betaSym, gammaSym, x_sizes, x_ref);
                }
                else if (args.IsCatched("auto_tgrid") && !args.IsCatched("grid") && !args.IsCatched("tgrid"))
                {
                    double tgridEps = args.GetDoubleValue("auto_tgrid", 0);
                    if (tgridEps <= 0) tgridEps = 0.05;
                    tracer->TraceAdaptiveTheta(nSobol, betaSym, gammaSym, tgridEps, 8);
                }
                else
                {
                    tracer->TraceFromSobol(nSobol, betaSym, gammaSym);
                }
            }

            delete handler;
        }
        else
        {
            cout << endl << "error";
        }
    }
    else // go
    {
        additionalSummary += "Geometrical optics";

        TracerGO tracer(particle, reflNum, dirName);
        tracer.m_scattering->m_wave = wave;
        ApplyTraceCutoffOptions(args, tracer.m_scattering);
        if (args.IsCatched("r"))
        {
            tracer.m_scattering->restriction = args.GetDoubleValue("r", 0);
        }
        tracer.m_summary = additionalSummary;
        tracer.SetIsOutputGroups(isOutputGroups);

        if (args.IsCatched("log"))
        {
            tracer.m_logTime = args.GetIntValue("log");
        }

        HandlerGO *handler;

        if (args.IsCatched("tr"))
        {
            handler = new HandlerTracksGO(particle, &tracer.m_incidentLight, nTheta, wave);
            handler->SetTracks(&trackGroups);
        }
        else
        {
            handler = new HandlerTotalGO(particle, &tracer.m_incidentLight, nTheta, wave);
        }

        tracer.m_summary = additionalSummary;
        ScatteringRange grid = SetConus(args);
        handler->isCoh = !args.IsCatched("incoh");
            handler->m_legacySign = args.IsCatched("legacy_sign");
        ApplyNphiOverride(args, grid);
        handler->SetScatteringSphere(grid);
        handler->SetAbsorptionAccounting(isAbs);
        ApplyAbsorptionPointOption(args, handler);
        ApplyBeamCutoffOptions(args, handler);
        tracer.SetHandler(handler);

        if (args.IsCatched("fixed"))
        {
            additionalSummary += ", fixed orientation, ";
            double beta  = args.GetDoubleValue("fixed", 0);
            double gamma = args.GetDoubleValue("fixed", 1);
            additionalSummary += "zenith " + to_string(beta) + "°, azimuth " + to_string(gamma) + "°" + "\n\n";
            cout << additionalSummary;
            tracer.m_summary = additionalSummary;
            tracer.TraceFixed(beta, gamma);
        }
        else if (args.IsCatched("random"))
        {
            additionalSummary += ", random orientation, ";
            AngleRange beta = GetRange(args, "b", particle);
            AngleRange gamma = GetRange(args, "g", particle);
            MakeMirrorGammaRangeConsistent(args, particle, gamma);
            additionalSummary += "grid: " + to_string(beta.number) + "x" + to_string(gamma.number) + "\n\n";
            cout << additionalSummary;
            tracer.m_summary = additionalSummary;
            tracer.TraceRandom(beta, gamma);
        }
        else if (args.IsCatched("montecarlo"))
        {
            additionalSummary += ", Monte Carlo method, ";
            double betaMax = particle->GetSymmetry().beta;
            double gammaMax = particle->GetSymmetry().gamma;
            if (args.IsCatched("mirror_gamma"))
                gammaMax *= 0.5;
            int nOr = args.GetIntValue("montecarlo", 0);
            if (args.IsCatched("mirror_gamma") && (nOr % 2 != 0))
            {
                std::cerr << "WARNING: --mirror_gamma rounded montecarlo "
                          << "gamma grid count to an even half-domain count: "
                          << nOr << " -> " << (nOr + 1) << std::endl;
                ++nOr;
            }
            AngleRange beta(0, betaMax, nOr);
            AngleRange gamma(0, gammaMax, nOr);
            cout << additionalSummary;
            tracer.m_summary = additionalSummary;
            tracer.TraceMonteCarlo(beta, gamma, nOr);
        }

        delete handler;
    }

    if (mpi_rank == 0)
    {
        if (!args.IsCatched("po") && !args.IsCatched("fixed")
            && !args.IsCatched("random") && !args.IsCatched("montecarlo")
            && !args.IsCatched("oldauto") && !args.IsCatched("sobol")
            && !args.IsCatched("sobol_seed")
            && !args.IsCatched("sobol_ring")
            && !args.IsCatched("hammersley") && !args.IsCatched("lattice")
            && !args.IsCatched("lattice_z")
            && !args.IsCatched("euler_quad") && !args.IsCatched("adaptive")
            && !args.IsCatched("auto") && !args.IsCatched("autofull")
            && !args.IsCatched("oldautofull"))
        {
            cerr << "\nWARNING: No computation mode specified. Nothing was computed.\n"
                 << "  Use --po for Physical Optics, or specify orientation mode:\n"
                 << "    --fixed BETA GAMMA    Fixed orientation\n"
                 << "    --random B G          Random orientation grid\n"
                 << "    --sobol N             Quasi-random orientations\n"
                 << "    --sobol_seed N S      Sobol with nested Owen scramble seed\n"
                 << "    --sobol_ring B G      Sobol beta x uniform gamma ring\n"
                 << "    --hammersley N        Hammersley orientations\n"
                 << "    --lattice N           Rank-1 lattice orientations\n"
                 << "    --lattice_z N Z       Rank-1 lattice with generator Z\n"
                 << "    --euler_quad B G      High-order Euler quadrature\n"
                 << "    --oldauto DIV         Physics-based auto grid\n"
                 << "    --auto EPS            Full auto mode\n"
                 << "  Example: mbs_po --po -p 1 10 5 --ri 1.3 0 -w 1 -n 4 --oldauto 8\n";
        }
        cout << endl << "done";
    }

    if (!args.IsCatched("close") && mpi_rank == 0)
    {
        getchar();
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
