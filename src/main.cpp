#include <iostream>
#include <fstream>
#include <assert.h>
#include <float.h>
#include <string>

#ifdef USE_MPI
#include <mpi.h>
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

#ifdef _OUTPUT_NRG_CONV
ofstream energyFile("energy.dat", ios::out);
double SSconfined=0;
int bcount=0;
#endif

using namespace std;
using namespace chrono;

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
    parser.AddRule("p", '+'); // particle (type, size, ...)
    parser.AddRule("ri", 2); // refractive index (Re and Im parts)
    parser.AddRule("n", 1, true); // number of internal reflection (optional for --autofull)
    parser.AddRule("pf", zero, true); // particle (filename)
    parser.AddRule("rs", 1, true, "pf"); // resize particle (new size)
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
    parser.AddRule("close", 0, true); // closing of program after calculation
    parser.AddRule("o", 1, true); // output folder name
    parser.AddRule("gr", zero, true);
    parser.AddRule("filter", 1, true); // scattering angle filter
    parser.AddRule("shadow", zero, true);
    parser.AddRule("incoh", zero, true);
    parser.AddRule("jones", zero, true);
    parser.AddRule("shadow_off", zero, true);
    parser.AddRule("forced_nonconvex", zero, true);
    parser.AddRule("forced_convex", zero, true);
    parser.AddRule("r", 1, true); // restriction ratio for small beams when intersection (100 by default)
    parser.AddRule("log", 1, true); // time of writing progress (in seconds)
    parser.AddRule("multigrid", 1, true); // multi-size: file with -p params (same type, diff size), one per line
    parser.AddRule("save_betas", 0, true); // save intermediate Mueller for each beta to betas/ subfolder
    parser.AddRule("tgrid", 1, true); // non-uniform theta grid file
    parser.AddRule("beam_cutoff", 1, true); // beam importance cutoff (relative to C_geo)
    parser.AddRule("sobol", 1, true); // Sobol quasi-random orientations (number, power of 2)
    parser.AddRule("auto_tgrid", 1, true); // adaptive theta grid (arg: tolerance, e.g. 0.05)
    parser.AddRule("auto_phi", 0, true); // auto-select N_phi based on size parameter
    parser.AddRule("adaptive", 1, true); // adaptive convergence (target relative accuracy)
    parser.AddRule("autofull", 1, true); // full 3D sequential: n → N_phi → N_orient
    parser.AddRule("auto", 1, true); // full auto: auto_tgrid + auto_phi + adaptive (one arg: eps)
    parser.AddRule("maxorient", 1, true); // max orientations for adaptive (power of 2)
    parser.AddRule("oldauto", 1, true); // physics-based: div2/div4/div8 of diffraction-limited grid
    parser.AddRule("coh_orient", 0, true); // coherent across orientations (legacy mode)
    parser.AddRule("legacy_sign", 0, true); // use old (+) Fresnel sign for forward direction
    parser.AddRule("sym", 2, true); // symmetry override: beta_factor gamma_factor (e.g. --sym 2 6)
}

ScatteringRange SetConus(ArgPP &parser)
{
    ScatteringRange range(0, M_PI, 1, 1); // placeholder

    if (parser.GetArgNumber("grid") == 3)
    {
        double radius = parser.GetDoubleValue("grid", 0);
        int nAz = parser.GetDoubleValue("grid", 1);
        int nZen = parser.GetDoubleValue("grid", 2);
        range = ScatteringRange(M_PI - DegToRad(radius), M_PI, nAz, nZen);
    }
    else if (parser.GetArgNumber("grid") == 4)
    {
        double zenStart = parser.GetDoubleValue("grid", 0);
        double zenEnd = parser.GetDoubleValue("grid", 1);

        if (zenStart < zenEnd)
        {
            int nAz = parser.GetDoubleValue("grid", 2);
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
            range.nAzimuth = 48;
            range.azinuthStep = 2.0 * M_PI / 48;
        }
        std::cout << "Non-uniform theta grid: " << (range.nZenith + 1)
                  << " points from " << RadToDeg(range.zenithStart)
                  << " to " << RadToDeg(range.zenithEnd) << " deg" << std::endl;
    }

    return range;
}

/// Parse --multigrid file: each line = "type L D [extra]" (same format as -p args).
/// Returns x_sizes (size parameters) and stores the reference particle index (largest).
/// All particles must have the same type.
std::vector<double> ParseMultigridFile(const std::string &filename, double wave,
                                        const complex &refrIndex, int &refIndex)
{
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "Error: cannot open multigrid file: " << filename << std::endl;
        throw std::exception();
    }

    struct PEntry { int type; double height, diameter, extra; int nArgs; };
    std::vector<PEntry> entries;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        PEntry e; e.extra = 0; e.nArgs = 0;
        if (!(iss >> e.type >> e.height >> e.diameter)) continue;
        e.nArgs = 3;
        if (iss >> e.extra) e.nArgs = 4;
        entries.push_back(e);
    }
    f.close();

    if (entries.empty()) {
        std::cerr << "Error: no entries in multigrid file" << std::endl;
        throw std::exception();
    }

    // All must be same type
    for (size_t i = 1; i < entries.size(); ++i) {
        if (entries[i].type != entries[0].type) {
            std::cerr << "Error: multigrid requires same particle type on all lines" << std::endl;
            throw std::exception();
        }
    }

    // Compute x = pi * Dmax / lambda for each entry using a temporary particle
    std::vector<double> x_sizes;
    double maxDmax = 0;
    refIndex = 0;
    for (size_t i = 0; i < entries.size(); ++i) {
        auto &e = entries[i];
        Particle *tmp = nullptr;
        switch ((ParticleType)e.type) {
        case ParticleType::Hexagonal:
            tmp = new Hexagonal(refrIndex, e.diameter, e.height); break;
        case ParticleType::ConcaveHexagonal: {
            double sup = (e.nArgs >= 4) ? e.extra : 0;
            tmp = new ConcaveHexagonal(refrIndex, e.diameter, e.height, sup); break;
        }
        case ParticleType::Bullet: {
            double sup = e.diameter * sqrt(3) * tan(DegToRad(62)) / 4;
            tmp = new Bullet(refrIndex, e.diameter, e.height, sup); break;
        }
        default:
            tmp = new Hexagonal(refrIndex, e.diameter, e.height); break;
        }
        double Dmax = tmp->MaximalDimention();
        double x = M_PI * Dmax / wave;
        x_sizes.push_back(x);
        if (Dmax > maxDmax) { maxDmax = Dmax; refIndex = i; }
        delete tmp;
    }

    return x_sizes;
}

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

    return AngleRange(min, max, number + 1);
}

int main(int argc, const char* argv[])
{
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

    if (argc <= 1) // no arguments
    {
        cout << "No arguments. Press any key to exit..." << endl;
        getchar();
        return 1;
    }

    ArgPP args;
    SetArgRules(args);
    args.Parse(argc, argv);

    bool isAbs = args.IsCatched("abs");

    double re = args.GetDoubleValue("ri", 0);
    double im = args.GetDoubleValue("ri", 1);
    complex refrIndex = complex(re, im);

    if (!isAbs)
    {
        refrIndex = complex(re, 0);
    }

    // TODO: AggregateBuilder

    Particle *particle = nullptr;


    additionalSummary += "Particle: ";

    if (args.IsCatched("pf"))
    {
        std::string filename = args.GetStringValue("p");
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

    additionalSummary += "\tRefractive index: " + to_string(re);

    if (fabs(im) > FLT_EPSILON)
    {
        additionalSummary +=  " + i\n";
    }

    additionalSummary += "\n";

    particle->Output("particle_for_check.dat");
    additionalSummary += "\tArea:" + to_string(particle->Area()) + "\n\n";

    int reflNum = args.IsCatched("n") ? (int)args.GetDoubleValue("n") : 6; // default n=6
    additionalSummary += "Number of secondary reflections: " + to_string(reflNum) + "\n";

    string dirName = (args.IsCatched("o")) ? args.GetStringValue("o")
                                           : "M";
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

    bool isOutputGroups = args.IsCatched("gr");
    double wave = args.IsCatched("w") ? args.GetDoubleValue("w") : 0;
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

    int nTheta = args.IsCatched("grid") ? (int)args.GetDoubleValue("grid", 2) : 1;

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

            HandlerPO *handler = new HandlerPO(particle, &tracer.m_incidentLight,
                                               nTheta, wave);
            if (args.IsCatched("r"))
            {
                tracer.m_scattering->restriction = args.GetDoubleValue("r", 0);
            }

            handler->isCoh = !args.IsCatched("incoh");
            handler->m_legacySign = args.IsCatched("legacy_sign");
            handler->useKarczewski = args.IsCatched("karczewski");
            handler->outputJones = args.IsCatched("jones");
            handler->SetScatteringSphere(bsCone);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);

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
            int points_per_ring = 3; // from Excel AB4

            // Angular step for orientation grid
            double orient_step = delta_theta_deg / points_per_ring;

            // Get symmetry from particle
            double betaSym_deg = RadToDeg(particle->GetSymmetry().beta);
            double gammaSym_deg = RadToDeg(particle->GetSymmetry().gamma);

            // Full grid: N = sym_range / step
            int N_beta_full = (int)ceil(betaSym_deg / orient_step);
            int N_gamma_full = (int)ceil(gammaSym_deg / orient_step);

            // Divide by div factor
            int N_beta = N_beta_full / div;
            int N_gamma = N_gamma_full / div;
            // Make odd (from Excel ODD function)
            if (N_beta % 2 == 0) N_beta++;
            if (N_gamma % 2 == 0) N_gamma++;
            if (N_beta < 3) N_beta = 3;
            if (N_gamma < 3) N_gamma = 3;

            // N_phi default 360, can override with --grid
            int N_phi = args.IsCatched("grid") ? (int)args.GetDoubleValue("grid", 2) : 360;

            cout << "=== oldauto (physics-based) ===" << endl;
            cout << "  L=" << L << " um, lambda=" << wave << " um" << endl;
            cout << "  Diffraction step: " << delta_theta_deg << " deg" << endl;
            cout << "  Orient step: " << orient_step << " deg (3 pts/ring)" << endl;
            cout << "  Full grid: " << N_beta_full << " x " << N_gamma_full
                 << " = " << (long long)N_beta_full * N_gamma_full << endl;
            cout << "  div" << div << ": " << N_beta << " x " << N_gamma
                 << " = " << N_beta * N_gamma << " orientations" << endl;
            cout << "  N_phi=" << N_phi << endl;
            cout << "  n=" << reflNum << endl;

            additionalSummary += ", oldauto div" + to_string(div) + "\n\n";

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
                conus.nAzimuth = N_phi;
                conus.azinuthStep = 2.0 * M_PI / N_phi;
            }

            TracerPOTotal *tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); } }
            tracer->m_scattering->m_wave = wave;
            tracer->shadowOff = args.IsCatched("shadow_off");
            trackGroups.push_back(TrackGroup());
            HandlerPOTotal *handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                         nTheta, wave);

            if (args.IsCatched("beam_cutoff"))
                handler->m_targetEps = args.GetDoubleValue("beam_cutoff", 0);

            cout << additionalSummary;
            tracer->m_summary = additionalSummary;

            handler->isCoh = !args.IsCatched("incoh");
            handler->m_legacySign = args.IsCatched("legacy_sign");
            handler->useKarczewski = args.IsCatched("karczewski");
            handler->SetScatteringSphere(conus);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);

            tracer->SetIsOutputGroups(isOutputGroups);
            tracer->SetHandler(handler);

            // Use --random with computed N_beta, N_gamma
            AngleRange betaRange(0, particle->GetSymmetry().beta, N_beta);
            AngleRange gammaRange(0, particle->GetSymmetry().gamma, N_gamma);

            tracer->TraceRandom(betaRange, gammaRange);
            delete handler;
        }
        else if (args.IsCatched("random"))
        {
            // Warn if --auto/--adaptive also specified (random takes priority)
            if (args.IsCatched("auto") || args.IsCatched("autofull") || args.IsCatched("adaptive"))
                std::cerr << "WARNING: --random overrides --auto/--adaptive. Use --sobol instead." << std::endl;
            additionalSummary += ", random orientation\n\n";
            AngleRange beta = GetRange(args, "b", particle);
            AngleRange gamma = GetRange(args, "g", particle);

            HandlerPO *handler;

            if (args.IsCatched("point"))
            {
//                 TracerBackScatterPoint tracer(particle, reflNum, dirName);
//                 tracer.m_scattering->m_wave = wave;
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
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); } }
                    tracer->m_scattering->m_wave = wave;
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
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); } }
                    tracer->m_scattering->m_wave = wave;
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
            handler->m_legacySign = args.IsCatched("legacy_sign");
                handler->useKarczewski = args.IsCatched("karczewski");
                handler->SetScatteringSphere(conus);
                handler->SetTracks(&trackGroups);
                handler->SetAbsorptionAccounting(isAbs);

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
                        int nProbe = std::min(256, (beta.number+1) * gamma.number);
                        double betaSym2 = particle->GetSymmetry().beta;
                        double gammaSym2 = particle->GetSymmetry().gamma;
                        tpt->TraceAdaptiveTheta(nProbe, betaSym2, gammaSym2, tgridEps, 8);
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
            int nOr = args.GetIntValue("montecarlo", 0);
            AngleRange beta(0, betaMax, nOr);
            AngleRange gamma(0, gammaMax, nOr);

            HandlerPO *handler;

            TracerPOTotal *tracer;
            ScatteringRange conus = SetConus(args);

            tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); } }
            tracer->m_scattering->m_wave = wave;
            tracer->shadowOff = args.IsCatched("shadow_off");
            if (args.IsCatched("r"))
                tracer->m_scattering->restriction = args.GetDoubleValue("r", 0);
            trackGroups.push_back(TrackGroup());
            handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                         nTheta, wave);

            cout << additionalSummary;
            tracer->m_summary = additionalSummary;

            handler->isCoh = !args.IsCatched("incoh");
            handler->m_legacySign = args.IsCatched("legacy_sign");
            handler->SetScatteringSphere(conus);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);

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
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); } }
            tracer->m_scattering->m_wave = wave;
            tracer->shadowOff = args.IsCatched("shadow_off");
            trackGroups.push_back(TrackGroup());
            HandlerPOTotal *handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                         nTheta, wave);

            cout << additionalSummary;
            tracer->m_summary = additionalSummary;

            handler->isCoh = !args.IsCatched("incoh");
            handler->m_legacySign = args.IsCatched("legacy_sign");
            handler->useKarczewski = args.IsCatched("karczewski");
            handler->SetScatteringSphere(conus);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);

            // Beam cutoff
            if (args.IsCatched("beam_cutoff"))
            {
                handler->m_targetEps = args.GetDoubleValue("beam_cutoff", 0);
            }

            tracer->SetIsOutputGroups(isOutputGroups);
            tracer->SetHandler(handler);

            std::string orientFileName = args.GetStringValue("orientfile", 0);

            tracer->TraceFromFile(orientFileName);

            delete handler;
        }
        else if (args.IsCatched("sobol") || args.IsCatched("adaptive") || args.IsCatched("auto") || args.IsCatched("autofull"))
        {
            // --auto EPS implies --adaptive EPS + --auto_tgrid + --auto_phi
            bool isAuto = args.IsCatched("auto") || args.IsCatched("autofull");
            bool isAutoFull = args.IsCatched("autofull");
            bool isAdaptive = args.IsCatched("adaptive") || isAuto;
            if (isAdaptive)
                additionalSummary += ", adaptive Sobol\n\n";
            else
                additionalSummary += ", Sobol quasi-random\n\n";

            TracerPOTotal *tracer;
            // --grid optional for --auto/--autofull/--tgrid
            ScatteringRange conus = (args.IsCatched("grid") || args.IsCatched("tgrid"))
                ? SetConus(args)
                : ScatteringRange(0, M_PI, 1, 1);

            tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) { tpt->SetMPI(mpi_rank, mpi_size); tpt->m_cohOrient = args.IsCatched("coh_orient"); tpt->m_saveBetas = args.IsCatched("save_betas"); } }
            tracer->m_scattering->m_wave = wave;
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
            handler->m_legacySign = args.IsCatched("legacy_sign");
            handler->useKarczewski = args.IsCatched("karczewski");
            handler->SetScatteringSphere(conus);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);

            // Beam importance cutoff
            // Beam cutoff: --beam_cutoff EPS sets relative accuracy for beam skipping
            // cutoff = EPS³ × totalBeamEnergy (computed per orientation in PrepareBeams)
            if (args.IsCatched("beam_cutoff"))
            {
                handler->m_targetEps = args.GetDoubleValue("beam_cutoff", 0);
            }

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

            // Pass target accuracy to handler for beam cutoff
            // --beam_cutoff has priority over --auto/--adaptive eps
            if (!args.IsCatched("beam_cutoff") && isAdaptive) {
                double epsForCutoff = args.IsCatched("autofull")
                    ? args.GetDoubleValue("autofull", 0)
                    : (args.IsCatched("auto") ? args.GetDoubleValue("auto", 0)
                    : (args.IsCatched("adaptive") ? args.GetDoubleValue("adaptive", 0) : 0.01));
                handler->m_targetEps = epsForCutoff;
            }

            if (isAutoFull)
            {
                double epsAdapt = args.GetDoubleValue("autofull", 0);
                int maxOrientUser = args.IsCatched("maxorient")
                    ? args.GetIntValue("maxorient", 0) : 0;
                tracer->TraceAutoFull(epsAdapt, betaSym, gammaSym, maxOrientUser,
                                      particle, wave, conus, handler);
            }
            else if (isAdaptive)
            {
                if (args.IsCatched("sobol"))
                    std::cerr << "WARNING: --sobol N ignored (--auto/--adaptive overrides with adaptive orientations)." << std::endl;
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
                    tracer->TraceAdaptiveTheta(256, betaSym, gammaSym, tgridEps, 8, true);
                    // Grid now set. TraceAdaptive will compute full Mueller.
                }

                tracer->TraceAdaptive(epsAdapt, betaSym, gammaSym, maxOrientUser);
            }
            else
            {
                int nSobol = args.GetIntValue("sobol", 0);

                if (args.IsCatched("multigrid"))
                {
                    int refIdx;
                    std::vector<double> x_sizes = ParseMultigridFile(
                        args.GetStringValue("multigrid", 0), wave, refrIndex, refIdx);
                    double x_ref = x_sizes[refIdx];
                    // Current particle is already the reference (largest from -p)
                    double D_current = particle->MaximalDimention();
                    double x_current = M_PI * D_current / wave;

                    cout << "Multigrid: " << x_sizes.size() << " sizes, x_ref=" << x_ref
                         << " (current x=" << x_current << ")" << endl;

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
        handler->SetScatteringSphere(grid);
        handler->SetAbsorptionAccounting(isAbs);
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
            int nOr = args.GetIntValue("montecarlo", 0);
            AngleRange beta(0, betaMax, nOr);
            AngleRange gamma(0, gammaMax, nOr);
            cout << additionalSummary;
            tracer.m_summary = additionalSummary;
            tracer.TraceMonteCarlo(beta, gamma, nOr);
        }

        delete handler;
    }

    if (mpi_rank == 0)
        cout << endl << "done";

    if (!args.IsCatched("close") && mpi_rank == 0)
    {
        getchar();
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
