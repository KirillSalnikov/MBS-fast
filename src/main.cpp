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
    parser.AddRule("sizefile", 1, true); // multi-size: file with size parameters (one per line)
    parser.AddRule("tgrid", 1, true); // non-uniform theta grid file
    parser.AddRule("beam_cutoff", 1, true); // beam importance cutoff (relative to C_geo)
    parser.AddRule("sobol", 1, true); // Sobol quasi-random orientations (number, power of 2)
    parser.AddRule("auto_tgrid", 0, true); // auto-generate theta grid based on size parameter
    parser.AddRule("auto_phi", 0, true); // auto-select N_phi based on size parameter
    parser.AddRule("adaptive", 1, true); // adaptive convergence (target relative accuracy)
    parser.AddRule("autofull", 1, true); // full 3D sequential: n → N_phi → N_orient
    parser.AddRule("auto", 1, true); // full auto: auto_tgrid + auto_phi + adaptive (one arg: eps)
    parser.AddRule("maxorient", 1, true); // max orientations for adaptive (power of 2)
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
    else
    {
        std::cerr << "ERROR!!!!!!! Wrong \"grid\" argument number";
        exit(1);
    }

    // Apply non-uniform theta grid if --tgrid is specified
    if (parser.IsCatched("tgrid"))
    {
        std::string tgridFile = parser.GetStringValue("tgrid", 0);
        if (!range.LoadThetaGrid(tgridFile))
        {
            std::cerr << "ERROR: Cannot load theta grid from file: " << tgridFile << std::endl;
            exit(1);
        }
        std::cout << "Non-uniform theta grid: " << (range.nZenith + 1)
                  << " points from " << RadToDeg(range.zenithStart)
                  << " to " << RadToDeg(range.zenithEnd) << " deg" << std::endl;
    }

    return range;
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
        else if (args.IsCatched("random"))
        {
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
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) tpt->SetMPI(mpi_rank, mpi_size); }
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
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) tpt->SetMPI(mpi_rank, mpi_size); }
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
                tracer->TraceRandom(beta, gamma);
            }

            delete handler;
        }
        else if (args.IsCatched("montecarlo"))
        {
            additionalSummary += ", Monte Carlo method\n\n";
            AngleRange beta = GetRange(args, "b", particle);
            AngleRange gamma = GetRange(args, "g", particle);

            HandlerPO *handler;

            TracerPOTotal *tracer;
            ScatteringRange conus = SetConus(args);

            tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) tpt->SetMPI(mpi_rank, mpi_size); }
            trackGroups.push_back(TrackGroup());
            handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                         nTheta, wave);

            cout << additionalSummary;
            tracer->m_summary = additionalSummary;

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

            int nOr = args.GetIntValue("montecarlo", 0);
            tracer->TraceMonteCarlo(beta, gamma, nOr);
        }
        else if (args.IsCatched("orientfile"))
        {
            additionalSummary += ", orientations from file\n\n";

            TracerPOTotal *tracer;
            ScatteringRange conus = SetConus(args);

            tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) tpt->SetMPI(mpi_rank, mpi_size); }
            tracer->m_scattering->m_wave = wave;
            tracer->shadowOff = args.IsCatched("shadow_off");
            trackGroups.push_back(TrackGroup());
            HandlerPOTotal *handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                         nTheta, wave);

            cout << additionalSummary;
            tracer->m_summary = additionalSummary;

            handler->isCoh = !args.IsCatched("incoh");
            handler->useKarczewski = args.IsCatched("karczewski");
            handler->SetScatteringSphere(conus);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);

            // Beam importance cutoff
            if (args.IsCatched("beam_cutoff"))
            {
                double eps = args.GetDoubleValue("beam_cutoff", 0);
                double D_particle = particle->MaximalDimention();
                double C_geo = M_PI * D_particle * D_particle / 4.0;
                handler->SetBeamCutoffRelative(eps, C_geo);
                cout << "Beam cutoff: eps=" << eps << ", C_geo=" << C_geo
                     << ", threshold=" << handler->m_beamCutoff << endl;
            }

            tracer->SetIsOutputGroups(isOutputGroups);
            tracer->SetHandler(handler);

            std::string orientFileName = args.GetStringValue("orientfile", 0);

            if (args.IsCatched("sizefile"))
            {
                std::string sizeFileName = args.GetStringValue("sizefile", 0);
                std::vector<double> x_sizes;
                std::ifstream sizeFile(sizeFileName);
                double xval;
                while (sizeFile >> xval) x_sizes.push_back(xval);
                sizeFile.close();

                // x_ref is determined from current particle size and wavelength
                // x = pi * D / lambda
                double D_current = particle->MaximalDimention();
                double x_ref = M_PI * D_current / wave;

                cout << "Multi-size mode: x_ref=" << x_ref
                     << ", computing " << x_sizes.size() << " sizes:";
                for (double xs : x_sizes) cout << " " << xs;
                cout << endl;

                tracer->TraceFromFileMultiSize(orientFileName, x_sizes, x_ref);
            }
            else
            {
                tracer->TraceFromFile(orientFileName);
            }

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
            // --auto/--autofull: --grid is optional (default 0→180°)
            ScatteringRange conus = args.IsCatched("grid")
                ? SetConus(args)
                : ScatteringRange(0, M_PI, 1, 1);

            tracer = new TracerPOTotal(particle, reflNum, dirName);
            { TracerPOTotal *tpt = dynamic_cast<TracerPOTotal*>(tracer); if(tpt) tpt->SetMPI(mpi_rank, mpi_size); }
            tracer->m_scattering->m_wave = wave;
            tracer->shadowOff = args.IsCatched("shadow_off");
            trackGroups.push_back(TrackGroup());
            HandlerPOTotal *handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                         nTheta, wave);

            // Apply auto_tgrid if requested (or implied by --auto)
            if (args.IsCatched("auto_tgrid") || isAuto)
            {
                double D = particle->MaximalDimention();
                ApplyAutoThetaGrid(conus, D, wave);
            }

            // Apply auto_phi if requested (or implied by --auto)
            if (args.IsCatched("auto_phi") || isAuto)
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
            handler->useKarczewski = args.IsCatched("karczewski");
            handler->SetScatteringSphere(conus);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);

            // Beam importance cutoff
            if (args.IsCatched("beam_cutoff"))
            {
                double eps_bc = args.GetDoubleValue("beam_cutoff", 0);
                double D_particle = particle->MaximalDimention();
                double C_geo = M_PI * D_particle * D_particle / 4.0;
                handler->SetBeamCutoffRelative(eps_bc, C_geo);
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
                double epsAdapt = isAuto ? args.GetDoubleValue("auto", 0)
                                         : args.GetDoubleValue("adaptive", 0);
                int maxOrientUser = args.IsCatched("maxorient")
                    ? args.GetIntValue("maxorient", 0) : 0;
                tracer->TraceAdaptive(epsAdapt, betaSym, gammaSym, maxOrientUser);
            }
            else
            {
                int nSobol = args.GetIntValue("sobol", 0);
                tracer->TraceFromSobol(nSobol, betaSym, gammaSym);
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
            AngleRange beta = GetRange(args, "b", particle);
            AngleRange gamma = GetRange(args, "g", particle);
            cout << additionalSummary;
            tracer.m_summary = additionalSummary;
            int nOr = args.GetIntValue("montecarlo", 0);
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
