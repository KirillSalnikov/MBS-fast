#include "CliOptions.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

#include "ArgPP.h"

namespace
{

std::string FlagToken(const std::string &name)
{
    return (name.size() == 1 ? "-" : "--") + name;
}

void PrintWrapped(std::ostream &out, const std::string &prefix,
                  const std::string &text, size_t width = 100)
{
    const std::string continuation(prefix.size(), ' ');
    std::istringstream words(text);
    std::string word;
    std::string line = prefix;
    bool haveWord = false;
    while (words >> word)
    {
        if (haveWord && line.size() + 1 + word.size() > width)
        {
            out << line << '\n';
            line = continuation + word;
        }
        else
        {
            if (haveWord)
                line += ' ';
            line += word;
        }
        haveWord = true;
    }
    out << line << '\n';
}

std::string AliasText(const CliOptionSpec &spec)
{
    std::vector<std::string> names;
    if (spec.key != spec.canonical)
        names.push_back(FlagToken(spec.key));
    for (const std::string &alias : spec.aliases)
    {
        if (alias != spec.canonical && alias != spec.key)
            names.push_back(FlagToken(alias));
    }
    if (names.empty())
        return std::string();

    std::ostringstream out;
    out << " Compatibility: ";
    for (size_t i = 0; i < names.size(); ++i)
    {
        if (i)
            out << ", ";
        out << names[i];
    }
    out << '.';
    return out.str();
}

} // namespace

CliOptionSpec::CliOptionSpec(const char *keyValue,
                             const char *canonicalValue,
                             char valueCountValue,
                             const char *valueSyntaxValue,
                             const char *sectionValue,
                             const char *descriptionValue,
                             bool debugOnlyValue,
                             std::initializer_list<const char *> extraAliases)
    : key(keyValue),
      canonical(canonicalValue),
      valueCount(valueCountValue),
      valueSyntax(valueSyntaxValue),
      section(sectionValue),
      description(descriptionValue),
      debugOnly(debugOnlyValue)
{
    for (const char *alias : extraAliases)
        aliases.push_back(alias);
}

const std::vector<CliOptionSpec> &GetCliOptionSpecs()
{
    static const std::vector<CliOptionSpec> options = {
        // Core method selectors. The zero-argument switches remain for existing scripts.
        {"method", "method", 1, "MODE", "Method",
         "Calculation method: po or go."},
        {"po", "po", 0, "", "Method",
         "Compatibility switch equivalent to --method po."},
        {"go", "go", 0, "", "Method",
         "Compatibility switch equivalent to --method go."},
        {"incoh", "incoherent", 0, "", "Method",
         "Accumulate a Mueller matrix per beam instead of summing Jones amplitudes coherently."},
        {"karczewski", "karczewski", 0, "", "Method",
         "Use the experimental Karczewski polarization rotation.", true},
        {"legacy_sign", "legacy-sign", 0, "", "Method",
         "Use the legacy positive Fresnel sign convention.", true},

        // Particle and optical parameters.
        {"p", "particle", '+', "TYPE PARAMETERS...", "Particle and physics",
         "Create a built-in particle: 1 L D hexagonal; 2 L D bullet; 3 L D [CAP] rosette; 4 SCALE droxtal; 10 L D CAVITY_DEG concave hexagonal; 12 L D 2 two-column aggregate; 999 SCALE fixed aggregate."},
        {"pf", "particle-file", 1, "FILE", "Particle and physics",
         "Load particle geometry from FILE."},
        {"rs", "resize-dmax-um", 1, "DMAX", "Particle and physics",
         "Resize a file particle to Dmax in micrometers; requires --particle-file."},
        {"k_eq", "k-eq", 1, "K", "Particle and physics",
         "Scale the particle to equivalent-volume size parameter K = 2*pi*r_eq/lambda."},
        {"ri", "refractive-index", 2, "REAL IMAG", "Particle and physics",
         "Complex refractive index. REAL must be positive; IMAG is the absorption term."},
        {"w", "wavelength-um", 1, "LAMBDA", "Particle and physics",
         "Positive wavelength in micrometers."},
        {"n", "max-reflections", 1, "N", "Particle and physics",
         "Maximum internal reflection/refraction depth in [0, 30]; default 6."},
        {"abs", "absorption", 0, "", "Particle and physics",
         "Enable absorption accounting explicitly. A nonzero imaginary refractive index also enables it."},
        {"abs_points", "abs-points", 1, "N|all", "Particle and physics",
         "Absorption samples per polygon: positive integer or all."},
        {"sym", "symmetry", 2, "BETA_FACTOR GAMMA_FACTOR", "Particle and physics",
         "Override particle symmetry using positive integer beta and gamma factors in supported averaged PO modes."},

        // Geometry classification.
        {"geometry", "geometry", 1, "MODE", "Geometry",
         "Geometry classification: auto, convex, or nonconvex."},
        {"forced_nonconvex", "force-nonconvex", 0, "", "Geometry",
         "Compatibility switch equivalent to --geometry nonconvex."},
        {"forced_convex", "force-convex", 0, "", "Geometry",
         "Compatibility switch equivalent to --geometry convex."},
        {"r", "beam-area-ratio", 1, "RATIO", "Geometry",
         "Small-fragment area ratio, at least 1. The legacy default is 100; use --cutoff-profile off to disable this simplification."},
        {"save_geometry", "save-geometry", 1, "FILE", "Geometry",
         "Save the effective, scaled particle in the native text format and continue the calculation.",
         false, {"save-particle"}},

        // Orientation modes. Exactly one primary mode is required.
        {"fixed", "fixed-orientation", 2, "BETA_DEG GAMMA_DEG", "Orientation",
         "Trace one fixed beta/gamma orientation in degrees."},
        {"random", "euler-grid", 2, "NBETA NGAMMA", "Orientation",
         "Regular beta by gamma orientation grid. This replaces the misleading legacy name --random."},
        {"montecarlo", "monte-carlo", 1, "N", "Orientation",
         "Monte Carlo orientation average with N samples."},
        {"orientfile", "orientation-file", 1, "FILE", "Orientation",
         "Read beta and gamma orientations from FILE."},
        {"oldauto", "diffraction-grid", 1, "DIVISOR", "Orientation",
         "Physics-based regular orientation grid divided from the diffraction-limit estimate."},
        {"sobol", "sobol", 1, "N", "Orientation",
         "Sobol quasi-random orientation set with N samples."},
        {"so3_quat", "so3-quaternion", 1, "N", "Orientation",
         "Hammersley samples directly on full SO(3), represented as quaternions."},
        {"sobol_seed", "sobol-seed", 2, "N SEED", "Orientation",
         "Sobol orientations with an explicit nested Owen scramble seed."},
        {"sobol_ring", "sobol-ring", 2, "NBETA NGAMMA", "Orientation",
         "Sobol beta samples combined with shifted uniform gamma rings."},
        {"hammersley", "hammersley", 1, "N", "Orientation",
         "Hammersley low-discrepancy orientation set."},
        {"lattice", "lattice", 1, "N", "Orientation",
         "Rank-1 lattice orientation set with an automatically selected generator."},
        {"lattice_z", "lattice-generator", 2, "N Z", "Orientation",
         "Rank-1 lattice orientation set with explicit generator Z."},
        {"euler_quad", "euler-quadrature", 2, "NBETA NGAMMA", "Orientation",
         "Gauss quadrature in cos(beta) with periodic gamma samples."},
        {"euler_adapt", "euler-adaptive", 2, "NBETA NGAMMA_MAX", "Orientation",
         "Gauss beta quadrature with an adaptive gamma count per beta ring."},
        {"adaptive_euler", "adaptive-euler-grid", 1, "EPS", "Orientation",
         "Converge N_beta and N_gamma independently, repeat coordinate refinement, and verify a joint refinement."},
        {"adaptive", "adaptive-orientations", 1, "EPS", "Orientation",
         "Converge only Sobol N to EPS; theta, phi, and reflection depth stay fixed. At its standalone cap it warns and writes the best result instead of failing."},
        {"auto", "auto", 1, "EPS", "Orientation",
         "Jointly converge phi, theta, and Sobol orientation sampling to EPS while keeping reflection depth fixed."},
        {"autofull", "autofull", 1, "EPS", "Orientation",
         "Jointly converge reflection depth, phi, theta, and Sobol orientations to EPS.",
         false, {"fullauto"}},
        {"oldautofull", "diffraction-autofull", 1, "EPS", "Orientation",
         "Jointly converge reflection depth, phi, theta, N_beta, and N_gamma; finish on a regular Euler quadrature grid.",
         false, {"oldfullauto"}},
        {"adaptive_config", "adaptive-config", 1, "FILE", "Orientation",
         "Load per-parameter adaptive minima, maxima, and tolerances from FILE. Explicit limit flags override the file."},
        {"b", "beta-range-deg", 2, "MIN MAX", "Orientation",
         "Beta range in degrees for --euler-grid."},
        {"g", "gamma-range-deg", 2, "MIN MAX", "Orientation",
         "Gamma range in degrees for --euler-grid."},
        {"owen_avg", "owen-average", 1, "K", "Orientation",
         "Average K nested Owen final seeds in autofull modes."},
        {"owen_seeds", "owen-seeds", '+', "SEED...", "Orientation",
         "Explicit nested Owen seeds for autofull final averaging."},
        {"maxorient", "max-orientations", 1, "N", "Orientation",
         "Maximum adaptive orientation count; use a power of two."},
        {"max_beta_points", "max-beta-points", 1, "N", "Orientation",
         "Hard upper limit for adaptive N_beta; default 1024."},
        {"max_gamma_points", "max-gamma-points", 1, "N", "Orientation",
         "Hard upper limit for adaptive N_gamma; default 2400."},
        {"chunk", "orientation-chunk", 1, "N", "Orientation",
         "Maximum orientations or gamma values held in one memory chunk for supported averaged PO modes."},
        {"ring_points", "ring-points", 1, "N", "Orientation",
         "Samples per diffraction ring used by physics-based orientation estimates; default 3."},
        {"mirror_gamma", "mirror-gamma", 0, "", "Orientation",
         "Use mirror symmetry in supported averages: halve default gamma domains and treat an explicit gamma range as already reduced."},
        {"coh_orient", "coherent-orientations", 0, "", "Orientation",
         "Sum different orientations coherently; retained for legacy research runs.", true},
        {"pole", "pole", 0, "", "Orientation",
         "At exact beta poles, trace one gamma value and apply the full pole weight."},

        // Scattering angle grid.
        {"grid", "scattering-grid", '+', "T1 T2 NPHI NTH | R NPHI NTH", "Scattering grid",
         "Uniform theta range or backscatter cone. NTH is the interval count, so a theta range writes NTH+1 rows."},
        {"tgrid", "theta-grid-file", 1, "FILE", "Scattering grid",
         "Read a non-uniform theta grid in degrees from FILE."},
        {"auto_tgrid", "auto-theta-grid", 1, "EPS", "Scattering grid",
         "Converge theta interpolation using all Mueller elements and no particle-specific halo angles.",
         false, {"autotgrid"}},
        {"auto_phi", "auto-phi", 0, "", "Scattering grid",
         "Compatibility switch that converges N_phi to the default relative tolerance 0.02.",
         false, {"autophi"}},
        {"adaptive_phi", "adaptive-phi", 1, "EPS", "Scattering grid",
         "Converge the laboratory alpha average, represented by scattering-azimuth N_phi, to relative tolerance EPS while keeping theta, reflections, and orientation rule fixed.",
         false, {"adaptive-alpha"}},
        {"adaptive_reflections", "adaptive-reflections", 1, "EPS", "Scattering grid",
         "Converge only the internal reflection depth n to relative tolerance EPS."},
        {"nphi", "phi-points", 1, "N", "Scattering grid",
         "Explicit scattering-azimuth count used for the equivalent laboratory-alpha average. Overrides NPHI embedded in --scattering-grid.",
         false, {"alpha-points"}},
        {"max_theta_points", "max-theta-points", 1, "N", "Scattering grid",
         "Hard upper limit for adaptive theta points; default 4097."},
        {"max_phi_points", "max-phi-points", 1, "N", "Scattering grid",
         "Hard upper limit for adaptive N_phi/N_alpha; rounded down to a multiple of six; default 2400.",
         false, {"max-alpha-points"}},
        {"convergence_passes", "convergence-passes", 1, "N", "Scattering grid",
         "Required pass streak for adaptive orientation, phi, reflection, and Euler searches; default 2. Theta uses interval-error refinement instead."},
        {"filter", "backscatter-filter-deg", 1, "ANGLE", "Scattering grid",
         "Restrict output to a backscattering cone with ANGLE in degrees."},
        {"point", "backscatter-point", 0, "", "Scattering grid",
         "Unavailable legacy single-point mode; use a narrow backscatter cone instead.", true},

        // Execution backend and performance controls.
        {"backend", "backend", 1, "MODE", "Backend",
         "Execution backend: auto, cpu, or cuda."},
        {"gpu", "gpu", 0, "", "Backend",
         "Compatibility switch equivalent to --backend cuda."},
        {"cpu", "cpu", 0, "", "Backend",
         "Compatibility switch equivalent to --backend cpu."},
        {"threads", "threads", 1, "N", "Backend",
         "Positive OpenMP host worker count; default is the physical core count."},
        {"fft", "fft", 0, "", "Backend",
         "Use cuFFT angular interpolation with an automatically selected phi compression factor; requires PO and CUDA."},
        {"fft_factor", "fft-factor", 1, "N", "Backend",
         "Enable FFT interpolation and evaluate approximately N times fewer direct phi samples; N must be 1..64."},
        {"fft_tolerance", "fft-tolerance", 1, "EPS", "Backend",
         "Enable nested-grid FFT validation; if the estimated relative error exceeds EPS, use the doubled-resolution result."},
        {"gpu_trace", "gpu-trace-prefilter", 0, "", "Backend",
         "Use the experimental CUDA candidate prefilter for nonconvex tracing.", true},
        {"trace_prefilter", "trace-prefilter", 0, "", "Backend",
         "Enable the CPU projected-AABB candidate prefilter for nonconvex tracing.", true},
        {"no_trace_prefilter", "no-trace-prefilter", 0, "", "Backend",
         "Disable the CPU projected-AABB tracing prefilter.", true},
        {"trace_prefilter_margin", "trace-prefilter-margin", 1, "MARGIN", "Backend",
         "Nonnegative projected-AABB prefilter margin.", true},
        {"trace_prefilter_stats", "trace-prefilter-stats", 0, "", "Backend",
         "Print tracing prefilter candidate counters.", true},

        // Accuracy/performance cutoffs.
        {"cutoff_profile", "cutoff-profile", 1, "MODE", "Cutoffs",
         "Coherent cutoff preset: off, safe, balanced, or fast. Do not combine it with individual cutoff thresholds."},
        {"beam_cutoff", "beam-cutoff", 1, "EPS", "Cutoffs",
         "Legacy shortcut: reject an output beam when either relative Jones intensity OR relative area is below EPS."},
        {"beam_cutoff_j", "beam-cutoff-jones", 1, "EPS", "Cutoffs",
         "Skip output beams with relative squared Jones magnitude below EPS."},
        {"beam_cutoff_area", "beam-cutoff-area", 1, "EPS", "Cutoffs",
         "Skip output beams with relative projected area below EPS."},
        {"beam_cutoff_importance", "beam-cutoff-importance", 1, "EPS", "Cutoffs",
         "Skip output beams with relative Jones-intensity times area below EPS."},
        {"trace_cutoff", "trace-cutoff", 1, "EPS", "Cutoffs",
         "Legacy shortcut: prune a trace branch when either relative Jones intensity OR relative area is below EPS."},
        {"trace_cutoff_j", "trace-cutoff-jones", 1, "EPS", "Cutoffs",
         "Prune internal beam branches by relative squared Jones magnitude."},
        {"trace_cutoff_area", "trace-cutoff-area", 1, "EPS", "Cutoffs",
         "Prune internal beam branches by relative area."},
        {"trace_cutoff_importance", "trace-cutoff-importance", 1, "EPS", "Cutoffs",
         "Prune internal beam branches by relative Jones-intensity times area."},
        {"trace_max_beams", "trace-max-beams", 1, "N", "Cutoffs",
         "Abort an orientation after N traced beam nodes; zero disables the limit."},

        // Multi-size execution.
        {"multigrid", "dmax-grid", 3, "DMIN DMAX N", "Multi-size",
         "Log-spaced Dmax scan. Serial PO supports selected regular/Sobol rules; use --scan-jobs for other modes, with a file particle."},
        {"multikeq", "k-eq-grid", 3, "KMIN KMAX N", "Multi-size",
         "Log-spaced equivalent-size scan. Serial PO supports Euler-grid, diffraction-grid, and diffraction-autofull; use --scan-jobs otherwise."},
        {"multikeq_list", "k-eq-list", 1, "FILE", "Multi-size",
         "Read exact positive k_eq values from FILE."},
        {"multigrid_parallel", "scan-jobs", 1, "N", "Multi-size",
         "Run a size scan in N child processes; zero selects the count automatically."},
        {"multigrid_threads", "scan-threads", 1, "N", "Multi-size",
         "Positive OpenMP thread count per scan child; do not combine it with --threads."},
        {"gpu_devices", "gpu-devices", 1, "LIST", "Multi-size",
         "Comma-separated CUDA device indices assigned to scan children."},
        {"multikeq_shared_batches", "shared-k-eq-batches", 0, "", "Multi-size",
         "Batch nearby k_eq values per GPU and reuse tracing from the largest value."},
        {"multikeq_batch_ratio", "k-eq-batch-ratio", 1, "RATIO", "Multi-size",
         "Maximum kmax/kmin ratio in a shared batch; must be at least 1."},

        // Trajectory diagnostics.
        {"tr", "trajectories", 1, "FILE", "Trajectories",
         "Load a trajectory-selection file.", true},
        {"all", "all-trajectories", 0, "", "Trajectories",
         "Calculate all loaded trajectories.", true},
        {"gr", "trajectory-groups", 0, "", "Trajectories",
         "Write trajectory group diagnostics.", true},

        // Output and diagnostics.
        {"o", "output", 1, "PATH", "Output",
         "Output directory/result prefix. Optional placeholders use %[INDEX]OPTION_, for example %0p_."},
        {"close", "close", 0, "", "Output",
         "Exit after calculation instead of waiting for input; retained for compatibility."},
        {"log", "progress-interval", 1, "SECONDS", "Output",
         "Nonnegative progress-reporting interval in seconds."},
        {"save_betas", "save-betas", 0, "", "Output",
         "Write intermediate Mueller results for each beta ring."},
        {"checkpoint", "checkpoint", 0, "", "Output",
         "Enable checkpoint save/resume for supported long orientation runs."},
        {"jones", "jones-output", 0, "", "Output",
         "Write Jones matrices where the selected PO mode supports them."},
        {"shadow_off", "no-shadow-beam", 0, "", "Output",
         "Disable the external shadow beam in tracing."},
        {"noshadow_output", "no-shadow-output", 0, "", "Output",
         "Also calculate and write the Mueller result without the shadow contribution."},
        {"full_only", "full-only", 0, "", "Output",
         "Write only the full Mueller output; this is already the default.", true},
        {"shadow", "shadow", 0, "", "Output",
         "Deprecated compatibility flag with no effect.", true},
        {"ot_phase_avg", "ot-phase-average", 0, "", "Output",
         "Average optical-theorem extinction over a far-reference phase period.", true},
        {"ot_phase_shift", "ot-phase-shift", 1, "WAVELENGTHS", "Output",
         "Diagnostic optical-theorem far-reference phase shift in wavelengths.", true},
        {"ot_ping", "ot-ping-distance", 1, "DISTANCE", "Output",
         "Legacy far-screen optical-theorem phase correction distance.", true},

        {"help", "help", 0, "", "Help",
         "Show the supported production options and examples.", false, {"h"}},
        {"version", "version", 0, "", "Help",
         "Show the source revision and compiled CPU/CUDA feature set."},
        {"help-debug", "help-debug", 0, "", "Help",
         "Show production, legacy, diagnostic, and experimental options."}
    };
    return options;
}

const CliOptionSpec *FindCliOptionSpec(const std::string &keyOrAlias)
{
    const std::vector<CliOptionSpec> &options = GetCliOptionSpecs();
    for (const CliOptionSpec &spec : options)
    {
        if (spec.key == keyOrAlias || spec.canonical == keyOrAlias)
            return &spec;
        if (std::find(spec.aliases.begin(), spec.aliases.end(), keyOrAlias)
            != spec.aliases.end())
            return &spec;
    }
    return nullptr;
}

std::string CliCanonicalFlag(const std::string &keyOrAlias)
{
    const CliOptionSpec *spec = FindCliOptionSpec(keyOrAlias);
    return FlagToken(spec ? spec->canonical : keyOrAlias);
}

void RegisterCliOptions(ArgPP &parser)
{
    const std::vector<CliOptionSpec> &options = GetCliOptionSpecs();
    for (const CliOptionSpec &spec : options)
    {
        parser.AddRule(spec.key, spec.valueCount, true);
        if (spec.canonical != spec.key)
            parser.AddAlias(spec.canonical, spec.key);
        for (const std::string &alias : spec.aliases)
            if (alias != spec.key && alias != spec.canonical)
                parser.AddAlias(alias, spec.key);
    }
}

void PrintGeneratedHelp(bool includeDebugOptions)
{
    using std::cout;
    cout << "MBS-fast: physical and geometrical optics for faceted particles\n\n"
         << "Usage:\n"
         << "  mbs_po --method po --particle ... --refractive-index REAL IMAG \\\n"
         << "      --wavelength-um LAMBDA ORIENTATION [GRID] [options]\n\n"
         << "Required selections:\n"
         << "  * exactly one method: --method po or --method go\n"
         << "  * exactly one particle source: --particle or --particle-file\n"
         << "  * exactly one orientation mode\n"
         << "  * positive wavelength and an explicit refractive index\n\n"
         << "Canonical names are shown first. Compatibility aliases remain accepted.\n";

    const char *sectionOrder[] = {
        "Method", "Particle and physics", "Geometry", "Orientation",
        "Scattering grid", "Backend", "Cutoffs", "Multi-size",
        "Trajectories", "Output", "Help"
    };
    const std::vector<CliOptionSpec> &options = GetCliOptionSpecs();
    for (const char *section : sectionOrder)
    {
        bool printedHeader = false;
        for (const CliOptionSpec &spec : options)
        {
            if (spec.section != section || (spec.debugOnly && !includeDebugOptions))
                continue;
            if (!printedHeader)
            {
                cout << "\n" << section << ":\n";
                printedHeader = true;
            }
            std::string signature = "  " + FlagToken(spec.canonical);
            if (!spec.valueSyntax.empty())
                signature += " " + spec.valueSyntax;
            if (signature.size() < 34)
            {
                signature.append(34 - signature.size(), ' ');
                PrintWrapped(cout, signature, spec.description + AliasText(spec));
            }
            else
            {
                cout << signature << '\n';
                PrintWrapped(cout, std::string(34, ' '),
                             spec.description + AliasText(spec));
            }
        }
    }

    cout << "\nExamples:\n"
         << "  mbs_po --method po --particle 1 100 70 --refractive-index 1.3116 0 \\\n"
         << "      --wavelength-um 0.532 --max-reflections 12 --diffraction-grid 2 \\\n"
         << "      --scattering-grid 0 180 600 180 --backend auto --threads 16 \\\n"
         << "      --output out --close\n\n"
         << "  mbs_po --method go --particle-file examples/cube.particle --k-eq 20 \\\n"
         << "      --refractive-index 1.3116 0 --wavelength-um 0.532 \\\n"
         << "      --fixed-orientation 37 19 --scattering-grid 0 180 360 180 \\\n"
         << "      --backend cpu --output out --close\n\n"
         << "  mbs_po --method po --particle 1 100 70 --refractive-index 1.3116 0 \\\n"
         << "      --wavelength-um 0.532 --autofull 0.02 \\\n"
         << "      --adaptive-config configs/adaptive.example.conf \\\n"
         << "      --backend cpu --threads 64 --output converged --close\n";
}
