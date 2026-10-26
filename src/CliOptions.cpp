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
         "Read one beta/gamma pair in degrees per data line from FILE."},
        {"deduplicate_orientations", "deduplicate-orientations", 0, "", "Orientation",
         "Trace exact duplicate CPU PO orientation-file angles once and sum their original weights."},
        {"parallel_trace", "parallel-trace", 0, "", "Orientation",
         "Trace CPU PO orientation-file entries in parallel using --threads; default is serial tracing."},
        {"analytic_backscatter", "analytic-backscatter", 0, "", "Orientation",
         "Average physical return-strip phases analytically and sample full coherent M11 minus that control at exactly 180 degrees; real-index convex prisms, CPU/CUDA Sobol/Hammersley/lattice/Euler quadrature."},
        {"analytic_return_order", "analytic-return-order", 1, "K", "Orientation",
         "Include return-image self terms p+q <= K+1 in the analytic backscatter control; 1..8, default 1. Cross-family interference remains in the full numerical residual."},
        {"haar_alpha", "haar-alpha", 0, "", "Orientation",
         "Sample laboratory alpha with a third Owen-Sobol dimension in PO CPU/CUDA Sobol runs, producing full Haar orientation samples."},
        {"haar_mirror_audit", "haar-mirror-audit", 0, "", "Orientation",
         "Trace both partners of each half-gamma Haar point: (beta,gamma,alpha) and (beta,Gamma-gamma,-alpha). GO supports nonabsorbing built-in types 1/10 with analytic alpha. N is base points, total traces 2*N; compare with --mirror-gamma at N."},
        {"analytic_facet_average", "analytic-facet-average", 0, "", "Orientation",
         "Apply analytic reflected-facet and projected-shadow orientation kernels to M11 at all scattering angles; PO CPU/CUDA with --haar-alpha or full SO(3) quaternion samples."},
        {"analytic_shadow_control", "analytic-shadow-control", 1, "MODEL", "Orientation",
         "Shadow control for analytic facet averaging: facets (default), circular, or off; both shadow models retain full MBS corrections in the residual."},
        {"analytic_facet_samples", "analytic-facet-samples", 0, "", "Output",
         "Write per-orientation, per-angle raw and physical control values for --analytic-facet-average."},
        {"analytic_mean_cache", "analytic-mean-cache", 1, "FILE", "Orientation",
         "Reuse verified physical analytic means for identical geometry/index/wavelength/grid/models; process-locked and input-hashed cache."},
        {"analytic_mean_reference", "analytic-mean-reference-grid", 1, "FILE", "Orientation",
         "Compute/cache facet means on a fixed reference theta grid, selecting current rows from it. Every current theta must be present; requires --analytic-facet-average."},
        {"analytic_mean_threads", "analytic-mean-threads", 1, "N", "Orientation",
         "CPU threads for analytic mean preparation, independent of tracing threads; positive integer. Fullauto selects a separate CPU budget unless explicitly fixed."},
        {"analytic_control_weights", "analytic-control-weights", 1, "FILE", "Orientation",
         "Frozen independent-pilot reflection/shadow coefficients per theta; requires --analytic-facet-average."},
        {"orientation_pipeline", "orientation-pipeline", 0, "", "Orientation",
         "Overlap CPU preparation of the next Haar/SO3 chunk with current CUDA diffraction; uses two host chunks."},
        {"profile_phases", "profile-phases", 0, "", "Output",
         "Write CPU preparation, diffraction, analytic controls and overlap timings for weighted orientation runs."},
        {"analytic_azimuth_gaussian", "analytic-azimuth-gaussian", 0, "", "Orientation",
         "Experimental internal-beam Gaussian aperture control with exact Bessel azimuth means; requires --analytic-facet-average and preserves the full coherent residual."},
        {"diffraction_sampling", "diffraction-sampling", 1, "Q", "Orientation",
         "Use Q angular intervals per diffraction scale xi=0.69*lambda/lmax for both grids; lmax is the longest particle-facet edge. Individual flags may override either grid."},
        {"orientation_diffraction_sampling", "orientation-diffraction-sampling", 1, "Q", "Orientation",
         "Use Q orientation intervals per diffraction scale: delta_beta and delta_gamma are approximately xi/Q."},
        {"oldauto", "diffraction-grid", 1, "DIVISOR", "Orientation",
         "Legacy physics-based orientation grid. Prefer --diffraction-sampling or --orientation-diffraction-sampling.", true},
        {"sobol", "sobol", 1, "N", "Orientation",
         "Sobol quasi-random orientation set with N samples."},
        {"so3_quat", "so3-quaternion", 1, "N", "Orientation",
         "Symmetry-reduced Hammersley sampling of isotropic SO(3): N beta/gamma particle orientations, with the laboratory alpha integral evaluated by scattering azimuth."},
        {"so3_mirror_audit", "so3-mirror-audit", 0, "", "Orientation",
         "For an even full-gamma --so3-quaternion N validation run, trace N/2 half-domain Hammersley points and their explicit reflected partners; compare with N/2 --mirror-gamma."},
        {"so3_full_quat", "so3-full-quaternion", 1, "N", "Orientation",
         "Audit mode: Hammersley samples directly on the full SO(3) with Shoemake unit quaternions; particle symmetry is intentionally not used."},
        {"sobol_seed", "sobol-seed", 2, "N SEED", "Orientation",
         "Sobol orientations with an explicit nested Owen scramble seed."},
        {"sobol_offset", "sobol-offset", 1, "INDEX", "Orientation",
         "Start a Haar Owen-Sobol block after INDEX preceding points; requires --sobol-seed and --haar-alpha."},
        {"so3_rule", "so3-rule", 1, "FILE", "Orientation",
         "Full SO(3) weighted Haar cubature file: BETA_DEG GAMMA_DEG ALPHA_DEG WEIGHT. Nonnegative weights sum to one; no assumed geometry symmetry."},
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
         "Legacy joint convergence of reflection depth, phi, theta, and Sobol orientations to EPS."},
        {"fullauto", "fullauto", 1, "EPS", "Orientation",
         "Self-contained adaptive controller with explicit-setting priority. PO uses Haar/MLMC and interpolation; GO uses CPU, fixed solid-angle bins and analytic alpha. Eight scrambles and durable recovery."},
        {"go_fullauto_worker", "go-fullauto-worker", 0, "", "Output",
         "Internal GO controller product: full-precision finite bins and separate full Mueller delta matrices.", true},
        {"fullauto_initial_model", "fullauto-initial-model", 1, "MODE", "Orientation",
         "Automatic initial budget: empirical (default, completed-data size model) or xi (legacy fine beta/gamma proxy). Explicit fullauto-initial wins."},
        {"fullauto_depth_pilot", "fullauto-depth-pilot", 1, "N", "Orientation",
         "Initial per-seed paired-depth pilot, power of two, default min(512, production initial). Independent of production budget."},
        {"fullauto_depth_order", "fullauto-depth-order", 1, "MODE", "Orientation",
         "Depth tuning order: auto (sampling-first for CUDA, depth-first for CPU), sampling-first, depth-first. GO uses depth-first."},
        {"fullauto_depth_search", "fullauto-depth-search", 1, "MODE", "Orientation",
         "Select reflection depth with paired comparisons: on by default when max-reflections is omitted; off fixes the start depth. Explicit max-reflections remains fixed."},
        {"fullauto_depth_start", "fullauto-depth-start", 1, "N", "Orientation",
         "Automatic depth search start, default 6; requires no explicit max-reflections."},
        {"fullauto_depth_max", "fullauto-depth-max", 1, "N", "Orientation",
         "Automatic depth search ceiling, default 18, at most 28. Validation also evaluates ceiling+2."},
        {"fullauto_cpu_workers", "fullauto-cpu-workers", 1, "N", "Orientation",
         "Fullauto CPU worker processes, 1..8, default 1. Threads are per process; requires CPU backend."},
        {"fullauto_pilot", "fullauto-pilot", 1, "N", "Orientation",
         "Fullauto independent training/validation samples per seed; power of two, default 512."},
        {"fullauto_kernel", "fullauto-kernel-count", 1, "N", "Orientation",
         "Fullauto GPU tuning sample count; power of two, default 256."},
        {"fullauto_initial", "fullauto-initial", 1, "N", "Orientation",
         "Fixed fullauto initial samples per seed, power of two. Default completed-data size model with a 512 floor; bounded by the orientation cap."},
        {"fullauto_correction", "fullauto-min-correction", 1, "N", "Orientation",
         "Fullauto minimum samples per correction and seed; power of two, default 256."},
        {"fullauto_rounds", "fullauto-max-rounds", 1, "N", "Orientation",
         "Fullauto maximum joint refinement rounds; default 32. Exhaustion returns a nonzero code."},
        {"fullauto_theta_range", "fullauto-theta-range", 2, "TH1 TH2", "Orientation",
         "Angular range in degrees for automatic fullauto theta selection, default 0..180. No explicit point count; incompatible with a supplied scattering/theta grid."},
        {"fullauto_refine_window", "fullauto-refine-window", 3, "TH1 TH2 EPS", "Orientation",
         "Resume an existing PO checkpoint into an isolated angular subdirectory at a new fractional tolerance. Preserve the original command/output; only the selected window is refined and certified. Outside rows remain in the unchanged parent result."},
        {"fullauto_theta_audit", "fullauto-theta-audit-count", 1, "N", "Orientation",
         "Independent dense theta pilot per seed before production. Automatic theta mode defaults to min(pilot,256); explicit grid mode defaults to 0. Zero disables this additional pilot."},
        {"fullauto_theta_start", "fullauto-theta-start", 1, "N", "Orientation",
         "Fullauto initial interpolation support points; default 33, including endpoints."},
        {"fullauto_short_packets", "fullauto-short-packets", 1, "auto|off", "Orientation",
         "CUDA pipeline: cap automatic prepared chunks for short theta packets; explicit chunk and disabled resource search win."},
        {"fullauto_audit_angular", "fullauto-audit-angular", 1, "auto|on|off", "Orientation",
         "Independent PO audit: frozen production-trained angular residual control; auto requires predicted savings, on forces testing, off keeps direct audit. Accuracy gates unchanged."},
        {"fullauto_dispatch_threads", "fullauto-dispatch-threads", 1, "N", "Orientation",
         "CUDA PO only: fix CPU preparation and analytic-mean pools to N without changing checkpoint identities or ordered accumulation. Conflicts with explicit thread flags."},
        {"fullauto_dispatch_gpus", "fullauto-dispatch-gpus", 1, "LIST", "Orientation",
         "Compatible CUDA resume only: dispatch saved work to visible devices. Repeating an index twice (0,0) shares a GPU; requires --fullauto-dispatch-threads per worker and divides automatic RAM/VRAM budgets. Maximum 8 slots. Saved sampling is unchanged."},
        {"fullauto_gpus", "fullauto-gpus", 1, "LIST", "Orientation",
         "Fullauto visible GPU indices or auto (benchmark visible devices); respects CUDA_VISIBLE_DEVICES. Default 0."},
        {"fullauto_resources", "fullauto-resources", 1, "MODE", "Orientation",
         "Benchmark host threads/chunks and estimate saved block size: auto (default), off. Explicit settings win."},
        {"fullauto_rotation", "fullauto-rotation-symmetry", 1, "MODE", "Orientation",
         "Verify fixed-axis body rotations and paired all-Mueller equivalence: auto (default), off. Explicit --symmetry wins."},
        {"fullauto_cutoffs", "fullauto-cutoffs", 1, "MODE", "Orientation",
         "Coherent PO: balanced beams with matched off-minus-balanced correction, auto (default) or off; no preset search. Explicit cutoffs win. GO accepts explicit off/balanced with its own audit."},
        {"verified_rotation", "verified-rotation-symmetry", 0, "", "Orientation",
         "Worker: detect facet-preserving rotations about body z and half turn about x; never trust file symmetry metadata."},
        {"haar_rotation_audit", "haar-rotation-audit", 1, "AXIS", "Orientation",
         "Worker: explicitly trace equivalent body-rotation partners, gamma or beta, for paired validation."},
        {"fullauto_mirror", "fullauto-mirror", 1, "MODE", "Orientation",
         "Fullauto mirror reduction: auto, on, off. PO verifies geometry and a paired pilot (default auto). GO defaults off; explicit on/auto supports nonabsorbing built-in types 1/10 with parity reconstruction, auto falls back otherwise. --mirror-gamma implies on."},
        {"fullauto_angular_cache", "fullauto-angular-cache", 1, "MODE", "Orientation",
         "Reuse identical per-seed angular rows when fullauto refines theta: on (default) or off. Recompute only missing angles and rebase M11 to full-grid analytic means."},
        {"fullauto_controls", "fullauto-analytic-controls", 1, "MODE", "Orientation",
         "Facet controls: auto (default, fall back to complete sampled PO on analytic mean limits), on, off. An explicit --analytic-facet-average implies on."},
        {"oldautofull", "diffraction-autofull", 1, "EPS", "Orientation",
         "Jointly converge reflection depth, phi, theta, N_beta, and N_gamma; finish on a regular Euler quadrature grid.",
         false, {"oldfullauto"}},
        {"fullauto_resume_compatible", "fullauto-resume-compatible", 1, "OLD_BINARY", "Orientation",
         "Resume an existing 8ae9e85, 89f83a3 or 336ef78 result after verified native-kernel compatible controller fixes. Requires the original binary, identical build profile, inputs and environment."},
        {"fullauto_workers", "fullauto-workers", 1, "MODE", "Orientation",
         "Worker lifetime: auto (persistent CUDA workers), persistent, process. Each completed block remains checkpointed; explicit process restores per-job exec."},
        {"fullauto_trace", "fullauto-trace-backend", 1, "MODE", "Orientation",
         "PO tracing: auto (validate and time CPU/resident/hybrid), cpu, resident, hybrid. Explicit choice is preserved; unsupported modes retain CPU preparation."},
        {"po_prepare_threads", "po-prepare-threads", 1, "N", "Orientation",
         "PO orientation preparation threads only; ordered matrix accumulation is unchanged. Explicit setting wins over fullauto selection."},
        {"trace_cache_mb", "trace-cache-mib", 1, "MIB", "Orientation",
         "Cap each process-local CPU prepared-beam and GPU raw-trace cache independently. Persistent defaults to at most 256 MiB each, standalone to 0. Zero disables both."},
        {"po_trace_backend", "po-trace-backend", 1, "MODE", "Orientation",
         "PO worker tracing: cpu, resident, hybrid. Resident/hybrid require CUDA; eligible FP64 convex coherent jobs use GPU tracing, others retain CPU preparation."},
        {"fullauto_work_queue", "fullauto-work-queue", 1, "MODE", "Orientation",
         "Dispatch independent orientation blocks, angle groups and levels in one worker queue: auto (default; CUDA, including one GPU, or multiple CPU workers), on, off. Sampling and tolerances are unchanged."},
        {"fullauto_incremental", "fullauto-incremental", 1, "MODE", "Orientation",
         "Reuse exact Owen-Sobol prefixes and checkpoint completed orientation blocks: on (default), off."},
        {"fullauto_orientation_block", "fullauto-orientation-block", 1, "N", "Orientation",
         "Maximum saved orientation block per production seed; positive power of two, default 65536."},
        {"fullauto_local", "fullauto-local-orientations", 1, "MODE", "Orientation",
         "Per-theta multilevel Sobol budgets: auto (default, cost-gated), on (force local), off (global reference)."},
        {"fullauto_ratios", "fullauto-ratio-accuracy", 1, "MODE", "Orientation",
         "Also bound normalized Mij/M11 using paired Fieller intervals: on (default), off."},
        {"fullauto_phi_recheck", "fullauto-phi-recheck", 1, "MODE", "Orientation",
         "Recheck phi after the dense pilot and, once per depth, on persistently noisy working angles with independent cost-gated confirmation: on (default), off. Explicit phi/pilot settings are retained."},
        {"fullauto_polar_controls", "fullauto-polarization-controls", 1, "MODE", "Orientation",
         "Frozen scalar analytic control variates for all Mueller entries, retained only after independent variance validation: on (default), off."},
        {"fullauto_final_audit", "fullauto-final-audit", 1, "MODE", "Orientation",
         "Independent final seed validation of the frozen estimator and interpolation: on (default), off."},
        {"fullauto_final_checks", "fullauto-final-checks", 1, "MODE", "Orientation",
         "Final verification policy: on (default) retains configured audits; off skips independent seeds and final depth+2/cutoff physics audits, including automatic-depth final verification. Production convergence and startup depth selection remain active. May change on an unfinished checkpoint; skipped checks are reported as unverified."},
        {"fullauto_physics_audit", "fullauto-physics-audit", 1, "MODE", "Orientation",
         "Compare target depth with depth+2 and requested cutoffs with off on representative angles: on, off (default). Does not replace requested physics."},
        {"fullauto_depth_share", "fullauto-depth-sharing", 1, "MODE", "Orientation",
         "Share convex beam tracing for paired depths after an all-Mueller equivalence probe: on (default), off. Unsupported geometry uses independent traces."},
        {"orientation_prefix_output", "orientation-prefix-output", 1, "M", "Orientation",
         "Also save the first M equally weighted seeded Haar orientations during this native PO pass. Internal balanced-correction worker option; no additional tracing."},
        {"paired_depth", "paired-depth", 1, "LOWER", "Orientation",
         "Also output LOWER-depth coherent matrices from a shared convex trace; requires seeded Haar Sobol, no branch pruning or diagnostic samples."},
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
         "Legacy density used by --diffraction-grid and adaptive orientation estimates; default 3. It does not affect the new diffraction-sampling flags.", true},
        {"mirror_gamma", "mirror-gamma", 0, "", "Orientation",
         "Assert verified numerical mirror symmetry in supported averages, including --so3-quaternion: sample half the default gamma domain and restore the omitted Mueller contribution with its Stokes parity."},
        {"coh_orient", "coherent-orientations", 0, "", "Orientation",
         "Disabled legacy option: coherent averaging across unrelated orientations has no defined relative phase.", true},
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
        {"diffraction_limit_grid", "diffraction-limit-grid", 1, "FACTOR", "Scattering grid",
         "Legacy scattering-grid factor: step=FACTOR*0.69*lambda/lmax, where lmax is the longest particle-facet edge. Prefer --diffraction-sampling or --scattering-diffraction-sampling.", true},
        {"scattering_diffraction_sampling", "scattering-diffraction-sampling", 1, "Q", "Scattering grid",
         "Use Q scattering-direction intervals per diffraction scale: delta_theta and equatorial delta_phi are approximately xi/Q."},
        {"latitude_phi_grid", "latitude-phi-grid", 0, "", "Scattering grid",
         "Use N_phi(theta) proportional to sin(theta), preserving the equatorial diffraction-limit spacing and avoiding pole oversampling."},
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
         "Force CUDA facet ordering and candidate prefiltering for nonconvex tracing; CUDA PO builds enable it automatically.", true},
        {"no_gpu_trace", "no-gpu-trace-prefilter", 0, "", "Backend",
         "Disable automatic CUDA facet ordering and candidate prefiltering.", true},
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
        {"trace_limit_retries", "trace-limit-retries", 1, "N", "Cutoffs",
         "Retry an orientation that reaches --trace-max-beams N times, preserving all accuracy cutoffs; later orientations reuse the lowest successful allowed retry level."},
        {"trace_retry_factor", "trace-retry-factor", 1, "FACTOR", "Cutoffs",
         "Multiply --trace-max-beams by FACTOR after each retry; must exceed 1, default 10."},

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
        {"extinction_only", "extinction-only", 0, "", "Output",
         "Compute CPU Hammersley forward-OT extinction only: reuse convex ray geometry for all sizes, omit the angular Mueller grid. Ray-energy albedo is explicitly approximate."},
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
        {"allow_experimental_environment", "allow-experimental-environment", 0, "", "Output",
         "Acknowledge and permit active MBS_* environment overrides that can change numerical or physical results; every active variable is recorded in the result log.", true},

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

const char *FullAutoOptionPolicy(const std::string &key)
{
    static const std::set<std::string> controller={"fullauto","fullauto_refine_window","fullauto_resume_compatible","fullauto_initial_model","fullauto_depth_pilot","fullauto_depth_order","fullauto_depth_search","fullauto_depth_start","fullauto_depth_max","fullauto_cpu_workers","fullauto_pilot","fullauto_kernel","fullauto_initial","fullauto_correction","fullauto_rounds","fullauto_theta_range","fullauto_theta_audit","fullauto_theta_start","fullauto_gpus","fullauto_dispatch_gpus","fullauto_dispatch_threads","fullauto_audit_angular","fullauto_short_packets","fullauto_resources","fullauto_rotation","fullauto_cutoffs","fullauto_mirror","fullauto_angular_cache","fullauto_controls","fullauto_workers","fullauto_trace","fullauto_work_queue","fullauto_incremental","fullauto_orientation_block","fullauto_local","fullauto_ratios","fullauto_phi_recheck","fullauto_polar_controls","fullauto_final_audit","fullauto_final_checks","fullauto_physics_audit","fullauto_depth_share"};
    if(controller.count(key))return "controller";
    static const std::set<std::string> fixed={"trace_cache_mb","po_prepare_threads","po_trace_backend","method","po","go","incoh","karczewski","legacy_sign","p","pf","rs","k_eq","ri","w","n","abs","abs_points","sym","geometry","forced_nonconvex","forced_convex","r","haar_alpha","analytic_facet_average","analytic_shadow_control","analytic_mean_cache","analytic_mean_threads","mirror_gamma","grid","tgrid","nphi","backend","gpu","cpu","threads","chunk","orientation_pipeline","gpu_trace","no_gpu_trace","trace_prefilter","no_trace_prefilter","trace_prefilter_margin","trace_prefilter_stats","cutoff_profile","beam_cutoff","beam_cutoff_j","beam_cutoff_area","beam_cutoff_importance","trace_cutoff","trace_cutoff_j","trace_cutoff_area","trace_cutoff_importance","trace_max_beams","trace_limit_retries","trace_retry_factor","o","close","log","shadow_off","full_only","save_geometry","ot_phase_avg","ot_phase_shift","ot_ping","allow_experimental_environment","owen_seeds"};
    static const std::set<std::string> limits={"maxorient","max_theta_points","max_phi_points","convergence_passes","adaptive_config"};
    if(fixed.count(key))return "fixed";
    if(limits.count(key))return "limit";
    if(key=="checkpoint"||key=="profile_phases")return "native equivalent";
    if(key=="help"||key=="help-debug"||key=="version")return "help";
    return "unsupported";
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
         << "      --wavelength-um 0.532 --max-reflections 12 --diffraction-sampling 2 \\\n"
         << "      --latitude-phi-grid --backend auto --threads 16 \\\n"
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
