#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ArgPP.h"
#include "CliOptions.h"
#include "IntegralCharacteristics.h"
#include "RunConfig.h"

namespace
{

int failures = 0;

RunConfig ParseConfig(const std::vector<std::string> &arguments,
                      bool gpuDefault = false,
                      bool cudaCompiled = false)
{
    std::vector<std::string> storage;
    storage.push_back("mbs_po_test");
    storage.insert(storage.end(), arguments.begin(), arguments.end());
    std::vector<const char *> argv;
    for (const std::string &arg : storage)
        argv.push_back(arg.c_str());

    ArgPP parser;
    RegisterCliOptions(parser);
    parser.Parse((int)argv.size(), argv.data());
    return RunConfig::FromCommandLine(parser, gpuDefault, cudaCompiled);
}

void Pass(const std::string &name)
{
    std::cout << "PASS: " << name << '\n';
}

void FailTest(const std::string &name, const std::string &message)
{
    ++failures;
    std::cerr << "FAIL: " << name << ": " << message << '\n';
}

void ExpectSuccess(const std::string &name,
                   const std::vector<std::string> &arguments)
{
    try
    {
        (void)ParseConfig(arguments);
        Pass(name);
    }
    catch (const std::exception &ex)
    {
        FailTest(name, ex.what());
    }
}

void ExpectFailure(const std::string &name,
                   const std::vector<std::string> &arguments,
                   const std::string &expectedText)
{
    try
    {
        (void)ParseConfig(arguments);
        FailTest(name, "command unexpectedly passed validation");
    }
    catch (const std::exception &ex)
    {
        const std::string message(ex.what());
        if (message.find(expectedText) == std::string::npos)
            FailTest(name, "expected '" + expectedText + "' in: " + message);
        else if (message.find("Fix:") == std::string::npos)
            FailTest(name, "error is not actionable: " + message);
        else
            Pass(name);
    }
}

void ExpectFailureForBuild(const std::string &name,
                           const std::vector<std::string> &arguments,
                           bool gpuDefault, bool cudaCompiled,
                           const std::string &expectedText)
{
    try
    {
        (void)ParseConfig(arguments, gpuDefault, cudaCompiled);
        FailTest(name, "command unexpectedly passed validation");
    }
    catch (const std::exception &ex)
    {
        const std::string message(ex.what());
        if (message.find(expectedText) == std::string::npos)
            FailTest(name, "expected '" + expectedText + "' in: " + message);
        else if (message.find("Fix:") == std::string::npos)
            FailTest(name, "error is not actionable: " + message);
        else
            Pass(name);
    }
}

void ExpectSuccessForBuild(const std::string &name,
                           const std::vector<std::string> &arguments,
                           bool gpuDefault, bool cudaCompiled)
{
    try
    {
        (void)ParseConfig(arguments, gpuDefault, cudaCompiled);
        Pass(name);
    }
    catch (const std::exception &ex)
    {
        FailTest(name, ex.what());
    }
}

std::string TestRoot()
{
    const char *root = std::getenv("MBS_TEST_ROOT");
    if (!root || !*root)
        throw std::runtime_error("MBS_TEST_ROOT is not set");
    return root;
}

std::string WriteAdaptiveTestFile(const std::string &name,
                                  const std::string &contents)
{
    const std::string path = TestRoot() + "/tests/.build/" + name;
    std::ofstream output(path.c_str());
    output << contents;
    output.close();
    return path;
}

std::vector<std::string> AdaptiveBase(const std::string &configPath)
{
    return {
        "--method", "po",
        "--particle", "1", "10", "10",
        "--refractive-index", "1.31", "0",
        "--wavelength-um", "0.532",
        "--autofull", "0.2",
        "--adaptive-config", configPath,
        "--backend", "cpu"
    };
}

std::vector<std::string> CanonicalBase();

void CheckAdaptiveConfig()
{
    const std::string example = TestRoot()
        + "/configs/adaptive.example.conf";
    try
    {
        RunConfig parsed = ParseConfig(AdaptiveBase(example));
        if (!parsed.adaptive.loadedFromFile
            || parsed.adaptive.minPhiPoints != 12
            || parsed.adaptive.maxPhiPoints != 2400
            || parsed.adaptive.minOrientations != 128
            || parsed.adaptive.maxOrientations != 262144
            || std::fabs(parsed.adaptive.thetaTolerance - 0.01) > 1e-12)
        {
            FailTest("adaptive config values", "loaded values do not match the example");
        }
        else
            Pass("adaptive config values");
    }
    catch (const std::exception &ex)
    {
        FailTest("adaptive config values", ex.what());
    }

    try
    {
        std::vector<std::string> command = AdaptiveBase(example);
        command.insert(command.end(), {
            "--max-reflections", "12",
            "--max-alpha-points", "600",
            "--max-orientations", "1024",
            "--convergence-passes", "1"
        });
        RunConfig parsed = ParseConfig(command);
        if (parsed.adaptive.maxReflections != 12
            || parsed.adaptive.maxPhiPoints != 600
            || parsed.adaptive.maxOrientations != 1024
            || parsed.adaptive.stablePasses != 1)
        {
            FailTest("adaptive CLI precedence", "explicit limits did not override the file");
        }
        else
            Pass("adaptive CLI precedence");
    }
    catch (const std::exception &ex)
    {
        FailTest("adaptive CLI precedence", ex.what());
    }

    const std::string typo = WriteAdaptiveTestFile(
        "adaptive_typo.conf",
        "alpha.min_points = 12\nalpha.max_pionts = 120\n");
    ExpectFailure("adaptive config typo", AdaptiveBase(typo),
                  "alpha.max_points");

    const std::string badRange = WriteAdaptiveTestFile(
        "adaptive_bad_range.conf",
        "theta.min_points = 65\ntheta.max_points = 33\n");
    ExpectFailure("adaptive config range", AdaptiveBase(badRange),
                  "cannot be refined");

    const std::string duplicate = WriteAdaptiveTestFile(
        "adaptive_duplicate.conf",
        "alpha.min_points = 12\nalpha.min_points = 18\n");
    ExpectFailure("adaptive config duplicate", AdaptiveBase(duplicate),
                  "repeats key");

    const std::string excessiveDepth = WriteAdaptiveTestFile(
        "adaptive_excessive_depth.conf",
        "reflections.min = 2\nreflections.max = 31\n");
    ExpectFailure("adaptive reflection history limit",
                  AdaptiveBase(excessiveDepth),
                  "supported limit of 30");

    ExpectFailure("missing adaptive config",
                  AdaptiveBase("/definitely/not/adaptive.conf"),
                  "cannot read adaptive config");

    std::vector<std::string> wrongMode = CanonicalBase();
    wrongMode.insert(wrongMode.end(), {"--adaptive-config", example});
    ExpectFailure("adaptive config wrong mode", wrongMode,
                  "does not affect the selected orientation mode");
}

std::vector<std::string> CanonicalBase()
{
    return {
        "--method", "po",
        "--particle", "1", "10", "10",
        "--refractive-index", "1.31", "0",
        "--wavelength-um", "0.532",
        "--max-reflections", "6",
        "--fixed-orientation", "0", "0",
        "--scattering-grid", "0", "180", "12", "24",
        "--backend", "cpu",
        "--close"
    };
}

std::vector<std::string> LegacyBase()
{
    return {
        "--po",
        "-p", "1", "10", "10",
        "--ri", "1.31", "0",
        "-w", "0.532",
        "-n", "6",
        "--fixed", "0", "0",
        "--grid", "0", "180", "12", "24",
        "--cpu",
        "--close"
    };
}

void CheckOptionSchema()
{
    const std::vector<CliOptionSpec> &specs = GetCliOptionSpecs();
    std::set<std::string> names;
    bool valid = !specs.empty();
    std::string problem;
    for (const CliOptionSpec &spec : specs)
    {
        if (spec.key.empty() || spec.canonical.empty() || spec.section.empty()
            || spec.description.empty())
        {
            valid = false;
            problem = "empty schema field for key '" + spec.key + "'";
            break;
        }

        std::vector<std::string> optionNames;
        optionNames.push_back(spec.key);
        if (spec.canonical != spec.key)
            optionNames.push_back(spec.canonical);
        optionNames.insert(optionNames.end(), spec.aliases.begin(), spec.aliases.end());
        for (const std::string &name : optionNames)
        {
            if (!names.insert(name).second)
            {
                valid = false;
                problem = "duplicate flag name '" + name + "'";
                break;
            }
        }
        if (!valid)
            break;
    }

    std::ostringstream generated;
    std::streambuf *oldBuffer = std::cout.rdbuf(generated.rdbuf());
    PrintGeneratedHelp(true);
    std::cout.rdbuf(oldBuffer);
    for (const CliOptionSpec &spec : specs)
    {
        const std::string flag = (spec.canonical.size() == 1 ? "-" : "--")
            + spec.canonical;
        if (generated.str().find(flag) == std::string::npos)
        {
            valid = false;
            problem = "generated help is missing '" + flag + "'";
            break;
        }
    }

    if (valid)
        Pass("option schema is unique and fully represented in help");
    else
        FailTest("option schema", problem);
}

} // namespace

int main()
{
    CheckOptionSchema();
    CheckAdaptiveConfig();
    ExpectSuccess("canonical command", CanonicalBase());
    ExpectSuccess("legacy aliases", LegacyBase());
    ExpectSuccess("adaptive Euler limits", {
        "--method", "po",
        "--particle", "1", "10", "10",
        "--refractive-index", "1.31", "0",
        "--wavelength-um", "0.532",
        "--max-reflections", "12",
        "--adaptive-euler-grid", "0.02",
        "--scattering-grid", "0", "180", "48", "180",
        "--max-beta-points", "512",
        "--max-gamma-points", "1200",
        "--max-orientations", "65536",
        "--backend", "cpu"
    });
    ExpectSuccess("fullauto compatibility alias", {
        "--method", "po",
        "--particle", "1", "10", "10",
        "--refractive-index", "1.31", "0",
        "--wavelength-um", "0.532",
        "--fullauto", "0.02",
        "--backend", "cpu"
    });
    ExpectSuccess("alpha quadrature aliases", {
        "--method", "po",
        "--particle", "1", "10", "10",
        "--refractive-index", "1.31", "0",
        "--wavelength-um", "0.532",
        "--max-reflections", "6",
        "--sobol", "64",
        "--scattering-grid", "0", "180", "12", "24",
        "--adaptive-alpha", "0.1",
        "--max-alpha-points", "600",
        "--backend", "cpu"
    });

    ExpectSuccess("dependency order is irrelevant", {
        "--method", "go",
        "--resize-dmax-um", "10",
        "--particle-file", "/dev/null",
        "--refractive-index", "1.31", "0",
        "--wavelength-um", "0.532",
        "--fixed-orientation", "0", "0",
        "--scattering-grid", "0", "180", "12", "24",
        "--backend", "cpu"
    });

    std::vector<std::string> missingMethod = CanonicalBase();
    missingMethod.erase(missingMethod.begin(), missingMethod.begin() + 2);
    ExpectFailure("missing method", missingMethod,
                  "calculation method is missing");

    std::vector<std::string> twoMethods = LegacyBase();
    twoMethods.insert(twoMethods.begin() + 1, "--go");
    ExpectFailure("conflicting methods", twoMethods,
                  "calculation method is ambiguous");

    std::vector<std::string> twoOrientations = CanonicalBase();
    twoOrientations.insert(twoOrientations.end(), {"--sobol", "64"});
    ExpectFailure("conflicting orientations", twoOrientations,
                  "orientation mode is ambiguous");

    std::vector<std::string> twoGeometries = CanonicalBase();
    twoGeometries.insert(twoGeometries.end(), {
        "--force-convex", "--forced_nonconvex"
    });
    ExpectFailure("conflicting geometry", twoGeometries,
                  "geometry classification is ambiguous");

    std::vector<std::string> twoThetaGrids = CanonicalBase();
    twoThetaGrids.insert(twoThetaGrids.end(), {
        "--auto-theta-grid", "0.05"
    });
    ExpectFailure("conflicting theta grids", twoThetaGrids,
                  "theta grid is ambiguous");

    std::vector<std::string> nanWavelength = CanonicalBase();
    for (size_t i = 0; i + 1 < nanWavelength.size(); ++i)
        if (nanWavelength[i] == "--wavelength-um")
            nanWavelength[i + 1] = "nan";
    ExpectFailure("NaN rejected", nanWavelength,
                  "expects a finite number");

    std::vector<std::string> badCutoff = CanonicalBase();
    badCutoff.insert(badCutoff.end(), {"--beam-cutoff", "1.1"});
    ExpectFailure("relative cutoff range", badCutoff,
                  "must be in [0, 1]");

    std::vector<std::string> badAreaRatio = CanonicalBase();
    badAreaRatio.insert(badAreaRatio.end(), {"--beam-area-ratio", "0.5"});
    ExpectFailure("beam area ratio range", badAreaRatio, "at least 1");

    std::vector<std::string> badProfile = CanonicalBase();
    badProfile.insert(badProfile.end(), {"--cutoff-profile", "turbo"});
    ExpectFailure("unknown cutoff profile", badProfile,
                  "unknown --cutoff-profile");

    std::vector<std::string> profileConflict = CanonicalBase();
    profileConflict.insert(profileConflict.end(), {
        "--cutoff-profile", "safe", "--trace-cutoff-importance", "1e-8"
    });
    ExpectFailure("cutoff profile conflict", profileConflict,
                  "cannot be combined");

    std::vector<std::string> profileOff = CanonicalBase();
    profileOff.insert(profileOff.end(), {"--cutoff-profile", "off"});
    ExpectSuccess("cutoff profile off", profileOff);

    std::vector<std::string> saveGeometry = CanonicalBase();
    saveGeometry.insert(saveGeometry.end(), {
        "--save-geometry", "saved.particle"
    });
    ExpectSuccess("save geometry option", saveGeometry);

    std::vector<std::string> badFftFactor = CanonicalBase();
    badFftFactor.insert(badFftFactor.end(), {"--fft-factor", "65"});
    ExpectFailure("FFT factor range", badFftFactor, "from 1 to 64");

    std::vector<std::string> badFftTolerance = CanonicalBase();
    badFftTolerance.insert(badFftTolerance.end(), {"--fft-tolerance", "1"});
    ExpectFailure("FFT tolerance range", badFftTolerance, "in (0, 1)");

    std::vector<std::string> fftCuda = CanonicalBase();
    for (size_t i = 0; i + 1 < fftCuda.size(); ++i)
    {
        if (fftCuda[i] == "--fixed-orientation")
        {
            fftCuda.erase(fftCuda.begin() + i, fftCuda.begin() + i + 3);
            fftCuda.insert(fftCuda.end(), {"--sobol", "16"});
            break;
        }
    }
    for (size_t i = 0; i + 1 < fftCuda.size(); ++i)
        if (fftCuda[i] == "--backend")
            fftCuda[i + 1] = "cuda";
    fftCuda.insert(fftCuda.end(), {
        "--fft-factor", "3", "--fft-tolerance", "0.02"
    });
    ExpectSuccessForBuild("explicit FFT controls", fftCuda, true, true);

    const IntegralCharacteristics absorbing = ComputeIntegralCharacteristics(
        IntegralMethod::PhysicalOptics, 100.0, 150.0, 170.0,
        true, true, 0.002);
    if (std::fabs(absorbing.cAbs - 20.0) > 1e-12
        || std::fabs(absorbing.qExt - 1.7) > 1e-12
        || std::fabs(absorbing.cScaErrorEstimateRel - 0.002) > 1e-12)
        FailTest("absorbing integral balance",
                 "shared characteristics formula is inconsistent");
    else
        Pass("absorbing integral balance");

    const IntegralCharacteristics goAbsorbing = ComputeIntegralCharacteristics(
        IntegralMethod::GeometricalOptics, 100.0, 80.0, 0.0,
        false, true, std::numeric_limits<double>::quiet_NaN());
    if (std::fabs(goAbsorbing.cExt - 100.0) > 1e-12
        || std::fabs(goAbsorbing.cSca - 80.0) > 1e-12
        || std::fabs(goAbsorbing.cAbs - 20.0) > 1e-12)
        FailTest("GO integral balance",
                 "GO characteristics do not use projected-area balance");
    else
        Pass("GO integral balance");

    const std::string compactLog = FormatIntegralCharacteristicsLog(
        absorbing, "full");
    if (compactLog.find("legacy") != std::string::npos
        || compactLog.find("C_abs_GO") != std::string::npos
        || compactLog.find("C_ext_OT") != std::string::npos)
        FailTest("compact integral log", "alternate definitions leaked into output");
    else
        Pass("compact integral log");

    const std::vector<double> theta = {
        0.0, M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0, M_PI
    };
    const std::vector<double> constantM11(theta.size(), 1.0);
    const double quadratureError = EstimateAngularIntegralRelativeError(
        theta, constantM11, 4.0*M_PI);
    if (!std::isfinite(quadratureError) || quadratureError > 1e-12)
        FailTest("Csca quadrature estimate",
                 "constant angular integral should be exact");
    else
        Pass("Csca quadrature estimate");

    std::vector<std::string> badAdaptiveTolerance = CanonicalBase();
    for (size_t i = 0; i + 2 < badAdaptiveTolerance.size(); ++i)
    {
        if (badAdaptiveTolerance[i] == "--fixed-orientation")
        {
            badAdaptiveTolerance.erase(badAdaptiveTolerance.begin() + i,
                                       badAdaptiveTolerance.begin() + i + 3);
            badAdaptiveTolerance.insert(badAdaptiveTolerance.end(),
                                        {"--adaptive-orientations", "1.2"});
            break;
        }
    }
    ExpectFailure("adaptive tolerance range", badAdaptiveTolerance,
                  "must be in (0, 1)");

    std::vector<std::string> smallThetaLimit = CanonicalBase();
    for (size_t i = 0; i + 2 < smallThetaLimit.size(); ++i)
    {
        if (smallThetaLimit[i] == "--fixed-orientation")
        {
            smallThetaLimit.erase(smallThetaLimit.begin() + i,
                                  smallThetaLimit.begin() + i + 3);
            smallThetaLimit.insert(smallThetaLimit.end(), {"--sobol", "16"});
            break;
        }
    }
    for (size_t i = 0; i + 4 < smallThetaLimit.size(); ++i)
    {
        if (smallThetaLimit[i] == "--scattering-grid")
        {
            smallThetaLimit.erase(smallThetaLimit.begin() + i,
                                  smallThetaLimit.begin() + i + 5);
            break;
        }
    }
    smallThetaLimit.insert(smallThetaLimit.end(), {
        "--auto-theta-grid", "0.1", "--max-theta-points", "16"
    });
    ExpectFailure("adaptive theta safety limit", smallThetaLimit,
                  "must be at least 17");

    std::vector<std::string> missingFile = CanonicalBase();
    missingFile.erase(missingFile.begin() + 2, missingFile.begin() + 6);
    missingFile.insert(missingFile.begin() + 2, {
        "--particle-file", "/definitely/not/a/particle.dat"
    });
    ExpectFailure("missing input file", missingFile,
                  "cannot read --particle-file");

    std::vector<std::string> duplicateAlias = CanonicalBase();
    duplicateAlias.insert(duplicateAlias.end(), {"--ri", "1.31", "0"});
    ExpectFailure("canonical plus legacy duplicate", duplicateAlias,
                  "duplicate option");

    std::vector<std::string> typo = CanonicalBase();
    typo.insert(typo.end(), {"--therads", "4"});
    ExpectFailure("unknown flag suggestion", typo,
                  "Did you mean");

    std::vector<std::string> emptyThetaFile = CanonicalBase();
    for (size_t i = 0; i + 4 < emptyThetaFile.size(); ++i)
    {
        if (emptyThetaFile[i] == "--scattering-grid")
        {
            emptyThetaFile.erase(emptyThetaFile.begin() + i,
                                 emptyThetaFile.begin() + i + 5);
            emptyThetaFile.insert(emptyThetaFile.end(), {
                "--theta-grid-file", "/dev/null"
            });
            break;
        }
    }
    ExpectFailure("empty theta file", emptyThetaFile,
                  "at least two distinct theta values");

    std::vector<std::string> goUnsupported = CanonicalBase();
    goUnsupported[1] = "go";
    for (size_t i = 0; i + 2 < goUnsupported.size(); ++i)
    {
        if (goUnsupported[i] == "--fixed-orientation")
        {
            goUnsupported.erase(goUnsupported.begin() + i,
                                goUnsupported.begin() + i + 3);
            goUnsupported.insert(goUnsupported.end(), {"--hammersley", "64"});
            break;
        }
    }
    ExpectFailure("unsupported GO orientation", goUnsupported,
                  "not implemented for GO");

    std::vector<std::string> goIncoherent = CanonicalBase();
    goIncoherent[1] = "go";
    goIncoherent.push_back("--incoherent");
    ExpectFailure("GO rejects ignored incoherent mode", goIncoherent,
                  "not consumed by the GO implementation");

    std::vector<std::string> goBeamCutoff = CanonicalBase();
    goBeamCutoff[1] = "go";
    goBeamCutoff.insert(goBeamCutoff.end(), {"--beam-cutoff", "0.01"});
    ExpectFailure("GO rejects ignored output-beam cutoff", goBeamCutoff,
                  "implemented only by PO diffraction");

    ExpectSuccess("fixed reflection depth does not constrain theta search", {
        "--method", "po", "--particle", "1", "1", "1",
        "--refractive-index", "1.31", "0", "--wavelength-um", "10",
        "--max-reflections", "1", "--sobol", "16",
        "--auto-theta-grid", "0.1", "--backend", "cpu"
    });

    std::vector<std::string> cudaOnCpu = CanonicalBase();
    for (size_t i = 0; i + 1 < cudaOnCpu.size(); ++i)
        if (cudaOnCpu[i] == "--backend")
            cudaOnCpu[i + 1] = "cuda";
    ExpectFailure("CUDA requested from CPU build", cudaOnCpu,
                  "built without CUDA support");

    ExpectFailureForBuild("fixed CUDA path is not silently ignored",
                          cudaOnCpu, true, true,
                          "currently uses only the CPU diffraction path");

    std::vector<std::string> hugePhi = CanonicalBase();
    hugePhi.insert(hugePhi.end(), {"--phi-points", "2147483647"});
    ExpectFailure("phi endpoint allocation overflow", hugePhi,
                  "required endpoint allocation");

    if (failures != 0)
    {
        std::cerr << failures << " CLI test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "All CLI tests passed\n";
    return EXIT_SUCCESS;
}
