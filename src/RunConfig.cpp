#include "RunConfig.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

#include "ArgPP.h"
#include "CliOptions.h"

namespace
{

struct SelectedMode
{
    SelectedMode(const char *keyValue, OrientationMode modeValue)
        : key(keyValue), mode(modeValue)
    {
    }

    std::string key;
    OrientationMode mode;
};

std::string JoinFlags(const std::vector<std::string> &keys)
{
    std::ostringstream out;
    for (size_t i = 0; i < keys.size(); ++i)
    {
        if (i)
            out << (i + 1 == keys.size() ? " and " : ", ");
        out << CliCanonicalFlag(keys[i]);
    }
    return out.str();
}

[[noreturn]] void Fail(const std::string &problem, const std::string &fix)
{
    throw std::runtime_error(problem + "\n  Fix: " + fix);
}

double DoubleValue(const ArgPP &args, const std::string &key, size_t index = 0)
{
    try
    {
        return args.GetDoubleValue(key, index);
    }
    catch (const std::exception &)
    {
        std::string raw = args.GetStringValue(key, index);
        Fail(CliCanonicalFlag(key) + " expects a finite number at value "
             + std::to_string(index + 1) + "; got '" + raw + "'.",
             "replace '" + raw + "' with a finite numeric value.");
    }
}

int IntValue(const ArgPP &args, const std::string &key, size_t index = 0)
{
    try
    {
        return args.GetIntValue(key, index);
    }
    catch (const std::exception &)
    {
        std::string raw = args.GetStringValue(key, index);
        Fail(CliCanonicalFlag(key) + " expects an integer at value "
             + std::to_string(index + 1) + "; got '" + raw + "'.",
             "replace '" + raw + "' with an integer without a decimal part.");
    }
}

void RequireReadableFile(const ArgPP &args, const std::string &key)
{
    if (!args.IsCatched(key))
        return;
    const std::string path = args.GetStringValue(key, 0);
    std::ifstream input(path.c_str(), std::ios::in);
    if (!input.good())
    {
        Fail("cannot read " + CliCanonicalFlag(key) + " file '" + path + "'.",
             "check that the path exists and is readable from the current directory.");
    }
}

bool ReadDataToken(std::istringstream &line, std::string &token)
{
    if (!(line >> token))
        return false;
    return token.empty() || token[0] != '#';
}

double ParseFileDouble(const std::string &token, const std::string &path,
                       int lineNumber, const std::string &what)
{
    char *end = nullptr;
    errno = 0;
    const double value = std::strtod(token.c_str(), &end);
    if (errno == ERANGE || end == token.c_str() || *end != '\0'
        || !std::isfinite(value))
    {
        Fail("invalid " + what + " in '" + path + "' at line "
             + std::to_string(lineNumber) + ": '" + token + "'.",
             "replace that value with a finite number.");
    }
    return value;
}

void RejectUnexpectedFileToken(std::istringstream &line,
                               const std::string &path, int lineNumber)
{
    std::string extra;
    if (line >> extra && (extra.empty() || extra[0] != '#'))
    {
        Fail("unexpected extra value '" + extra + "' in '" + path
             + "' at line " + std::to_string(lineNumber) + ".",
             "keep only the documented columns, followed optionally by a # comment.");
    }
}

void ValidateThetaGridFile(const ArgPP &args)
{
    RequireReadableFile(args, "tgrid");
    if (!args.IsCatched("tgrid"))
        return;
    const std::string path = args.GetStringValue("tgrid", 0);
    std::ifstream input(path.c_str());
    std::vector<double> values;
    std::string text;
    int lineNumber = 0;
    while (std::getline(input, text))
    {
        ++lineNumber;
        std::istringstream line(text);
        std::string token;
        if (!ReadDataToken(line, token))
            continue;
        const double theta = ParseFileDouble(token, path, lineNumber, "theta angle");
        if (theta < 0.0 || theta > 180.0)
            Fail("theta angle in '" + path + "' at line "
                 + std::to_string(lineNumber) + " is outside [0, 180] degrees.",
                 "correct or remove that theta value.");
        RejectUnexpectedFileToken(line, path, lineNumber);
        values.push_back(theta);
    }
    std::sort(values.begin(), values.end());
    for (size_t i = 1; i < values.size(); ++i)
    {
        if (values[i] == values[i - 1])
        {
            std::ostringstream value;
            value << values[i];
            Fail("--theta-grid-file repeats theta=" + value.str() + " degrees.",
                 "remove duplicate theta rows; every grid point must be distinct.");
        }
    }
    if (values.size() < 2)
        Fail("--theta-grid-file must contain at least two distinct theta values.",
             "add two or more angles in [0, 180], one per line.");
}

void ValidateOrientationFile(const ArgPP &args)
{
    RequireReadableFile(args, "orientfile");
    if (!args.IsCatched("orientfile"))
        return;
    const std::string path = args.GetStringValue("orientfile", 0);
    std::ifstream input(path.c_str());
    std::string text;
    int lineNumber = 0;
    int orientations = 0;
    while (std::getline(input, text))
    {
        ++lineNumber;
        std::istringstream line(text);
        std::string betaToken;
        if (!ReadDataToken(line, betaToken))
            continue;
        std::string gammaToken;
        if (!(line >> gammaToken) || (!gammaToken.empty() && gammaToken[0] == '#'))
            Fail("orientation file '" + path + "' line "
                 + std::to_string(lineNumber) + " is missing gamma.",
                 "write two finite values per data line: BETA GAMMA.");
        (void)ParseFileDouble(betaToken, path, lineNumber, "beta angle");
        (void)ParseFileDouble(gammaToken, path, lineNumber, "gamma angle");
        RejectUnexpectedFileToken(line, path, lineNumber);
        ++orientations;
    }
    if (orientations == 0)
        Fail("orientation file '" + path + "' contains no orientations.",
             "add at least one BETA GAMMA data line.");
}

void ValidatePositiveListFile(const ArgPP &args)
{
    RequireReadableFile(args, "multikeq_list");
    if (!args.IsCatched("multikeq_list"))
        return;
    const std::string path = args.GetStringValue("multikeq_list", 0);
    std::ifstream input(path.c_str());
    std::string text;
    int lineNumber = 0;
    int values = 0;
    while (std::getline(input, text))
    {
        ++lineNumber;
        std::istringstream line(text);
        std::string token;
        if (!ReadDataToken(line, token))
            continue;
        const double value = ParseFileDouble(token, path, lineNumber, "k_eq value");
        if (value <= 0.0)
            Fail("k_eq value in '" + path + "' at line "
                 + std::to_string(lineNumber) + " must be positive.",
                 "replace it with a value greater than zero.");
        RejectUnexpectedFileToken(line, path, lineNumber);
        ++values;
    }
    if (values == 0)
        Fail("--k-eq-list file '" + path + "' contains no sizes.",
             "add at least one positive k_eq value.");
}

void RequirePositiveInt(const ArgPP &args, const std::string &key, size_t index = 0)
{
    if (args.IsCatched(key) && IntValue(args, key, index) < 1)
    {
        Fail(CliCanonicalFlag(key) + " value " + std::to_string(index + 1)
             + " must be at least 1.",
             "pass a positive integer for this value.");
    }
}

void RequireNonnegativeInt(const ArgPP &args, const std::string &key, size_t index = 0)
{
    if (args.IsCatched(key) && IntValue(args, key, index) < 0)
    {
        Fail(CliCanonicalFlag(key) + " value " + std::to_string(index + 1)
             + " must not be negative.",
             "pass zero or a positive integer.");
    }
}

void RequirePositiveDouble(const ArgPP &args, const std::string &key, size_t index = 0)
{
    if (args.IsCatched(key) && DoubleValue(args, key, index) <= 0.0)
    {
        Fail(CliCanonicalFlag(key) + " value " + std::to_string(index + 1)
             + " must be greater than zero.",
             "pass a positive finite number.");
    }
}

void RequireRelativeTolerance(const ArgPP &args, const std::string &key)
{
    if (!args.IsCatched(key))
        return;
    const double value = DoubleValue(args, key, 0);
    if (!(value > 0.0 && value < 1.0))
    {
        Fail(CliCanonicalFlag(key)
             + " is a relative tolerance and must be in (0, 1); got "
             + args.GetStringValue(key, 0) + ".",
             "use a fraction such as 0.02 for 2% relative accuracy.");
    }
}

void RequireNonnegativeDouble(const ArgPP &args, const std::string &key,
                              size_t index = 0)
{
    if (args.IsCatched(key) && DoubleValue(args, key, index) < 0.0)
    {
        Fail(CliCanonicalFlag(key) + " value " + std::to_string(index + 1)
             + " must not be negative.",
             "pass zero or a positive finite number.");
    }
}

void RequireRelativeCutoff(const ArgPP &args, const std::string &key)
{
    if (!args.IsCatched(key))
        return;
    const double value = DoubleValue(args, key, 0);
    if (value < 0.0 || value > 1.0)
    {
        Fail(CliCanonicalFlag(key) + " is a relative threshold and must be in [0, 1]; got "
             + args.GetStringValue(key, 0) + ".",
             "use 0 to disable the cutoff or a fraction such as 1e-6.");
    }
}

bool IsPowerOfTwo(int value)
{
    return value > 0 && (value & (value - 1)) == 0;
}

int GreatestCommonDivisor(int a, int b)
{
    a = std::abs(a);
    b = std::abs(b);
    while (b != 0)
    {
        const int remainder = a % b;
        a = b;
        b = remainder;
    }
    return a;
}

void RequirePowerOfTwo(const ArgPP &args, const std::string &key, size_t index = 0)
{
    if (!args.IsCatched(key))
        return;
    const int value = IntValue(args, key, index);
    if (!IsPowerOfTwo(value))
    {
        Fail(CliCanonicalFlag(key) + " value " + std::to_string(index + 1)
             + " must be a power of two; got " + std::to_string(value) + ".",
             "use 64, 128, 256, 512, 1024, or another power of two.");
    }
}

void RejectTogether(const ArgPP &args, const std::string &a,
                    const std::string &b, const std::string &replacement)
{
    if (args.IsCatched(a) && args.IsCatched(b))
    {
        Fail(CliCanonicalFlag(a) + " and " + CliCanonicalFlag(b)
             + " cannot be used together.", replacement);
    }
}

void RequireProductFitsInt(const ArgPP &args, const std::string &key,
                           size_t firstIndex, size_t secondIndex,
                           const std::string &what)
{
    if (!args.IsCatched(key))
        return;
    const long long first = IntValue(args, key, firstIndex);
    const long long second = IntValue(args, key, secondIndex);
    if (first > 0 && second > 0
        && first * second > std::numeric_limits<int>::max())
    {
        Fail(what + " exceeds the supported 32-bit index range.",
             "reduce one or both counts so their product is at most "
                 + std::to_string(std::numeric_limits<int>::max()) + ".");
    }
}

void ValidateOutputTemplate(const ArgPP &args)
{
    if (!args.IsCatched("o"))
        return;

    const std::string path = args.GetStringValue("o", 0);
    if (!path.empty() && path[path.size() - 1] == '/')
        Fail("--output must include a result-name prefix after the final '/'.",
             "use a path such as results/run rather than results/.");
    const size_t finalSlash = path.find_last_of('/');
    const std::string finalComponent = finalSlash == std::string::npos
        ? path : path.substr(finalSlash + 1);
    if (finalComponent == "." || finalComponent == "..")
        Fail("--output cannot use '.' or '..' as its result-name prefix.",
             "append a result name, for example ./results/run.");
    size_t percent = 0;
    while ((percent = path.find('%', percent)) != std::string::npos)
    {
        const size_t end = path.find('_', percent + 1);
        if (end == std::string::npos)
            Fail("--output contains an unterminated '%' placeholder in '" + path + "'.",
                 "terminate it with '_', for example results/%0p_run, or remove the '%' character.");

        std::string placeholder = path.substr(percent + 1, end - percent - 1);
        size_t digitCount = 0;
        while (digitCount < placeholder.size()
               && std::isdigit(static_cast<unsigned char>(placeholder[digitCount])))
            ++digitCount;
        size_t valueIndex = 0;
        if (digitCount != 0)
        {
            const std::string rawIndex = placeholder.substr(0, digitCount);
            char *indexEnd = nullptr;
            errno = 0;
            const unsigned long parsed = std::strtoul(rawIndex.c_str(), &indexEnd, 10);
            if (errno == ERANGE || !indexEnd || *indexEnd != '\0'
                || parsed > std::numeric_limits<unsigned>::max())
                Fail("--output placeholder index '" + rawIndex + "' is too large.",
                     "use a zero-based value index that exists for the referenced option.");
            valueIndex = static_cast<size_t>(parsed);
            placeholder.erase(0, digitCount);
        }

        const CliOptionSpec *spec = FindCliOptionSpec(placeholder);
        if (placeholder.empty() || !spec)
            Fail("--output references unknown placeholder '%" + path.substr(percent + 1, end - percent - 1) + "_'.",
                 "use an option name from --help, such as %0p_, or remove the placeholder.");
        if (!args.IsCatched(spec->key))
            Fail("--output placeholder '%" + path.substr(percent + 1, end - percent - 1)
                 + "_' references " + CliCanonicalFlag(spec->key) + ", which is not present.",
                 "add that option to the command or remove its output placeholder.");
        if (spec->key == "o")
            Fail("--output cannot expand a placeholder that references --output itself.",
                 "remove the %o_/%output_ placeholder or reference another option.");
        if (valueIndex >= args.GetArgNumber(spec->key))
            Fail("--output placeholder for " + CliCanonicalFlag(spec->key)
                 + " requests value index " + std::to_string(valueIndex)
                 + ", but the option has " + std::to_string(args.GetArgNumber(spec->key))
                 + " value(s).",
                 "use an existing zero-based value index or remove the placeholder.");
        percent = end + 1;
    }
}

void ValidateBuiltinParticle(const ArgPP &args, std::vector<std::string> &warnings)
{
    const unsigned count = args.GetArgNumber("p");
    if (count == 0)
    {
        Fail("--particle requires a numeric particle type and its dimensions.",
             "for example, use --particle 1 100 70 for a 100 by 70 micrometer hexagonal column.");
    }

    const int type = IntValue(args, "p", 0);
    unsigned expectedMin = 0;
    unsigned expectedMax = 0;
    switch (type)
    {
    case 1:
    case 2:
        expectedMin = expectedMax = 3;
        break;
    case 3:
        expectedMin = 3;
        expectedMax = 4;
        break;
    case 10:
    case 12:
        expectedMin = expectedMax = 4;
        break;
    case 4:
        expectedMin = 2;
        expectedMax = 4;
        break;
    case 999:
        expectedMin = expectedMax = 2;
        break;
    default:
        Fail("unknown built-in particle type " + std::to_string(type) + ".",
             "use one of 1, 2, 3, 4, 10, 12, 999, or load geometry with --particle-file FILE.");
    }

    const bool invalidDroxtalArity = type == 4 && count != 2 && count != 4;
    if (count < expectedMin || count > expectedMax || invalidDroxtalArity)
    {
        std::ostringstream expected;
        if (type == 4)
            expected << "2 (TYPE SCALE) or legacy 4 (TYPE UNUSED UNUSED SCALE)";
        else if (expectedMin == expectedMax)
            expected << expectedMin;
        else
            expected << expectedMin << " or " << expectedMax;
        Fail("--particle type " + std::to_string(type) + " expects "
             + expected.str() + " values including TYPE; got "
             + std::to_string(count) + ".",
             "check --help for the dimensions required by this particle type.");
    }

    for (unsigned i = 1; i < count; ++i)
        (void)DoubleValue(args, "p", i);

    if (type == 1 || type == 2 || type == 3 || type == 10 || type == 12)
    {
        if (DoubleValue(args, "p", 1) <= 0.0 || DoubleValue(args, "p", 2) <= 0.0)
        {
            Fail("built-in particle length and diameter must both be positive.",
                 "pass --particle TYPE LENGTH_UM DIAMETER_UM with values greater than zero.");
        }
    }
    if (type == 3 && count == 4 && DoubleValue(args, "p", 3) <= 0.0)
        Fail("bullet-rosette cap height must be positive.",
             "omit the optional cap height or pass a value greater than zero.");
    if (type == 4)
    {
        const size_t scaleIndex = count == 2 ? 1 : 3;
        if (DoubleValue(args, "p", scaleIndex) <= 0.0)
            Fail("droxtal scale must be positive.",
                 "use --particle 4 SCALE with SCALE greater than zero.");
        if (count == 4)
            warnings.push_back(
                "legacy --particle 4 L D SCALE syntax ignores L and D; use --particle 4 SCALE.");
    }
    if (type == 10)
    {
        const double height = DoubleValue(args, "p", 1);
        const double diameter = DoubleValue(args, "p", 2);
        const double angle = DoubleValue(args, "p", 3);
        const double maxAngle = std::atan(height / diameter) * 180.0 / M_PI;
        if (!(angle > 0.0 && angle < maxAngle))
        {
            std::ostringstream limit;
            limit << maxAngle;
            Fail("concave-hexagonal cavity angle must be in (0, "
                     + limit.str() + ") degrees for the requested L/D.",
                 "reduce the angle so the top and bottom cavities remain distinct and inside the particle.");
        }
    }
    if (type == 12 && IntValue(args, "p", 3) != 2)
        Fail("built-in particle type 12 currently defines exactly two hexagonal components.",
             "pass --particle 12 LENGTH_UM DIAMETER_UM 2, or load a general aggregate with --particle-file FILE.");
    if (type == 999 && DoubleValue(args, "p", 1) <= 0.0)
        Fail("built-in aggregate scale must be positive.", "pass a positive second value to --particle 999.");
}

OrientationMode ResolveOrientation(const ArgPP &args)
{
    static const SelectedMode modes[] = {
        {"fixed", OrientationMode::Fixed},
        {"random", OrientationMode::EulerGrid},
        {"montecarlo", OrientationMode::MonteCarlo},
        {"orientfile", OrientationMode::File},
        {"oldauto", OrientationMode::DiffractionGrid},
        {"sobol", OrientationMode::Sobol},
        {"so3_quat", OrientationMode::SO3Quaternion},
        {"sobol_seed", OrientationMode::SobolSeed},
        {"sobol_ring", OrientationMode::SobolRing},
        {"hammersley", OrientationMode::Hammersley},
        {"lattice", OrientationMode::Lattice},
        {"lattice_z", OrientationMode::LatticeGenerator},
        {"euler_quad", OrientationMode::EulerQuadrature},
        {"euler_adapt", OrientationMode::EulerAdaptive},
        {"adaptive_euler", OrientationMode::EulerConvergence},
        {"adaptive", OrientationMode::Adaptive},
        {"auto", OrientationMode::Auto},
        {"autofull", OrientationMode::AutoFull},
        {"oldautofull", OrientationMode::DiffractionAutoFull}
    };

    std::vector<std::string> selected;
    OrientationMode resolved = OrientationMode::Fixed;
    for (const SelectedMode &mode : modes)
    {
        if (args.IsCatched(mode.key))
        {
            selected.push_back(mode.key);
            resolved = mode.mode;
        }
    }
    if (selected.empty())
    {
        Fail("orientation mode is missing.",
             "add exactly one mode, for example --fixed-orientation 0 0, --diffraction-grid 2, or --sobol 1024.");
    }
    if (selected.size() > 1)
    {
        Fail("orientation mode is ambiguous: " + JoinFlags(selected)
             + " were provided together.",
             "keep exactly one primary orientation mode; modifiers such as --pole and --mirror-gamma may remain.");
    }
    return resolved;
}

void ValidateOrientationValues(const ArgPP &args, OrientationMode mode)
{
    switch (mode)
    {
    case OrientationMode::Fixed:
        (void)DoubleValue(args, "fixed", 0);
        (void)DoubleValue(args, "fixed", 1);
        break;
    case OrientationMode::EulerGrid:
        RequirePositiveInt(args, "random", 0);
        RequirePositiveInt(args, "random", 1);
        RequireProductFitsInt(args, "random", 0, 1, "the Euler-grid orientation count");
        if (args.IsCatched("mirror_gamma")
            && IntValue(args, "random", 1) == std::numeric_limits<int>::max())
            Fail("--mirror-gamma cannot round an INT_MAX Euler gamma count to an even value.",
                 "reduce NGAMMA by at least one.");
        break;
    case OrientationMode::MonteCarlo:
        RequirePositiveInt(args, "montecarlo");
        if (args.IsCatched("mirror_gamma")
            && IntValue(args, "montecarlo", 0) == std::numeric_limits<int>::max())
            Fail("--mirror-gamma cannot round an INT_MAX Monte Carlo count to an even value.",
                 "reduce the orientation count by at least one.");
        break;
    case OrientationMode::File:
        ValidateOrientationFile(args);
        break;
    case OrientationMode::DiffractionGrid:
        RequirePositiveInt(args, "oldauto");
        break;
    case OrientationMode::Sobol:
        RequirePositiveInt(args, "sobol");
        RequirePowerOfTwo(args, "sobol");
        break;
    case OrientationMode::SO3Quaternion:
        RequirePositiveInt(args, "so3_quat");
        break;
    case OrientationMode::SobolSeed:
        RequirePositiveInt(args, "sobol_seed", 0);
        RequirePowerOfTwo(args, "sobol_seed", 0);
        RequireNonnegativeInt(args, "sobol_seed", 1);
        break;
    case OrientationMode::SobolRing:
        RequirePositiveInt(args, "sobol_ring", 0);
        RequirePositiveInt(args, "sobol_ring", 1);
        RequireProductFitsInt(args, "sobol_ring", 0, 1, "the Sobol-ring orientation count");
        break;
    case OrientationMode::Hammersley:
        RequirePositiveInt(args, "hammersley");
        break;
    case OrientationMode::Lattice:
        RequirePositiveInt(args, "lattice");
        break;
    case OrientationMode::LatticeGenerator:
        RequirePositiveInt(args, "lattice_z", 0);
        RequirePositiveInt(args, "lattice_z", 1);
        if (GreatestCommonDivisor(IntValue(args, "lattice_z", 0),
                                  IntValue(args, "lattice_z", 1)) != 1)
        {
            Fail("--lattice-generator requires Z to be coprime with N.",
                 "choose Z with gcd(N, Z) = 1, for example --lattice-generator 64 35.");
        }
        break;
    case OrientationMode::EulerQuadrature:
        RequirePositiveInt(args, "euler_quad", 0);
        RequirePositiveInt(args, "euler_quad", 1);
        RequireProductFitsInt(args, "euler_quad", 0, 1, "the Euler-quadrature orientation count");
        break;
    case OrientationMode::EulerAdaptive:
        RequirePositiveInt(args, "euler_adapt", 0);
        RequirePositiveInt(args, "euler_adapt", 1);
        RequireProductFitsInt(args, "euler_adapt", 0, 1, "the adaptive-Euler upper orientation count");
        break;
    case OrientationMode::EulerConvergence:
        RequireRelativeTolerance(args, "adaptive_euler");
        break;
    case OrientationMode::Adaptive:
        RequireRelativeTolerance(args, "adaptive");
        break;
    case OrientationMode::Auto:
        RequireRelativeTolerance(args, "auto");
        break;
    case OrientationMode::AutoFull:
        RequireRelativeTolerance(args, "autofull");
        break;
    case OrientationMode::DiffractionAutoFull:
        RequireRelativeTolerance(args, "oldautofull");
        break;
    }

    if (args.IsCatched("b"))
    {
        if (mode != OrientationMode::EulerGrid)
            Fail("--beta-range-deg is only valid with --euler-grid.",
                 "remove --beta-range-deg or select --euler-grid NBETA NGAMMA.");
        const double lo = DoubleValue(args, "b", 0);
        const double hi = DoubleValue(args, "b", 1);
        if (lo >= hi)
            Fail("--beta-range-deg requires MIN < MAX.", "swap or correct the beta range endpoints.");
        if (lo < 0.0 || hi > 180.0)
            Fail("--beta-range-deg must stay within [0, 180] degrees.",
                 "use a physical Euler beta interval with 0 <= MIN < MAX <= 180.");
    }
    if (args.IsCatched("g"))
    {
        if (mode != OrientationMode::EulerGrid)
            Fail("--gamma-range-deg is only valid with --euler-grid.",
                 "remove --gamma-range-deg or select --euler-grid NBETA NGAMMA.");
        const double lo = DoubleValue(args, "g", 0);
        const double hi = DoubleValue(args, "g", 1);
        if (lo >= hi)
            Fail("--gamma-range-deg requires MIN < MAX.", "swap or correct the gamma range endpoints.");
        if (hi - lo > 360.0)
            Fail("--gamma-range-deg spans more than one full period.",
                 "use an interval no wider than 360 degrees.");
    }
}

bool UsesAdaptiveOrientationSearch(OrientationMode mode)
{
    return mode == OrientationMode::EulerConvergence
        || mode == OrientationMode::Adaptive
        || mode == OrientationMode::Auto
        || mode == OrientationMode::AutoFull
        || mode == OrientationMode::DiffractionAutoFull;
}

bool UsesAutoFull(OrientationMode mode)
{
    return mode == OrientationMode::AutoFull
        || mode == OrientationMode::DiffractionAutoFull;
}

bool UsesUnifiedThetaSearch(OrientationMode mode)
{
    return mode == OrientationMode::Auto
        || mode == OrientationMode::AutoFull
        || mode == OrientationMode::DiffractionAutoFull;
}

bool UsesUnifiedPhiSearch(OrientationMode mode)
{
    return UsesUnifiedThetaSearch(mode);
}

bool UsesAdaptiveEulerSearch(OrientationMode mode)
{
    return mode == OrientationMode::EulerConvergence
        || mode == OrientationMode::DiffractionAutoFull;
}

bool SupportsSymmetryOverride(OrientationMode mode)
{
    return mode == OrientationMode::Sobol
        || mode == OrientationMode::SobolSeed
        || mode == OrientationMode::SobolRing
        || mode == OrientationMode::Hammersley
        || mode == OrientationMode::Lattice
        || mode == OrientationMode::LatticeGenerator
        || mode == OrientationMode::EulerQuadrature
        || mode == OrientationMode::EulerAdaptive
        || mode == OrientationMode::EulerConvergence
        || mode == OrientationMode::Adaptive
        || mode == OrientationMode::Auto
        || mode == OrientationMode::AutoFull
        || mode == OrientationMode::DiffractionAutoFull;
}

bool SupportsOrientationChunk(OrientationMode mode)
{
    return mode != OrientationMode::Fixed
        && mode != OrientationMode::MonteCarlo
        && mode != OrientationMode::File;
}

bool SupportsSerialDmaxScan(OrientationMode mode)
{
    return mode == OrientationMode::EulerGrid
        || mode == OrientationMode::DiffractionGrid
        || mode == OrientationMode::Sobol
        || mode == OrientationMode::Lattice
        || mode == OrientationMode::LatticeGenerator
        || mode == OrientationMode::EulerQuadrature
        || mode == OrientationMode::DiffractionAutoFull;
}

bool SupportsSerialKeqScan(OrientationMode mode)
{
    return mode == OrientationMode::EulerGrid
        || mode == OrientationMode::DiffractionGrid
        || mode == OrientationMode::DiffractionAutoFull;
}

bool SupportsAdaptiveTheta(OrientationMode mode)
{
    return mode == OrientationMode::EulerGrid
        || mode == OrientationMode::Sobol
        || mode == OrientationMode::SO3Quaternion
        || mode == OrientationMode::SobolSeed
        || mode == OrientationMode::SobolRing
        || mode == OrientationMode::Hammersley
        || mode == OrientationMode::Lattice
        || mode == OrientationMode::LatticeGenerator
        || mode == OrientationMode::EulerQuadrature
        || mode == OrientationMode::EulerAdaptive
        || mode == OrientationMode::EulerConvergence
        || UsesAdaptiveOrientationSearch(mode);
}

bool SupportsAutoPhi(OrientationMode mode)
{
    return mode == OrientationMode::Sobol
        || mode == OrientationMode::SO3Quaternion
        || mode == OrientationMode::SobolSeed
        || mode == OrientationMode::SobolRing
        || mode == OrientationMode::Hammersley
        || mode == OrientationMode::Lattice
        || mode == OrientationMode::LatticeGenerator
        || mode == OrientationMode::EulerQuadrature
        || mode == OrientationMode::EulerAdaptive
        || mode == OrientationMode::EulerConvergence
        || UsesAdaptiveOrientationSearch(mode);
}

void ValidateOrientationModifiers(const ArgPP &args, const RunConfig &config)
{
    const OrientationMode mode = config.orientation;
    const bool usesAdaptiveSettings = UsesAdaptiveOrientationSearch(mode)
        || args.IsCatched("auto_tgrid")
        || args.IsCatched("auto_phi")
        || args.IsCatched("adaptive_phi")
        || args.IsCatched("adaptive_reflections");
    if (args.IsCatched("adaptive_config") && !usesAdaptiveSettings)
        Fail("--adaptive-config does not affect the selected orientation mode.",
             "select an adaptive mode or add an individual adaptive theta, alpha, or reflection search.");
    if (args.IsCatched("maxorient") && !UsesAdaptiveOrientationSearch(mode))
        Fail("--max-orientations is only used by adaptive orientation modes.",
             "select --adaptive-orientations, --auto, --autofull, or --diffraction-autofull, or remove the limit.");
    const bool thetaSearch = args.IsCatched("auto_tgrid")
        || (UsesUnifiedThetaSearch(mode)
            && config.thetaGrid == ThetaGridMode::Default);
    const bool phiSearch = args.IsCatched("auto_phi")
        || args.IsCatched("adaptive_phi")
        || (UsesUnifiedPhiSearch(mode) && !args.IsCatched("nphi"));
    const bool eulerSearch = UsesAdaptiveEulerSearch(mode);
    const bool stablePassSearch = UsesAdaptiveOrientationSearch(mode)
        || args.IsCatched("auto_phi")
        || args.IsCatched("adaptive_phi")
        || args.IsCatched("adaptive_reflections");
    if (args.IsCatched("max_theta_points") && !thetaSearch)
        Fail("--max-theta-points was provided without an adaptive theta search.",
             "add --auto-theta-grid EPS/use a unified auto mode without an explicit theta grid, or remove this limit.");
    if (args.IsCatched("max_phi_points") && !phiSearch)
        Fail("--max-phi-points was provided without an adaptive phi search.",
             "add --adaptive-phi EPS/use a unified auto mode, or remove this limit.");
    if (args.IsCatched("max_beta_points") && !eulerSearch)
        Fail("--max-beta-points is only used by adaptive Euler searches.",
             "select --adaptive-euler-grid or --diffraction-autofull, or remove this limit.");
    if (args.IsCatched("max_gamma_points") && !eulerSearch)
        Fail("--max-gamma-points is only used by adaptive Euler searches.",
             "select --adaptive-euler-grid or --diffraction-autofull, or remove this limit.");
    if (args.IsCatched("convergence_passes") && !stablePassSearch)
        Fail("--convergence-passes was provided without a search that uses pass streaks.",
             "select an adaptive orientation/phi/reflection/Euler mode, or remove this setting.");
    if ((args.IsCatched("owen_avg") || args.IsCatched("owen_seeds"))
        && !UsesAutoFull(mode))
        Fail("Owen final-seed averaging is only used by autofull modes.",
             "select --autofull or --diffraction-autofull, or remove --owen-average/--owen-seeds.");
    if (args.IsCatched("owen_avg") && args.IsCatched("owen_seeds"))
        Fail("--owen-average and --owen-seeds define the same final averaging set.",
             "keep the generated seed count or the explicit seed list, not both.");
    if (args.IsCatched("ring_points"))
    {
        const bool used = mode == OrientationMode::DiffractionGrid
            || UsesAdaptiveOrientationSearch(mode)
            || (mode == OrientationMode::EulerGrid && args.IsCatched("auto_tgrid"));
        if (!used)
            Fail("--ring-points does not affect the selected orientation mode.",
                 "remove it or use --diffraction-grid/an adaptive orientation mode.");
    }
    if (args.IsCatched("pole"))
    {
        const bool supported = mode == OrientationMode::DiffractionGrid
            || mode == OrientationMode::EulerGrid
            || mode == OrientationMode::DiffractionAutoFull;
        if (!supported)
            Fail("--pole does not apply to the selected orientation mode.",
                 "remove it or use --diffraction-grid, --euler-grid, or --diffraction-autofull.");
    }
    if (args.IsCatched("mirror_gamma")
        && (mode == OrientationMode::Fixed
            || mode == OrientationMode::File
            || mode == OrientationMode::SO3Quaternion))
        Fail("--mirror-gamma does not reduce the selected orientation mode.",
             "remove it or select a beta/gamma orientation-average mode.");
    if (args.IsCatched("sym")
        && (config.method != RunMethod::PhysicalOptics
            || !SupportsSymmetryOverride(mode)))
        Fail("--symmetry is not consumed by the selected orientation implementation.",
             "remove it or select a Sobol, lattice, quadrature, or adaptive PO orientation mode.");
    if (args.IsCatched("auto_tgrid") && !SupportsAdaptiveTheta(mode))
        Fail("--auto-theta-grid is not implemented for the selected orientation mode.",
             "use --scattering-grid/--theta-grid-file, or select a Sobol, Euler-grid, lattice, quadrature, or auto mode.");
    if (args.IsCatched("auto_phi") && !SupportsAutoPhi(mode))
        Fail("--auto-phi is not implemented for the selected orientation mode.",
             "set --phi-points explicitly or select a Sobol, lattice, quadrature, or auto mode.");
    if (args.IsCatched("adaptive_phi") && !SupportsAutoPhi(mode))
        Fail("--adaptive-phi is not implemented for the selected orientation mode.",
             "set --phi-points explicitly or select a Sobol, lattice, quadrature, or auto mode.");
    if (args.IsCatched("auto_phi") && args.IsCatched("adaptive_phi"))
        Fail("--auto-phi and --adaptive-phi request the same search with different tolerances.",
             "keep --adaptive-phi EPS, or keep --auto-phi to use its 0.02 default.");
    if ((args.IsCatched("auto_phi") || args.IsCatched("adaptive_phi"))
        && args.IsCatched("nphi"))
        Fail("adaptive phi search and --phi-points select phi sampling in different ways.",
             "keep exactly one of them.");
    if (args.IsCatched("adaptive_reflections")
        && !SupportsAutoPhi(mode))
        Fail("--adaptive-reflections is not implemented for the selected orientation mode.",
             "select a Sobol, lattice, quadrature, or auto mode, or set --max-reflections explicitly.");
    if ((UsesAutoFull(mode) || mode == OrientationMode::Auto)
        && (args.IsCatched("auto_tgrid")
            || args.IsCatched("auto_phi")
            || args.IsCatched("adaptive_phi")
            || args.IsCatched("adaptive_reflections")))
        Fail("--auto/--autofull already include the requested individual adaptive search.",
             "remove the redundant --auto-theta-grid/--auto-phi/--adaptive-phi/--adaptive-reflections modifier.");
    if (args.IsCatched("checkpoint"))
    {
        const bool supported = mode == OrientationMode::File
            || mode == OrientationMode::DiffractionGrid
            || mode == OrientationMode::EulerGrid;
        if (config.method != RunMethod::PhysicalOptics || !supported)
            Fail("--checkpoint is not implemented for the selected orientation mode.",
                 "remove it or use PO with --orientation-file, --diffraction-grid, or --euler-grid.");
    }
    if (args.IsCatched("save_betas")
        && (config.method != RunMethod::PhysicalOptics
            || (mode != OrientationMode::EulerGrid
                && mode != OrientationMode::DiffractionGrid)))
        Fail("--save-betas is not available for the selected method/orientation mode.",
             "remove it or use PO with --diffraction-grid or --euler-grid.");
    if (args.IsCatched("coh_orient")
        && (config.method != RunMethod::PhysicalOptics
            || (mode != OrientationMode::EulerGrid
                && mode != OrientationMode::DiffractionGrid)))
        Fail("--coherent-orientations is implemented only by regular PO beta/gamma grids.",
             "use PO with --diffraction-grid or --euler-grid, or remove this legacy option.");
    if (args.IsCatched("coh_orient") && args.IsCatched("incoh"))
        Fail("--coherent-orientations conflicts with --incoherent.",
             "choose coherent cross-orientation Jones accumulation or incoherent per-beam Mueller accumulation.");
    if (args.IsCatched("chunk")
        && (config.method != RunMethod::PhysicalOptics
            || !SupportsOrientationChunk(mode)))
        Fail("--orientation-chunk is not consumed by the selected method/orientation mode.",
             "remove it or use an averaged PO regular-grid/Sobol/lattice/quadrature/adaptive mode.");
    if (args.IsCatched("filter")
        && (config.method != RunMethod::PhysicalOptics
            || mode != OrientationMode::EulerGrid))
        Fail("--backscatter-filter-deg is implemented only for the PO Euler-grid path.",
             "use --method po --euler-grid NBETA NGAMMA, or select the desired cone with --scattering-grid.");
    if (args.IsCatched("jones"))
    {
        if (config.method != RunMethod::PhysicalOptics
            || mode != OrientationMode::Fixed)
            Fail("--jones-output is implemented only for fixed-orientation PO.",
                 "use --method po --fixed-orientation BETA GAMMA, or remove Jones output.");
        if (args.IsCatched("incoh"))
            Fail("--jones-output is undefined with incoherent beam accumulation.",
                 "remove --incoherent or remove --jones-output.");
    }
    if (args.IsCatched("point"))
        Fail("--backscatter-point is no longer supported by the optimized calculation paths.",
             "use --scattering-grid 0 NPHI 1 for the exact backscatter direction, or a small cone radius for nearby angles.");
}

ThetaGridMode ValidateThetaGrid(const ArgPP &args, std::vector<std::string> &warnings)
{
    std::vector<std::string> selected;
    if (args.IsCatched("grid")) selected.push_back("grid");
    if (args.IsCatched("tgrid")) selected.push_back("tgrid");
    if (args.IsCatched("auto_tgrid")) selected.push_back("auto_tgrid");
    if (selected.size() > 1)
    {
        Fail("theta grid is ambiguous: " + JoinFlags(selected)
             + " were provided together.",
             "keep one theta source: --scattering-grid, --theta-grid-file, or --auto-theta-grid.");
    }

    if (args.IsCatched("grid"))
    {
        const unsigned count = args.GetArgNumber("grid");
        if (count != 3 && count != 4)
        {
            Fail("--scattering-grid expects either 3 values (R NPHI NTH) or 4 values (T1 T2 NPHI NTH); got "
                 + std::to_string(count) + ".",
                 "use, for example, --scattering-grid 0 180 600 180.");
        }
        if (count == 3)
        {
            const double radius = DoubleValue(args, "grid", 0);
            if (radius < 0.0 || radius > 180.0)
                Fail("backscatter cone radius must be in [0, 180] degrees.",
                     "correct the first --scattering-grid value.");
            RequirePositiveInt(args, "grid", 1);
            RequirePositiveInt(args, "grid", 2);
            if ((args.IsCatched("mirror_gamma") || args.IsCatched("fft"))
                && IntValue(args, "grid", 1) > std::numeric_limits<int>::max() - 5)
                Fail("NPHI is too large to round safely to a multiple of six.",
                     "reduce NPHI to at most "
                         + std::to_string(std::numeric_limits<int>::max() - 5) + ".");
            if ((long long)IntValue(args, "grid", 1)
                    * ((long long)IntValue(args, "grid", 2) + 1)
                > std::numeric_limits<int>::max())
                Fail("the scattering grid exceeds the supported 32-bit cell index range.",
                     "reduce NPHI or NTH so NPHI*(NTH+1) is at most "
                         + std::to_string(std::numeric_limits<int>::max()) + ".");
        }
        else
        {
            const double t1 = DoubleValue(args, "grid", 0);
            const double t2 = DoubleValue(args, "grid", 1);
            if (t1 < 0.0 || t1 > 180.0 || t2 < 0.0 || t2 > 180.0 || t1 >= t2)
                Fail("theta range must satisfy 0 <= T1 < T2 <= 180 degrees.",
                     "use distinct ordered endpoints, or use the three-value R NPHI NTH form for a backscatter cone.");
            RequirePositiveInt(args, "grid", 2);
            RequirePositiveInt(args, "grid", 3);
            if ((args.IsCatched("mirror_gamma") || args.IsCatched("fft"))
                && IntValue(args, "grid", 2) > std::numeric_limits<int>::max() - 5)
                Fail("NPHI is too large to round safely to a multiple of six.",
                     "reduce NPHI to at most "
                         + std::to_string(std::numeric_limits<int>::max() - 5) + ".");
            if ((long long)IntValue(args, "grid", 2)
                    * ((long long)IntValue(args, "grid", 3) + 1)
                > std::numeric_limits<int>::max())
                Fail("the scattering grid exceeds the supported 32-bit cell index range.",
                     "reduce NPHI or NTH so NPHI*(NTH+1) is at most "
                         + std::to_string(std::numeric_limits<int>::max()) + ".");
        }
        if (args.IsCatched("nphi"))
            warnings.push_back("--phi-points overrides NPHI from --scattering-grid.");
        return ThetaGridMode::Uniform;
    }
    if (args.IsCatched("tgrid"))
    {
        ValidateThetaGridFile(args);
        return ThetaGridMode::File;
    }
    if (args.IsCatched("auto_tgrid"))
    {
        RequireRelativeTolerance(args, "auto_tgrid");
        return ThetaGridMode::Adaptive;
    }
    return ThetaGridMode::Default;
}

void ValidateGpuDeviceList(const ArgPP &args)
{
    if (!args.IsCatched("gpu_devices"))
        return;
    const std::string value = args.GetStringValue("gpu_devices", 0);
    std::istringstream input(value);
    std::string token;
    std::set<int> devices;
    while (std::getline(input, token, ','))
    {
        if (token.empty())
            Fail("--gpu-devices contains an empty device entry in '" + value + "'.",
                 "use a comma-separated list such as --gpu-devices 0,1,2,3.");
        char *end = nullptr;
        errno = 0;
        long device = std::strtol(token.c_str(), &end, 10);
        if (errno == ERANGE || end == token.c_str() || *end != '\0'
            || device < 0 || device > std::numeric_limits<int>::max())
        {
            Fail("--gpu-devices entry '" + token + "' is not a nonnegative CUDA device index.",
                 "use a comma-separated list such as --gpu-devices 0,1.");
        }
        if (!devices.insert((int)device).second)
            Fail("--gpu-devices repeats device " + token + ".",
                 "list each CUDA device once.");
    }
    if (devices.empty())
        Fail("--gpu-devices is empty.", "provide at least one CUDA device index, for example --gpu-devices 0.");
}

void ValidateScan(const ArgPP &args, const RunConfig &config)
{
    const int maxScanSizes = 100000;
    std::vector<std::string> scans;
    if (args.IsCatched("multigrid")) scans.push_back("multigrid");
    if (args.IsCatched("multikeq")) scans.push_back("multikeq");
    if (args.IsCatched("multikeq_list")) scans.push_back("multikeq_list");
    if (scans.size() > 1)
    {
        Fail("size scan is ambiguous: " + JoinFlags(scans) + " were provided together.",
             "keep exactly one of --dmax-grid, --k-eq-grid, or --k-eq-list.");
    }

    if (args.IsCatched("multigrid"))
    {
        RequirePositiveDouble(args, "multigrid", 0);
        RequirePositiveDouble(args, "multigrid", 1);
        RequirePositiveInt(args, "multigrid", 2);
        if (IntValue(args, "multigrid", 2) < 2)
            Fail("--dmax-grid requires at least two sizes.",
                 "use --resize-dmax-um for one size, or set N >= 2.");
        if (DoubleValue(args, "multigrid", 0) >= DoubleValue(args, "multigrid", 1))
            Fail("--dmax-grid requires DMIN < DMAX.", "use two distinct ordered endpoints.");
        if (IntValue(args, "multigrid", 2) > maxScanSizes)
            Fail("--dmax-grid requests more than " + std::to_string(maxScanSizes) + " sizes.",
                 "split the scan into smaller ranges or reduce N.");
    }
    if (args.IsCatched("multikeq"))
    {
        RequirePositiveDouble(args, "multikeq", 0);
        RequirePositiveDouble(args, "multikeq", 1);
        RequirePositiveInt(args, "multikeq", 2);
        if (IntValue(args, "multikeq", 2) < 2)
            Fail("--k-eq-grid requires at least two sizes.",
                 "use --k-eq for one size, or set N >= 2.");
        if (DoubleValue(args, "multikeq", 0) >= DoubleValue(args, "multikeq", 1))
            Fail("--k-eq-grid requires KMIN < KMAX.", "use two distinct ordered endpoints.");
        if (IntValue(args, "multikeq", 2) > maxScanSizes)
            Fail("--k-eq-grid requests more than " + std::to_string(maxScanSizes) + " sizes.",
                 "split the scan into smaller ranges or reduce N.");
    }
    ValidatePositiveListFile(args);

    const bool haveScan = !scans.empty();
    const bool parallelScan = args.IsCatched("multigrid_parallel")
        || (config.method == RunMethod::GeometricalOptics && haveScan);
    if (args.IsCatched("multigrid_parallel") && !haveScan)
        Fail("--scan-jobs requires a size scan.", "add --dmax-grid, --k-eq-grid, or --k-eq-list.");
    if (args.IsCatched("multigrid_threads") && !haveScan)
        Fail("--scan-threads requires a size scan.", "add --dmax-grid, --k-eq-grid, or --k-eq-list.");
    if (args.IsCatched("multigrid_threads") && !parallelScan)
        Fail("--scan-threads applies only to a process-parallel size scan.",
             "add --scan-jobs N, or use --threads N for a supported serial scan.");
    if (args.IsCatched("multigrid_threads") && args.IsCatched("threads"))
        Fail("--threads and --scan-threads both select child OpenMP workers.",
             "keep one thread-count option; use --scan-threads for an explicit parallel scan.");
    if (args.IsCatched("gpu_devices") && !args.IsCatched("multigrid_parallel"))
        Fail("--gpu-devices is used only by parallel size scans.", "add --scan-jobs N or remove --gpu-devices.");
    if (args.IsCatched("gpu_devices") && !config.useGpu)
        Fail("--gpu-devices requires the CUDA backend.", "select --backend cuda or remove the device list.");
    if (args.IsCatched("multikeq_shared_batches"))
    {
        if (!args.IsCatched("multikeq") && !args.IsCatched("multikeq_list"))
            Fail("--shared-k-eq-batches requires a k_eq scan.", "add --k-eq-grid or --k-eq-list.");
        if (!args.IsCatched("multigrid_parallel"))
            Fail("--shared-k-eq-batches requires --scan-jobs.", "add --scan-jobs 0 for automatic parallelism.");
        if (!config.useGpu)
            Fail("--shared-k-eq-batches requires the CUDA backend.", "select --backend cuda or remove shared batching.");
        if (!SupportsSerialKeqScan(config.orientation))
            Fail("--shared-k-eq-batches cannot execute the selected orientation mode in a batch child.",
                 "use --diffraction-grid, --euler-grid, or --diffraction-autofull, or remove --shared-k-eq-batches to run singleton scan children.");
    }
    if (args.IsCatched("multikeq_batch_ratio"))
    {
        if (!args.IsCatched("multikeq_shared_batches"))
            Fail("--k-eq-batch-ratio requires --shared-k-eq-batches.", "enable shared batching or remove the ratio.");
        if (DoubleValue(args, "multikeq_batch_ratio", 0) < 1.0)
            Fail("--k-eq-batch-ratio must be at least 1.", "use a value such as 1.05.");
    }
    RequireNonnegativeInt(args, "multigrid_parallel");
    RequirePositiveInt(args, "multigrid_threads");
    if (args.IsCatched("multigrid_threads")
        && IntValue(args, "multigrid_threads", 0) > 65536)
        Fail("--scan-threads exceeds the supported child-worker limit of 65536.",
             "use the per-child logical CPU count or a smaller value.");
    ValidateGpuDeviceList(args);

    if (haveScan && (args.IsCatched("rs") || args.IsCatched("k_eq")))
        Fail("a single-particle size override cannot be combined with a size scan.",
             "remove --resize-dmax-um/--k-eq and let the selected scan define every size.");

    if (config.method == RunMethod::PhysicalOptics && haveScan
        && !args.IsCatched("multigrid_parallel"))
    {
        if (args.IsCatched("multigrid")
            && !SupportsSerialDmaxScan(config.orientation))
            Fail("--dmax-grid is not implemented by the selected serial PO orientation path.",
                 "add --scan-jobs N with --particle-file FILE, or use --diffraction-grid, --euler-grid, --sobol, --lattice, --euler-quadrature, or --diffraction-autofull.");
        if ((args.IsCatched("multikeq") || args.IsCatched("multikeq_list"))
            && !SupportsSerialKeqScan(config.orientation))
            Fail("the selected serial PO orientation path does not consume a k_eq size scan.",
                 "add --scan-jobs N, or use --diffraction-grid, --euler-grid, or --diffraction-autofull.");
    }
}

} // namespace

RunConfig::RunConfig()
    : method(RunMethod::PhysicalOptics),
      backend(RunBackend::Auto),
      geometry(GeometryClassification::Auto),
      particleSource(ParticleSourceMode::Builtin),
      orientation(OrientationMode::Fixed),
      thetaGrid(ThetaGridMode::Default),
      useGpu(false),
      useFft(false),
      refractiveReal(0.0),
      refractiveImag(0.0),
      wavelengthUm(0.0),
      maxReflections(6),
      threads(0)
{
}

RunConfig RunConfig::FromCommandLine(const ArgPP &args,
                                     bool gpuDefaultBuild,
                                     bool cudaCompiled)
{
    RunConfig config;

    std::vector<std::string> methodSelectors;
    if (args.IsCatched("method")) methodSelectors.push_back("method");
    if (args.IsCatched("po")) methodSelectors.push_back("po");
    if (args.IsCatched("go")) methodSelectors.push_back("go");
    if (methodSelectors.empty())
        Fail("calculation method is missing.", "add --method po or --method go.");
    if (methodSelectors.size() > 1)
        Fail("calculation method is ambiguous: " + JoinFlags(methodSelectors)
             + " were provided together.", "keep one method selector, preferably --method po or --method go.");
    if (args.IsCatched("method"))
    {
        const std::string value = args.GetStringValue("method", 0);
        if (value == "po" || value == "physical-optics")
            config.method = RunMethod::PhysicalOptics;
        else if (value == "go" || value == "geometrical-optics")
            config.method = RunMethod::GeometricalOptics;
        else
            Fail("unknown --method value '" + value + "'.", "use --method po or --method go.");
    }
    else
        config.method = args.IsCatched("po") ? RunMethod::PhysicalOptics
                                             : RunMethod::GeometricalOptics;

    std::vector<std::string> backendSelectors;
    if (args.IsCatched("backend")) backendSelectors.push_back("backend");
    if (args.IsCatched("gpu")) backendSelectors.push_back("gpu");
    if (args.IsCatched("cpu")) backendSelectors.push_back("cpu");
    if (backendSelectors.size() > 1)
        Fail("execution backend is ambiguous: " + JoinFlags(backendSelectors)
             + " were provided together.", "keep one backend selector, preferably --backend auto, cpu, or cuda.");
    if (args.IsCatched("backend"))
    {
        const std::string value = args.GetStringValue("backend", 0);
        if (value == "auto") config.backend = RunBackend::Auto;
        else if (value == "cpu") config.backend = RunBackend::Cpu;
        else if (value == "cuda" || value == "gpu") config.backend = RunBackend::Cuda;
        else Fail("unknown --backend value '" + value + "'.", "use --backend auto, --backend cpu, or --backend cuda.");
    }
    else if (args.IsCatched("gpu")) config.backend = RunBackend::Cuda;
    else if (args.IsCatched("cpu")) config.backend = RunBackend::Cpu;

    if (config.backend == RunBackend::Cuda && !cudaCompiled)
        Fail("the CUDA backend was requested, but this binary was built without CUDA support.",
             "run a gpu/bin/mbs_po_gpu_* binary or select --backend cpu.");
    if (config.method == RunMethod::GeometricalOptics
        && config.backend == RunBackend::Cuda)
        Fail("geometrical optics does not use the CUDA diffraction backend.",
             "select --backend cpu for GO, or select --method po to use CUDA diffraction.");
    config.useGpu = config.method == RunMethod::PhysicalOptics
        && (config.backend == RunBackend::Cuda
            || (config.backend == RunBackend::Auto && gpuDefaultBuild));
    config.useFft = args.IsCatched("fft");
    if (config.useFft && config.method != RunMethod::PhysicalOptics)
        Fail("--fft is available only for physical optics.", "select --method po or remove --fft.");
    if (config.useFft && !config.useGpu)
        Fail("--fft requires the CUDA backend.", "select --backend cuda or remove --fft.");

    std::vector<std::string> geometrySelectors;
    if (args.IsCatched("geometry")) geometrySelectors.push_back("geometry");
    if (args.IsCatched("forced_convex")) geometrySelectors.push_back("forced_convex");
    if (args.IsCatched("forced_nonconvex")) geometrySelectors.push_back("forced_nonconvex");
    if (geometrySelectors.size() > 1)
        Fail("geometry classification is ambiguous: " + JoinFlags(geometrySelectors)
             + " were provided together.", "keep one selector, preferably --geometry auto, convex, or nonconvex.");
    if (args.IsCatched("geometry"))
    {
        const std::string value = args.GetStringValue("geometry", 0);
        if (value == "auto") config.geometry = GeometryClassification::Auto;
        else if (value == "convex") config.geometry = GeometryClassification::Convex;
        else if (value == "nonconvex" || value == "non-convex") config.geometry = GeometryClassification::Nonconvex;
        else Fail("unknown --geometry value '" + value + "'.", "use --geometry auto, convex, or nonconvex.");
    }
    else if (args.IsCatched("forced_convex")) config.geometry = GeometryClassification::Convex;
    else if (args.IsCatched("forced_nonconvex")) config.geometry = GeometryClassification::Nonconvex;

    if (args.IsCatched("p") == args.IsCatched("pf"))
        Fail("particle source is not unique.", "provide exactly one of --particle TYPE ... or --particle-file FILE.");
    config.particleSource = args.IsCatched("p") ? ParticleSourceMode::Builtin
                                                : ParticleSourceMode::File;
    if (config.particleSource == ParticleSourceMode::Builtin)
        ValidateBuiltinParticle(args, config.warnings);
    else
        RequireReadableFile(args, "pf");
    if (args.IsCatched("rs") && !args.IsCatched("pf"))
        Fail("--resize-dmax-um requires a file particle.", "add --particle-file FILE or remove --resize-dmax-um.");
    RejectTogether(args, "rs", "k_eq", "choose either Dmax scaling or k_eq scaling, not both.");
    RequirePositiveDouble(args, "rs");
    RequirePositiveDouble(args, "k_eq");

    if (!args.IsCatched("ri"))
        Fail("refractive index is missing.", "add --refractive-index REAL IMAG, for example --refractive-index 1.3116 0.");
    config.refractiveReal = DoubleValue(args, "ri", 0);
    config.refractiveImag = DoubleValue(args, "ri", 1);
    if (config.refractiveReal <= 0.0)
        Fail("the real refractive-index component must be positive.", "pass --refractive-index REAL IMAG with REAL > 0.");
    if (!args.IsCatched("w"))
        Fail("wavelength is missing.", "add --wavelength-um LAMBDA, for example --wavelength-um 0.532.");
    config.wavelengthUm = DoubleValue(args, "w", 0);
    if (config.wavelengthUm <= 0.0)
        Fail("wavelength must be greater than zero.", "pass a positive value in micrometers, for example --wavelength-um 0.532.");
    if (args.IsCatched("n"))
    {
        config.maxReflections = IntValue(args, "n", 0);
        if (config.maxReflections < 0)
            Fail("--max-reflections must not be negative.", "use zero or a positive integer.");
        if (config.maxReflections > 30)
            Fail("--max-reflections exceeds the supported limit of 30.",
                 "use a value from 0 to 30; deeper paths cannot be represented by the beam location history.");
    }

    config.orientation = ResolveOrientation(args);
    ValidateOrientationValues(args, config.orientation);
    config.thetaGrid = ValidateThetaGrid(args, config.warnings);
    ValidateOrientationModifiers(args, config);

    if (config.useGpu
        && (config.orientation == OrientationMode::Fixed
            || config.orientation == OrientationMode::MonteCarlo
            || config.orientation == OrientationMode::File))
    {
        Fail(std::string(OrientationModeName(config.orientation))
                 + " currently uses only the CPU diffraction path.",
             "select --backend cpu, or use a CUDA-enabled averaged mode such as --sobol, --euler-grid, or --diffraction-grid.");
    }

    if (config.method == RunMethod::GeometricalOptics)
    {
        const bool supported = config.orientation == OrientationMode::Fixed
            || config.orientation == OrientationMode::EulerGrid
            || config.orientation == OrientationMode::MonteCarlo
            || config.orientation == OrientationMode::DiffractionGrid
            || config.orientation == OrientationMode::Sobol
            || config.orientation == OrientationMode::SobolSeed;
        if (!supported)
            Fail(std::string(OrientationModeName(config.orientation)) + " orientation mode is not implemented for GO.",
                 "use --fixed-orientation, --euler-grid, --monte-carlo, --diffraction-grid, --sobol, or --sobol-seed.");
        if (args.IsCatched("pole"))
            Fail("--pole is not used by the GO implementation.", "remove --pole or select --method po.");
        if (args.IsCatched("auto_tgrid") || args.IsCatched("auto_phi")
            || args.IsCatched("adaptive_phi")
            || args.IsCatched("adaptive_reflections"))
            Fail("automatic scattering-grid selection is not implemented for GO.",
                 "use --scattering-grid/--theta-grid-file and --phi-points explicitly, or select --method po.");
        if (args.IsCatched("incoh"))
            Fail("--incoherent is not consumed by the GO implementation.",
                 "remove --incoherent, or select --method po for incoherent per-beam Mueller accumulation.");
        if (args.IsCatched("beam_cutoff") || args.IsCatched("beam_cutoff_j")
            || args.IsCatched("beam_cutoff_area")
            || args.IsCatched("beam_cutoff_importance"))
            Fail("output-beam cutoffs are implemented only by PO diffraction.",
                 "remove the --beam-cutoff* options, or select --method po; GO tracing cutoffs use --trace-cutoff*.");
    }

    if (args.IsCatched("point"))
        Fail("--backscatter-point is not implemented in the optimized solver.",
             "remove it and use --scattering-grid R NPHI NTH with a small backscatter-cone radius R.");
    if (args.IsCatched("karczewski") && config.method != RunMethod::PhysicalOptics)
        Fail("--karczewski is available only for PO.", "select --method po or remove this option.");
    if (args.IsCatched("jones") && config.method != RunMethod::PhysicalOptics)
        Fail("--jones-output is available only for PO.", "select --method po or remove this option.");
    if (args.IsCatched("noshadow_output") && config.method != RunMethod::PhysicalOptics)
        Fail("--no-shadow-output is available only for PO.", "select --method po or remove this option.");
    if (args.IsCatched("shadow_off") && config.method != RunMethod::PhysicalOptics)
        Fail("--no-shadow-beam is available only for PO.", "select --method po or remove this option.");
    if ((args.IsCatched("ot_phase_avg") || args.IsCatched("ot_phase_shift")
         || args.IsCatched("ot_ping"))
        && config.method != RunMethod::PhysicalOptics)
        Fail("optical-theorem phase diagnostics are available only for PO.",
             "select --method po or remove the --ot-* options.");
    if (args.IsCatched("gpu_trace") && !cudaCompiled)
        Fail("--gpu-trace-prefilter requires a CUDA-enabled binary.",
             "run a gpu/bin/mbs_po_gpu_* binary or remove this experimental option.");

    RejectTogether(args, "trace_prefilter", "no_trace_prefilter",
                   "keep one prefilter setting or remove both to use the default.");
    RejectTogether(args, "full_only", "noshadow_output",
                   "remove --full-only to request both full and no-shadow output.");
    RejectTogether(args, "shadow_off", "noshadow_output",
                   "keep --no-shadow-beam to calculate only without the shadow beam, or keep --no-shadow-output to write both variants.");

    RequirePositiveInt(args, "threads");
    if (args.IsCatched("threads") && IntValue(args, "threads", 0) > 65536)
        Fail("--threads exceeds the supported host-worker limit of 65536.",
             "use the number of available logical CPUs, or a smaller explicit worker count.");
    config.threads = args.IsCatched("threads") ? IntValue(args, "threads", 0) : 0;
    RequirePositiveInt(args, "ring_points");
    RequirePositiveInt(args, "nphi");
    if (args.IsCatched("nphi")
        && IntValue(args, "nphi", 0) >= std::numeric_limits<int>::max())
        Fail("--phi-points leaves no room for the required endpoint allocation.",
             "use a value below 2147483647 and size the grid to available memory.");
    if (args.IsCatched("nphi")
        && (args.IsCatched("mirror_gamma") || args.IsCatched("fft"))
        && IntValue(args, "nphi", 0) > std::numeric_limits<int>::max() - 5)
        Fail("--phi-points is too large to round safely to a multiple of six.",
             "reduce it to at most "
                 + std::to_string(std::numeric_limits<int>::max() - 5) + ".");
    RequirePositiveInt(args, "chunk");
    RequirePositiveInt(args, "owen_avg");
    RequirePositiveInt(args, "maxorient");
    RequirePowerOfTwo(args, "maxorient");
    RequireRelativeTolerance(args, "adaptive_phi");
    RequireRelativeTolerance(args, "adaptive_reflections");
    RequirePositiveInt(args, "max_theta_points");
    RequirePositiveInt(args, "max_phi_points");
    RequirePositiveInt(args, "max_beta_points");
    RequirePositiveInt(args, "max_gamma_points");
    RequirePositiveInt(args, "convergence_passes");
    if (args.IsCatched("max_theta_points")
        && IntValue(args, "max_theta_points", 0) < 17)
        Fail("--max-theta-points must be at least 17.",
             "use 17 or a larger safety limit; the default is 4097.");
    if (args.IsCatched("max_phi_points")
        && IntValue(args, "max_phi_points", 0) < 12)
        Fail("--max-phi-points must be at least 12.",
             "use 12 or a larger safety limit; the default is 2400.");
    if (args.IsCatched("max_phi_points")
        && IntValue(args, "max_phi_points", 0)
            > (std::numeric_limits<int>::max() - 6) / 2)
        Fail("--max-phi-points exceeds the safe adaptive candidate range.",
             "use at most 1073741820; the practical default is 2400.");
    if (args.IsCatched("max_beta_points")
        && IntValue(args, "max_beta_points", 0) < 2)
        Fail("--max-beta-points must be at least 2.",
             "use 2 or a larger safety limit; the default is 1024.");
    if (args.IsCatched("max_gamma_points")
        && IntValue(args, "max_gamma_points", 0) < 6)
        Fail("--max-gamma-points must be at least 6.",
             "use 6 or a larger safety limit; the default is 2400.");
    if (args.IsCatched("convergence_passes")
        && IntValue(args, "convergence_passes", 0) > 8)
        Fail("--convergence-passes above 8 is not supported.",
             "use 1 to 8 consecutive passing refinements; the default is 2.");

    const bool usesAdaptiveSettings = UsesAdaptiveOrientationSearch(config.orientation)
        || args.IsCatched("auto_tgrid")
        || args.IsCatched("auto_phi")
        || args.IsCatched("adaptive_phi")
        || args.IsCatched("adaptive_reflections");
    if (usesAdaptiveSettings)
    {
        if (args.IsCatched("adaptive_config"))
            config.adaptive = LoadAdaptiveConvergenceConfig(
                args.GetStringValue("adaptive_config", 0));
        const bool tunesReflections = UsesAutoFull(config.orientation)
            || args.IsCatched("adaptive_reflections");
        if (args.IsCatched("n") && tunesReflections)
            config.adaptive.maxReflections = config.maxReflections;
        if (args.IsCatched("max_theta_points"))
            config.adaptive.maxThetaPoints = IntValue(
                args, "max_theta_points", 0);
        if (args.IsCatched("max_phi_points"))
        {
            const int requested = IntValue(args, "max_phi_points", 0);
            config.adaptive.maxPhiPoints = (requested / 6) * 6;
            if (config.adaptive.maxPhiPoints != requested)
                config.warnings.push_back(
                    "--max-phi-points was rounded down to "
                    + std::to_string(config.adaptive.maxPhiPoints)
                    + " because adaptive alpha grids use multiples of six.");
        }
        if (args.IsCatched("max_beta_points"))
            config.adaptive.maxBetaPoints = IntValue(
                args, "max_beta_points", 0);
        if (args.IsCatched("max_gamma_points"))
        {
            const int requested = IntValue(args, "max_gamma_points", 0);
            config.adaptive.maxGammaPoints = (requested / 6) * 6;
            if (config.adaptive.maxGammaPoints != requested)
                config.warnings.push_back(
                    "--max-gamma-points was rounded down to "
                    + std::to_string(config.adaptive.maxGammaPoints)
                    + " because periodic gamma grids use multiples of six.");
        }
        if (args.IsCatched("maxorient"))
            config.adaptive.maxOrientations = IntValue(args, "maxorient", 0);
        if (args.IsCatched("convergence_passes"))
            config.adaptive.stablePasses = IntValue(
                args, "convergence_passes", 0);
        ValidateAdaptiveConvergenceLimits(
            config.adaptive,
            args.IsCatched("adaptive_config")
                ? "effective adaptive settings after command-line overrides"
                : "adaptive command-line settings");
    }
    RequireNonnegativeInt(args, "trace_max_beams");
    RequireNonnegativeInt(args, "log");
    RequirePositiveDouble(args, "r");
    RequireNonnegativeDouble(args, "trace_prefilter_margin");
    if (args.IsCatched("ot_phase_shift"))
        (void)DoubleValue(args, "ot_phase_shift", 0);
    RequireNonnegativeDouble(args, "ot_ping");

    if (args.IsCatched("sym"))
    {
        RequirePositiveInt(args, "sym", 0);
        RequirePositiveInt(args, "sym", 1);
    }
    if (args.IsCatched("filter"))
    {
        const double angle = DoubleValue(args, "filter", 0);
        if (angle < 0.0 || angle > 180.0)
            Fail("--backscatter-filter-deg must be in [0, 180].", "pass an angle in degrees within that range.");
    }
    if (args.IsCatched("abs_points"))
    {
        if (!args.IsCatched("abs") && config.refractiveImag == 0.0)
            Fail("--abs-points was provided while absorption is disabled.",
                 "add --absorption or a nonzero imaginary refractive index, or remove --abs-points.");
        const std::string value = args.GetStringValue("abs_points", 0);
        if (value != "all" && IntValue(args, "abs_points", 0) < 1)
            Fail("--abs-points must be 'all' or a positive integer.", "use --abs-points all or --abs-points 1.");
    }
    if (args.IsCatched("owen_seeds"))
    {
        std::set<int> seeds;
        for (unsigned i = 0; i < args.GetArgNumber("owen_seeds"); ++i)
        {
            const int seed = IntValue(args, "owen_seeds", i);
            if (seed < 0)
                Fail("--owen-seeds values must not be negative.", "replace negative seeds with nonnegative integers.");
            if (!seeds.insert(seed).second)
                Fail("--owen-seeds repeats seed " + std::to_string(seed) + ".",
                     "list each Owen seed once so every final sample is independent.");
        }
    }

    const char *cutoffKeys[] = {
        "beam_cutoff", "beam_cutoff_j", "beam_cutoff_area",
        "beam_cutoff_importance", "trace_cutoff", "trace_cutoff_j",
        "trace_cutoff_area", "trace_cutoff_importance"
    };
    for (const char *key : cutoffKeys)
        RequireRelativeCutoff(args, key);

    if (args.IsCatched("all") && !args.IsCatched("tr"))
        Fail("--all-trajectories requires --trajectories FILE.", "add a trajectory file or remove --all-trajectories.");
    if (args.IsCatched("gr") && !args.IsCatched("tr"))
        Fail("--trajectory-groups requires --trajectories FILE.", "add a trajectory file or remove --trajectory-groups.");
    RequireReadableFile(args, "tr");

    if (args.IsCatched("o") && args.GetStringValue("o", 0).empty())
        Fail("--output path is empty.", "provide a nonempty output path.");
    ValidateOutputTemplate(args);

    ValidateScan(args, config);
    const bool haveScan = args.IsCatched("multigrid")
        || args.IsCatched("multikeq") || args.IsCatched("multikeq_list");
    const bool parallelScan = args.IsCatched("multigrid_parallel")
        || (config.method == RunMethod::GeometricalOptics && haveScan);
    if (parallelScan && args.IsCatched("multigrid")
        && config.particleSource != ParticleSourceMode::File)
        Fail("parallel --dmax-grid currently requires a file particle.",
             "use --particle-file FILE so each child can apply --resize-dmax-um, or run a supported serial scan.");
    if (args.IsCatched("shadow"))
        config.warnings.push_back("--shadow is deprecated and has no effect; remove it from the command.");
    if (args.IsCatched("full_only"))
        config.warnings.push_back("--full-only is redundant because full-only output is already the default.");
    if ((config.orientation == OrientationMode::Auto
         || config.orientation == OrientationMode::AutoFull
         || config.orientation == OrientationMode::DiffractionAutoFull)
        && config.thetaGrid != ThetaGridMode::Default)
        config.warnings.push_back("the explicit theta grid overrides theta-grid selection from the auto orientation mode.");
    if (UsesUnifiedPhiSearch(config.orientation) && args.IsCatched("nphi"))
        config.warnings.push_back("--phi-points fixes N_phi and overrides phi selection from the auto orientation mode.");

    return config;
}

const char *RunMethodName(RunMethod method)
{
    return method == RunMethod::PhysicalOptics ? "physical optics" : "geometrical optics";
}

const char *RunBackendName(RunBackend backend)
{
    switch (backend)
    {
    case RunBackend::Auto: return "auto";
    case RunBackend::Cpu: return "cpu";
    case RunBackend::Cuda: return "cuda";
    }
    return "unknown";
}

const char *GeometryClassificationName(GeometryClassification geometry)
{
    switch (geometry)
    {
    case GeometryClassification::Auto: return "auto";
    case GeometryClassification::Convex: return "convex";
    case GeometryClassification::Nonconvex: return "nonconvex";
    }
    return "unknown";
}

const char *OrientationModeName(OrientationMode orientation)
{
    switch (orientation)
    {
    case OrientationMode::Fixed: return "fixed";
    case OrientationMode::EulerGrid: return "Euler grid";
    case OrientationMode::MonteCarlo: return "Monte Carlo";
    case OrientationMode::File: return "orientation file";
    case OrientationMode::DiffractionGrid: return "diffraction grid";
    case OrientationMode::Sobol: return "Sobol";
    case OrientationMode::SO3Quaternion: return "SO(3) quaternion";
    case OrientationMode::SobolSeed: return "seeded Sobol";
    case OrientationMode::SobolRing: return "Sobol ring";
    case OrientationMode::Hammersley: return "Hammersley";
    case OrientationMode::Lattice: return "lattice";
    case OrientationMode::LatticeGenerator: return "explicit-generator lattice";
    case OrientationMode::EulerQuadrature: return "Euler quadrature";
    case OrientationMode::EulerAdaptive: return "adaptive Euler";
    case OrientationMode::EulerConvergence: return "converged Euler grid";
    case OrientationMode::Adaptive: return "adaptive orientations";
    case OrientationMode::Auto: return "auto";
    case OrientationMode::AutoFull: return "autofull";
    case OrientationMode::DiffractionAutoFull: return "diffraction autofull";
    }
    return "unknown";
}
