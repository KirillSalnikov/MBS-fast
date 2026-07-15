#include "AdaptiveConfig.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

typedef int AdaptiveConvergenceLimits::* IntMember;
typedef double AdaptiveConvergenceLimits::* DoubleMember;

std::string Trim(const std::string &value)
{
    const std::string whitespace = " \t\r\n";
    const size_t first = value.find_first_not_of(whitespace);
    if (first == std::string::npos)
        return std::string();
    const size_t last = value.find_last_not_of(whitespace);
    return value.substr(first, last - first + 1);
}

[[noreturn]] void ConfigError(const std::string &problem,
                              const std::string &fix)
{
    throw std::runtime_error(problem + "\n  Fix: " + fix);
}

const std::map<std::string, IntMember> &IntFields()
{
    static const std::map<std::string, IntMember> fields = {
        {"reflections.min", &AdaptiveConvergenceLimits::minReflections},
        {"reflections.max", &AdaptiveConvergenceLimits::maxReflections},
        {"alpha.min_points", &AdaptiveConvergenceLimits::minPhiPoints},
        {"alpha.max_points", &AdaptiveConvergenceLimits::maxPhiPoints},
        {"theta.min_points", &AdaptiveConvergenceLimits::minThetaPoints},
        {"theta.max_points", &AdaptiveConvergenceLimits::maxThetaPoints},
        {"orientations.min", &AdaptiveConvergenceLimits::minOrientations},
        {"orientations.max", &AdaptiveConvergenceLimits::maxOrientations},
        {"beta.min_points", &AdaptiveConvergenceLimits::minBetaPoints},
        {"beta.max_points", &AdaptiveConvergenceLimits::maxBetaPoints},
        {"gamma.min_points", &AdaptiveConvergenceLimits::minGammaPoints},
        {"gamma.max_points", &AdaptiveConvergenceLimits::maxGammaPoints},
        {"controller.stable_passes", &AdaptiveConvergenceLimits::stablePasses},
        {"controller.min_pilot_orientations", &AdaptiveConvergenceLimits::minPilotOrientations},
        {"controller.max_pilot_orientations", &AdaptiveConvergenceLimits::maxPilotOrientations},
        {"controller.max_joint_sweeps", &AdaptiveConvergenceLimits::maxJointSweeps},
        {"controller.max_euler_sweeps", &AdaptiveConvergenceLimits::maxEulerSweeps}
    };
    return fields;
}

const std::map<std::string, DoubleMember> &DoubleFields()
{
    static const std::map<std::string, DoubleMember> fields = {
        {"reflections.tolerance", &AdaptiveConvergenceLimits::reflectionTolerance},
        {"alpha.tolerance", &AdaptiveConvergenceLimits::phiTolerance},
        {"theta.tolerance", &AdaptiveConvergenceLimits::thetaTolerance},
        {"orientations.tolerance", &AdaptiveConvergenceLimits::orientationTolerance},
        {"beta.tolerance", &AdaptiveConvergenceLimits::betaTolerance},
        {"gamma.tolerance", &AdaptiveConvergenceLimits::gammaTolerance}
    };
    return fields;
}

int ParseInteger(const std::string &raw, const std::string &path,
                 int lineNumber, const std::string &key)
{
    char *end = nullptr;
    errno = 0;
    const long value = std::strtol(raw.c_str(), &end, 10);
    if (errno == ERANGE || end == raw.c_str() || *end != '\0'
        || value < std::numeric_limits<int>::min()
        || value > std::numeric_limits<int>::max())
    {
        ConfigError("adaptive config '" + path + "' line "
                    + std::to_string(lineNumber) + ": '" + key
                    + "' expects an integer, got '" + raw + "'.",
                    "replace the value with an integer without units or a decimal point.");
    }
    return static_cast<int>(value);
}

double ParseTolerance(const std::string &raw, const std::string &path,
                      int lineNumber, const std::string &key)
{
    char *end = nullptr;
    errno = 0;
    const double value = std::strtod(raw.c_str(), &end);
    if (errno == ERANGE || end == raw.c_str() || *end != '\0'
        || !std::isfinite(value))
    {
        ConfigError("adaptive config '" + path + "' line "
                    + std::to_string(lineNumber) + ": '" + key
                    + "' expects a finite number, got '" + raw + "'.",
                    "write a relative tolerance such as 0.02 for 2 percent.");
    }
    if (!(value > 0.0 && value < 1.0))
    {
        ConfigError("adaptive config '" + path + "' line "
                    + std::to_string(lineNumber) + ": '" + key
                    + "' must be in (0, 1), got '" + raw + "'.",
                    "write a fraction, for example 0.01 for 1 percent.");
    }
    return value;
}

int EditDistance(const std::string &left, const std::string &right)
{
    std::vector<int> previous(right.size() + 1);
    std::vector<int> current(right.size() + 1);
    for (size_t j = 0; j <= right.size(); ++j)
        previous[j] = static_cast<int>(j);
    for (size_t i = 1; i <= left.size(); ++i)
    {
        current[0] = static_cast<int>(i);
        for (size_t j = 1; j <= right.size(); ++j)
        {
            const int substitution = previous[j - 1]
                + (left[i - 1] == right[j - 1] ? 0 : 1);
            current[j] = std::min(std::min(previous[j] + 1,
                                           current[j - 1] + 1),
                                  substitution);
        }
        previous.swap(current);
    }
    return previous[right.size()];
}

std::string ClosestKey(const std::string &key)
{
    std::string closest;
    int best = std::numeric_limits<int>::max();
    for (const auto &field : IntFields())
    {
        const int distance = EditDistance(key, field.first);
        if (distance < best)
        {
            best = distance;
            closest = field.first;
        }
    }
    for (const auto &field : DoubleFields())
    {
        const int distance = EditDistance(key, field.first);
        if (distance < best)
        {
            best = distance;
            closest = field.first;
        }
    }
    return best <= 5 ? closest : std::string();
}

bool IsPowerOfTwo(int value)
{
    return value > 0 && (value & (value - 1)) == 0;
}

void RequireRange(const std::string &name, int minimum, int maximum,
                  int hardMinimum, const std::string &source)
{
    if (minimum < hardMinimum)
        ConfigError(source + ": " + name + " minimum "
                    + std::to_string(minimum) + " is below "
                    + std::to_string(hardMinimum) + ".",
                    "increase " + name + " minimum in the adaptive config.");
    if (maximum <= minimum)
        ConfigError(source + ": " + name + " range "
                    + std::to_string(minimum) + ".."
                    + std::to_string(maximum) + " cannot be refined.",
                    "set the maximum above the minimum so at least one convergence comparison is possible.");
}

std::string ToleranceText(double value)
{
    if (value <= 0.0)
        return "mode EPS";
    std::ostringstream out;
    out << value << " (" << value * 100.0 << "%)";
    return out.str();
}

} // namespace

AdaptiveConvergenceLimits::AdaptiveConvergenceLimits()
    : minThetaPoints(33),
      maxThetaPoints(4097),
      minPhiPoints(12),
      maxPhiPoints(2400),
      minBetaPoints(2),
      maxBetaPoints(1024),
      minGammaPoints(6),
      maxGammaPoints(2400),
      minReflections(2),
      maxReflections(30),
      minOrientations(64),
      maxOrientations(0),
      thetaTolerance(0.0),
      phiTolerance(0.0),
      betaTolerance(0.0),
      gammaTolerance(0.0),
      reflectionTolerance(0.0),
      orientationTolerance(0.0),
      stablePasses(2),
      minPilotOrientations(32),
      maxPilotOrientations(2048),
      maxJointSweeps(4),
      maxEulerSweeps(3),
      loadedFromFile(false)
{
}

AdaptiveConvergenceLimits LoadAdaptiveConvergenceConfig(
    const std::string &path)
{
    std::ifstream input(path.c_str());
    if (!input.good())
    {
        ConfigError("cannot read adaptive config file '" + path + "'.",
                    "check --adaptive-config path and file permissions.");
    }

    AdaptiveConvergenceLimits limits;
    limits.loadedFromFile = true;
    limits.sourcePath = path;
    std::set<std::string> assigned;
    std::string text;
    int lineNumber = 0;
    while (std::getline(input, text))
    {
        ++lineNumber;
        if (lineNumber == 1 && text.size() >= 3
            && static_cast<unsigned char>(text[0]) == 0xef
            && static_cast<unsigned char>(text[1]) == 0xbb
            && static_cast<unsigned char>(text[2]) == 0xbf)
        {
            text.erase(0, 3);
        }
        const size_t comment = text.find_first_of("#;");
        if (comment != std::string::npos)
            text.erase(comment);
        text = Trim(text);
        if (text.empty())
            continue;

        const size_t equals = text.find('=');
        if (equals == std::string::npos || text.find('=', equals + 1) != std::string::npos)
        {
            ConfigError("adaptive config '" + path + "' line "
                        + std::to_string(lineNumber)
                        + " must have exactly one '='.",
                        "use the form key = value, followed optionally by a # comment.");
        }
        const std::string key = Trim(text.substr(0, equals));
        const std::string value = Trim(text.substr(equals + 1));
        if (key.empty() || value.empty())
        {
            ConfigError("adaptive config '" + path + "' line "
                        + std::to_string(lineNumber)
                        + " has an empty key or value.",
                        "use the form key = value.");
        }
        if (!assigned.insert(key).second)
        {
            ConfigError("adaptive config '" + path + "' line "
                        + std::to_string(lineNumber) + " repeats key '"
                        + key + "'.",
                        "keep one assignment for each key.");
        }

        const auto intField = IntFields().find(key);
        if (intField != IntFields().end())
        {
            limits.*(intField->second) = ParseInteger(
                value, path, lineNumber, key);
            continue;
        }
        const auto doubleField = DoubleFields().find(key);
        if (doubleField != DoubleFields().end())
        {
            limits.*(doubleField->second) = ParseTolerance(
                value, path, lineNumber, key);
            continue;
        }

        const std::string closest = ClosestKey(key);
        ConfigError("adaptive config '" + path + "' line "
                    + std::to_string(lineNumber) + " has unknown key '"
                    + key + "'.",
                    closest.empty()
                        ? "use only keys documented in configs/adaptive.example.conf."
                        : "replace it with '" + closest + "'.");
    }

    if (assigned.empty())
    {
        ConfigError("adaptive config file '" + path + "' contains no settings.",
                    "start from configs/adaptive.example.conf and keep the settings you need.");
    }
    ValidateAdaptiveConvergenceLimits(limits,
        "adaptive config '" + path + "'");
    return limits;
}

void ValidateAdaptiveConvergenceLimits(
    const AdaptiveConvergenceLimits &limits,
    const std::string &sourceDescription)
{
    const std::pair<const char *, double> tolerances[] = {
        {"reflections.tolerance", limits.reflectionTolerance},
        {"alpha.tolerance", limits.phiTolerance},
        {"theta.tolerance", limits.thetaTolerance},
        {"orientations.tolerance", limits.orientationTolerance},
        {"beta.tolerance", limits.betaTolerance},
        {"gamma.tolerance", limits.gammaTolerance}
    };
    for (const auto &entry : tolerances)
    {
        if (entry.second != 0.0
            && !(entry.second > 0.0 && entry.second < 1.0))
        {
            ConfigError(sourceDescription + ": " + entry.first
                        + " must be in (0, 1).",
                        "use a fraction such as 0.01 for 1 percent.");
        }
    }

    RequireRange("reflections", limits.minReflections,
                 limits.maxReflections, 1, sourceDescription);
    if (limits.maxReflections > 30)
        ConfigError(sourceDescription
                    + ": reflections.max exceeds the supported limit of 30.",
                    "use reflections.max = 30 or a smaller value.");
    RequireRange("alpha points", limits.minPhiPoints,
                 limits.maxPhiPoints, 6, sourceDescription);
    RequireRange("theta points", limits.minThetaPoints,
                 limits.maxThetaPoints, 17, sourceDescription);
    RequireRange("beta points", limits.minBetaPoints,
                 limits.maxBetaPoints, 2, sourceDescription);
    RequireRange("gamma points", limits.minGammaPoints,
                 limits.maxGammaPoints, 6, sourceDescription);

    if (limits.minPhiPoints % 6 != 0 || limits.maxPhiPoints % 6 != 0)
        ConfigError(sourceDescription
                    + ": alpha.min_points and alpha.max_points must be multiples of 6.",
                    "use values such as 12, 18, 24, ..., 2400.");
    if (limits.minGammaPoints % 6 != 0 || limits.maxGammaPoints % 6 != 0)
        ConfigError(sourceDescription
                    + ": gamma.min_points and gamma.max_points must be multiples of 6.",
                    "use values such as 6, 12, 24, ..., 2400.");
    if ((long long)limits.maxBetaPoints * limits.maxGammaPoints
        > std::numeric_limits<int>::max())
        ConfigError(sourceDescription
                    + ": beta.max_points * gamma.max_points exceeds the supported orientation index range.",
                    "reduce one or both Euler limits so their product is at most "
                        + std::to_string(std::numeric_limits<int>::max()) + ".");

    if (limits.minOrientations < 16
        || !IsPowerOfTwo(limits.minOrientations))
        ConfigError(sourceDescription
                    + ": orientations.min must be a power of two and at least 16.",
                    "use 16, 32, 64, 128, and so on.");
    if (limits.maxOrientations < 0)
        ConfigError(sourceDescription
                    + ": orientations.max cannot be negative.",
                    "use a positive power of two, or omit the key for the physics-based cap.");
    if (limits.maxOrientations > 0
        && (!IsPowerOfTwo(limits.maxOrientations)
            || limits.maxOrientations <= limits.minOrientations))
        ConfigError(sourceDescription
                    + ": orientations.max must be a power of two above orientations.min.",
                    "increase it to the next power of two, or omit the key to use the physics-based cap.");

    if (limits.stablePasses < 1 || limits.stablePasses > 8)
        ConfigError(sourceDescription
                    + ": controller.stable_passes must be in [1, 8].",
                    "use 2 for production convergence or 1 for a quick exploratory run.");
    if (!IsPowerOfTwo(limits.minPilotOrientations)
        || limits.minPilotOrientations < 16)
        ConfigError(sourceDescription
                    + ": controller.min_pilot_orientations must be a power of two and at least 16.",
                    "use 32, 64, 128, and so on.");
    if (!IsPowerOfTwo(limits.maxPilotOrientations)
        || limits.maxPilotOrientations < limits.minPilotOrientations)
        ConfigError(sourceDescription
                    + ": controller.max_pilot_orientations must be a power of two not below the pilot minimum.",
                    "increase the maximum pilot count.");
    if (limits.maxJointSweeps < 1 || limits.maxJointSweeps > 12)
        ConfigError(sourceDescription
                    + ": controller.max_joint_sweeps must be in [1, 12].",
                    "use 4 unless the coupled search demonstrably needs more passes.");
    if (limits.maxEulerSweeps < 1 || limits.maxEulerSweeps > 12)
        ConfigError(sourceDescription
                    + ": controller.max_euler_sweeps must be in [1, 12].",
                    "use 3 unless beta/gamma coordinate refinement needs more passes.");
}

std::string DescribeAdaptiveConvergenceLimits(
    const AdaptiveConvergenceLimits &limits)
{
    std::ostringstream out;
    out << "Adaptive config: " << limits.sourcePath << '\n'
        << "  reflections n: " << limits.minReflections << ".."
        << limits.maxReflections << ", eps="
        << ToleranceText(limits.reflectionTolerance) << '\n'
        << "  alpha/phi points: " << limits.minPhiPoints << ".."
        << limits.maxPhiPoints << ", eps="
        << ToleranceText(limits.phiTolerance) << '\n'
        << "  theta points: " << limits.minThetaPoints << ".."
        << limits.maxThetaPoints << ", eps="
        << ToleranceText(limits.thetaTolerance) << '\n'
        << "  Sobol orientations: " << limits.minOrientations << "..";
    if (limits.maxOrientations > 0)
        out << limits.maxOrientations;
    else
        out << "physics cap";
    out << ", eps=" << ToleranceText(limits.orientationTolerance) << '\n'
        << "  beta points: " << limits.minBetaPoints << ".."
        << limits.maxBetaPoints << ", eps="
        << ToleranceText(limits.betaTolerance) << '\n'
        << "  gamma points: " << limits.minGammaPoints << ".."
        << limits.maxGammaPoints << ", eps="
        << ToleranceText(limits.gammaTolerance) << '\n'
        << "  controller: stable passes=" << limits.stablePasses
        << ", pilot=" << limits.minPilotOrientations << ".."
        << limits.maxPilotOrientations
        << ", joint sweeps=" << limits.maxJointSweeps
        << ", Euler sweeps=" << limits.maxEulerSweeps << '\n';
    return out.str();
}

double ResolveAdaptiveTolerance(double configured, double fallback)
{
    return configured > 0.0 ? configured : fallback;
}
