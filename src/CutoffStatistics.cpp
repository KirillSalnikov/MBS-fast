#include "CutoffStatistics.h"

#include <cmath>
#include <iomanip>
#include <sstream>

BeamCutoffStatistics::BeamCutoffStatistics()
    : prepareCalls(0), candidates(0), kept(0), rejected(0),
      rejectedJones(0), rejectedArea(0), rejectedImportance(0)
{
}

TraceCutoffStatistics::TraceCutoffStatistics()
    : evaluated(0), rejected(0), rejectedJones(0), rejectedArea(0),
      rejectedImportance(0), smallFragmentSimplifications(0),
      configuredBeamLimitHits(0), hardBeamLimitHits(0)
{
}

std::string FormatBeamCutoffReport(
    const BeamCutoffStatistics &statistics,
    const std::string &profile,
    double jonesRelative,
    double areaRelative,
    double importanceRelative)
{
    const long long candidates = statistics.candidates.load();
    const long long rejected = statistics.rejected.load();
    const bool enabled = jonesRelative > 0.0 || areaRelative > 0.0
        || importanceRelative > 0.0;
    std::ostringstream out;
    out << "\n===== OUTPUT BEAM CUTOFF REPORT =====\n";
    out << "Profile: " << profile << "\n";
    out << std::scientific << std::setprecision(3);
    out << "Relative thresholds: Jones=" << jonesRelative
        << ", area=" << areaRelative
        << ", importance=" << importanceRelative << "\n";
    out << std::defaultfloat;
    out << "Status: " << (enabled ? "enabled" : "disabled") << "\n";
    out << "Prepare calls: " << statistics.prepareCalls.load() << "\n";
    out << "Candidate output beams: " << candidates << "\n";
    out << "Kept output beams: " << statistics.kept.load() << "\n";
    out << "Rejected output beams: " << rejected;
    if (candidates > 0)
        out << " (" << std::fixed << std::setprecision(4)
            << 100.0 * rejected / candidates << "%)";
    out << "\n";
    out << "Reason counters (overlap allowed): Jones="
        << statistics.rejectedJones.load()
        << ", area=" << statistics.rejectedArea.load()
        << ", importance=" << statistics.rejectedImportance.load() << "\n";
    out << "=====================================\n";
    return out.str();
}

std::string FormatTraceCutoffReport(
    const TraceCutoffStatistics &statistics,
    const std::string &profile,
    double jonesRelative,
    double areaRelative,
    double importanceRelative,
    double areaRatio,
    int maximumBeams)
{
    const long long evaluated = statistics.evaluated.load();
    const long long rejected = statistics.rejected.load();
    const bool relativeEnabled = jonesRelative > 0.0 || areaRelative > 0.0
        || importanceRelative > 0.0;
    std::ostringstream out;
    out << "\n===== TRACE CUTOFF REPORT =====\n";
    out << "Profile: " << profile << "\n";
    out << std::scientific << std::setprecision(3);
    out << "Relative thresholds: Jones=" << jonesRelative
        << ", area=" << areaRelative
        << ", importance=" << importanceRelative << "\n";
    out << std::defaultfloat;
    out << "Relative pruning: " << (relativeEnabled ? "enabled" : "disabled") << "\n";
    if (std::isfinite(areaRatio))
        out << "Small-fragment area ratio: " << areaRatio << "\n";
    else
        out << "Small-fragment area ratio: disabled\n";
    out << "Maximum traced beams: "
        << (maximumBeams > 0 ? std::to_string(maximumBeams) : "disabled") << "\n";
    out << "Branches evaluated by relative cutoff: " << evaluated << "\n";
    out << "Branches rejected by relative cutoff: " << rejected;
    if (evaluated > 0)
        out << " (" << std::fixed << std::setprecision(4)
            << 100.0 * rejected / evaluated << "%)";
    out << "\n";
    out << "Reason counters (overlap allowed): Jones="
        << statistics.rejectedJones.load()
        << ", area=" << statistics.rejectedArea.load()
        << ", importance=" << statistics.rejectedImportance.load() << "\n";
    out << "Small fragments simplified: "
        << statistics.smallFragmentSimplifications.load() << "\n";
    out << "Configured beam-limit hits: "
        << statistics.configuredBeamLimitHits.load() << "\n";
    out << "Hard tree-limit hits: " << statistics.hardBeamLimitHits.load() << "\n";
    out << "===============================\n";
    return out.str();
}
