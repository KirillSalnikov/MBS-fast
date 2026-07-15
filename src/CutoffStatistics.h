#pragma once

#include <atomic>
#include <string>

struct BeamCutoffStatistics
{
    BeamCutoffStatistics();

    std::atomic<long long> prepareCalls;
    std::atomic<long long> candidates;
    std::atomic<long long> kept;
    std::atomic<long long> rejected;
    std::atomic<long long> rejectedJones;
    std::atomic<long long> rejectedArea;
    std::atomic<long long> rejectedImportance;
};

struct TraceCutoffStatistics
{
    TraceCutoffStatistics();

    std::atomic<long long> evaluated;
    std::atomic<long long> rejected;
    std::atomic<long long> rejectedJones;
    std::atomic<long long> rejectedArea;
    std::atomic<long long> rejectedImportance;
    std::atomic<long long> smallFragmentSimplifications;
    std::atomic<long long> configuredBeamLimitHits;
    std::atomic<long long> hardBeamLimitHits;
};

std::string FormatBeamCutoffReport(
    const BeamCutoffStatistics &statistics,
    const std::string &profile,
    double jonesRelative,
    double areaRelative,
    double importanceRelative);

std::string FormatTraceCutoffReport(
    const TraceCutoffStatistics &statistics,
    const std::string &profile,
    double jonesRelative,
    double areaRelative,
    double importanceRelative,
    double areaRatio,
    int maximumBeams);
