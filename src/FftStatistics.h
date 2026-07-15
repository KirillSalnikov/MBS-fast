#pragma once

#include <mutex>
#include <string>

struct FftInterpolationStatistics
{
    FftInterpolationStatistics();

    void RecordCall(int outputPhi, int directPhi, int factor);
    void RecordValidation(double relativeError, bool refined);

    mutable std::mutex mutex;
    long long calls;
    long long validationCalls;
    long long refinedCalls;
    int lastOutputPhi;
    int lastDirectPhi;
    int lastFactor;
    double worstEstimatedRelativeError;
};

std::string FormatFftInterpolationReport(
    const FftInterpolationStatistics &statistics,
    bool enabled,
    int requestedFactor,
    double tolerance);
