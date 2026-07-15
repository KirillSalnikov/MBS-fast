#include "FftStatistics.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

FftInterpolationStatistics::FftInterpolationStatistics()
    : calls(0), validationCalls(0), refinedCalls(0), lastOutputPhi(0),
      lastDirectPhi(0), lastFactor(0), worstEstimatedRelativeError(0.0)
{
}

void FftInterpolationStatistics::RecordCall(int outputPhi, int directPhi,
                                            int factor)
{
    std::lock_guard<std::mutex> lock(mutex);
    ++calls;
    lastOutputPhi = outputPhi;
    lastDirectPhi = directPhi;
    lastFactor = factor;
}

void FftInterpolationStatistics::RecordValidation(double relativeError,
                                                   bool refined)
{
    std::lock_guard<std::mutex> lock(mutex);
    ++validationCalls;
    if (refined)
        ++refinedCalls;
    worstEstimatedRelativeError = std::max(
        worstEstimatedRelativeError, relativeError);
}

std::string FormatFftInterpolationReport(
    const FftInterpolationStatistics &statistics,
    bool enabled,
    int requestedFactor,
    double tolerance)
{
    std::lock_guard<std::mutex> lock(statistics.mutex);
    std::ostringstream out;
    out << "\n===== FFT INTERPOLATION REPORT =====\n";
    out << "Status: " << (enabled ? "enabled" : "disabled") << "\n";
    out << "Requested phi compression factor: ";
    if (!enabled)
        out << "not requested\n";
    else if (requestedFactor > 0)
        out << requestedFactor << "\n";
    else
        out << "auto\n";
    out << "Validation tolerance: ";
    if (tolerance > 0.0)
        out << std::scientific << std::setprecision(3) << tolerance << "\n";
    else
        out << "disabled\n";
    out << std::defaultfloat;
    out << "FFT-path calls: " << statistics.calls << "\n";
    if (statistics.calls > 0)
    {
        out << "Last grid: direct N_phi=" << statistics.lastDirectPhi
            << ", output N_phi=" << statistics.lastOutputPhi
            << ", effective factor=" << statistics.lastFactor;
        if (statistics.lastFactor <= 1)
            out << " (direct, no interpolation)";
        out << "\n";
    }
    out << "Validated calls: " << statistics.validationCalls << "\n";
    out << "Calls refined to doubled direct N_phi: "
        << statistics.refinedCalls << "\n";
    if (statistics.validationCalls > 0)
    {
        out << std::scientific << std::setprecision(6)
            << "Worst nested-grid relative error estimate: "
            << statistics.worstEstimatedRelativeError << "\n";
    }
    out << "Estimator: coarse FFT result versus doubled direct-N_phi result; "
           "it is an a-posteriori estimate, not a strict error bound.\n";
    out << "====================================\n";
    return out.str();
}
