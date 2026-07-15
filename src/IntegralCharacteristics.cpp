#include "IntegralCharacteristics.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{

double CleanRoundoff(double value, double scale)
{
    const double tolerance = std::max(1.0, std::fabs(scale)) * 1e-10;
    return std::fabs(value) < tolerance ? 0.0 : value;
}

double RelativeDifference(double first, double second, double scale)
{
    const double denominator = std::max(
        std::max(std::fabs(first), std::fabs(second)),
        std::max(std::fabs(scale), 1.0) * 1e-15);
    return std::fabs(first - second) / denominator;
}

std::string LogNameForResult(const std::string &destName)
{
    const std::string suffix = "_noshadow";
    if (destName.size() >= suffix.size()
        && destName.compare(destName.size() - suffix.size(), suffix.size(), suffix) == 0)
    {
        return destName.substr(0, destName.size() - suffix.size()) + "_log.txt";
    }
    return destName + "_log.txt";
}

double CoarseAngularIntegral(const std::vector<double> &theta,
                             const std::vector<double> &values)
{
    std::vector<size_t> rows;
    rows.reserve((theta.size() + 1) / 2 + 1);
    for (size_t row = 0; row < theta.size(); row += 2)
        rows.push_back(row);
    if (rows.back() != theta.size() - 1)
        rows.push_back(theta.size() - 1);

    const double twoPi = 2.0 * std::acos(-1.0);
    double integral = 0.0;
    for (size_t k = 0; k < rows.size(); ++k)
    {
        const size_t row = rows[k];
        const double lower = (k == 0)
            ? theta.front()
            : 0.5 * (theta[rows[k - 1]] + theta[row]);
        const double upper = (k + 1 == rows.size())
            ? theta.back()
            : 0.5 * (theta[row] + theta[rows[k + 1]]);
        integral += values[row] * twoPi * (std::cos(lower) - std::cos(upper));
    }
    return integral;
}

} // namespace

const char *IntegralMethodName(IntegralMethod method)
{
    return method == IntegralMethod::PhysicalOptics ? "PO" : "GO";
}

double EstimateAngularIntegralRelativeError(
    const std::vector<double> &thetaRadians,
    const std::vector<double> &sampleValues,
    double fineIntegral)
{
    if (thetaRadians.size() != sampleValues.size()
        || thetaRadians.size() < 5
        || !std::isfinite(fineIntegral))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    for (size_t i = 0; i < thetaRadians.size(); ++i)
    {
        if (!std::isfinite(thetaRadians[i]) || !std::isfinite(sampleValues[i])
            || (i > 0 && !(thetaRadians[i] > thetaRadians[i - 1])))
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    const double coarseIntegral = CoarseAngularIntegral(thetaRadians, sampleValues);
    return RelativeDifference(fineIntegral, coarseIntegral, fineIntegral);
}

IntegralCharacteristics ComputeIntegralCharacteristics(
    IntegralMethod method,
    double projectedArea,
    double cScaComputed,
    double cExtOt,
    bool hasOpticalTheorem,
    bool hasAbsorption,
    double cScaQuadratureErrorRel)
{
    IntegralCharacteristics values;
    values.method = method;
    values.validProjectedArea = std::isfinite(projectedArea) && projectedArea > 0.0;
    values.hasOpticalTheorem = hasOpticalTheorem && std::isfinite(cExtOt);
    values.hasAbsorption = hasAbsorption;
    values.projectedArea = projectedArea;

    if (method == IntegralMethod::PhysicalOptics)
    {
        values.cExt = values.hasOpticalTheorem
            ? cExtOt : std::numeric_limits<double>::quiet_NaN();
    }
    else
    {
        // In geometrical optics extinction is the intercepted projected area.
        values.cExt = values.validProjectedArea
            ? projectedArea : std::numeric_limits<double>::quiet_NaN();
    }

    if (hasAbsorption)
    {
        values.cSca = cScaComputed;
        values.cAbs = CleanRoundoff(values.cExt - values.cSca, projectedArea);
    }
    else
    {
        // Enforce the physical nonabsorbing balance and use the independently
        // computed scattering value below only to estimate its numerical error.
        values.cSca = values.cExt;
        values.cAbs = 0.0;
    }

    if (std::isfinite(cScaQuadratureErrorRel)
        && cScaQuadratureErrorRel >= 0.0)
    {
        values.hasCScaErrorEstimate = true;
        values.cScaErrorEstimateRel = cScaQuadratureErrorRel;
    }

    if (!hasAbsorption && std::isfinite(values.cSca)
        && std::isfinite(cScaComputed))
    {
        const double balanceError = RelativeDifference(
            cScaComputed, values.cSca, values.cSca);
        values.cScaErrorEstimateRel = values.hasCScaErrorEstimate
            ? std::max(values.cScaErrorEstimateRel, balanceError)
            : balanceError;
        values.hasCScaErrorEstimate = true;
    }

    if (values.validProjectedArea)
    {
        const double inverseArea = 1.0 / projectedArea;
        values.qExt = values.cExt * inverseArea;
        values.qSca = values.cSca * inverseArea;
        values.qAbs = values.cAbs * inverseArea;
    }
    else
    {
        values.qExt = values.qSca = values.qAbs =
            std::numeric_limits<double>::quiet_NaN();
    }
    return values;
}

std::string IntegralCharacteristicsStatus(
    const IntegralCharacteristics &values,
    const std::string &label)
{
    if (!values.validProjectedArea)
        return "invalid_projected_area";
    if (label != "full")
        return "diagnostic_no_shadow";
    if (values.method == IntegralMethod::PhysicalOptics
        && !values.hasOpticalTheorem)
    {
        return "optical_theorem_unavailable";
    }
    if (values.hasCScaErrorEstimate && values.cScaErrorEstimateRel > 1e-2)
        return "check_csca_integration";
    return "ok";
}

std::string FormatIntegralCharacteristicsLog(
    const IntegralCharacteristics &values,
    const std::string &label)
{
    std::ostringstream log;
    log << std::fixed << std::setprecision(6);
    log << "\n===== INTEGRAL CHARACTERISTICS ("
        << IntegralMethodName(values.method) << "): " << label << " =====\n";
    log << "C_ext = " << values.cExt << "\n";
    log << "C_sca = " << values.cSca << "\n";
    log << "C_abs = " << values.cAbs << "\n";
    log << "Q_ext = " << values.qExt << "\n";
    log << "Q_sca = " << values.qSca << "\n";
    log << "Q_abs = " << values.qAbs << "\n";
    if (values.hasCScaErrorEstimate)
    {
        log << "C_sca "
            << (values.method == IntegralMethod::PhysicalOptics
                ? "integral" : "beam-sum closure")
            << " relative error estimate = "
            << values.cScaErrorEstimateRel * 100.0 << "%\n";
    }
    else
    {
        log << "C_sca relative error estimate = unavailable\n";
    }

    const std::string status = IntegralCharacteristicsStatus(values, label);
    log << "integral_status = " << status << "\n";
    log << "EFFICIENCY_SUMMARY "
        << "Cext=" << values.cExt << ' '
        << "Csca=" << values.cSca << ' '
        << "Cabs=" << values.cAbs << ' '
        << "Qext=" << values.qExt << ' '
        << "Qsca=" << values.qSca << ' '
        << "Qabs=" << values.qAbs << ' ';
    if (values.hasCScaErrorEstimate)
        log << "Csca_error_rel=" << values.cScaErrorEstimateRel << ' ';
    log << "integral_status=" << status << "\n";

    if (status == "check_csca_integration")
    {
        log << "WARNING: the C_sca numerical error estimate is "
            << values.cScaErrorEstimateRel * 100.0 << "%.\n"
            << "  Fix: refine the theta/phi grid, orientation sampling, and "
               "reflection depth before using the integral characteristics.\n";
    }
    if (status == "optical_theorem_unavailable")
    {
        log << "WARNING: PO extinction is unavailable because the forward "
               "optical-theorem value was not computed.\n"
            << "  Fix: use a full PO run with the shadow beam enabled.\n";
    }
    log << "Characteristics file: <result>_integrals.tsv\n";
    log << "==============================================\n";
    return log.str();
}

void WriteIntegralCharacteristicsTsv(
    const std::string &destName,
    const std::string &label,
    const IntegralCharacteristics &values)
{
    const std::string path = destName + "_integrals.tsv";
    std::ofstream out(path.c_str(), std::ios::out);
    if (!out.is_open())
        throw std::runtime_error(
            "cannot open integral-characteristics output '" + path
            + "'.\n  Fix: verify output permissions and free disk space.");

    out << "method\tlabel\tstatus"
        << "\tCext\tCsca\tCabs\tQext\tQsca\tQabs"
        << "\tCsca_error_estimate_rel\n";
    out << std::setprecision(17)
        << IntegralMethodName(values.method)
        << '\t' << label
        << '\t' << IntegralCharacteristicsStatus(values, label)
        << '\t' << values.cExt
        << '\t' << values.cSca
        << '\t' << values.cAbs
        << '\t' << values.qExt
        << '\t' << values.qSca
        << '\t' << values.qAbs
        << '\t';
    if (values.hasCScaErrorEstimate)
        out << values.cScaErrorEstimateRel;
    else
        out << "nan";
    out << '\n';
    out.close();
    if (!out)
        throw std::runtime_error(
            "failed while writing integral-characteristics output '" + path
            + "'.\n  Fix: verify free disk space and filesystem health.");
}

void AppendIntegralCharacteristicsLog(
    const std::string &destName,
    const std::string &text)
{
    const std::string path = LogNameForResult(destName);
    std::ofstream out(path.c_str(), std::ios::app);
    if (!out.is_open())
        throw std::runtime_error(
            "cannot open output log '" + path
            + "'.\n  Fix: verify output permissions and free disk space.");
    out << text;
    out.close();
    if (!out)
        throw std::runtime_error(
            "failed while writing output log '" + path
            + "'.\n  Fix: verify free disk space and filesystem health.");
}
