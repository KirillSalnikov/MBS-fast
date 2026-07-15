#pragma once

#include <string>
#include <vector>

enum class IntegralMethod
{
    PhysicalOptics,
    GeometricalOptics
};

struct IntegralCharacteristics
{
    IntegralMethod method = IntegralMethod::PhysicalOptics;
    bool validProjectedArea = false;
    bool hasOpticalTheorem = false;
    bool hasAbsorption = false;
    bool hasCScaErrorEstimate = false;

    double projectedArea = 0.0;
    double cExt = 0.0;
    double cSca = 0.0;
    double cAbs = 0.0;
    double qExt = 0.0;
    double qSca = 0.0;
    double qAbs = 0.0;
    double cScaErrorEstimateRel = 0.0;
};

IntegralCharacteristics ComputeIntegralCharacteristics(
    IntegralMethod method,
    double projectedArea,
    double cScaComputed,
    double cExtOt,
    bool hasOpticalTheorem,
    bool hasAbsorption,
    double cScaQuadratureErrorRel);

double EstimateAngularIntegralRelativeError(
    const std::vector<double> &thetaRadians,
    const std::vector<double> &sampleValues,
    double fineIntegral);

const char *IntegralMethodName(IntegralMethod method);

std::string IntegralCharacteristicsStatus(
    const IntegralCharacteristics &values,
    const std::string &label);

std::string FormatIntegralCharacteristicsLog(
    const IntegralCharacteristics &values,
    const std::string &label);

void WriteIntegralCharacteristicsTsv(
    const std::string &destName,
    const std::string &label,
    const IntegralCharacteristics &values);

void AppendIntegralCharacteristicsLog(
    const std::string &destName,
    const std::string &text);
