#pragma once

#include <string>

struct AdaptiveConvergenceLimits
{
    AdaptiveConvergenceLimits();

    int minThetaPoints;
    int maxThetaPoints;
    int minPhiPoints;
    int maxPhiPoints;
    int minBetaPoints;
    int maxBetaPoints;
    int minGammaPoints;
    int maxGammaPoints;
    int minReflections;
    int maxReflections;
    int minOrientations;
    int maxOrientations;

    double thetaTolerance;
    double phiTolerance;
    double betaTolerance;
    double gammaTolerance;
    double reflectionTolerance;
    double orientationTolerance;

    int stablePasses;
    int minPilotOrientations;
    int maxPilotOrientations;
    int maxJointSweeps;
    int maxEulerSweeps;

    bool loadedFromFile;
    std::string sourcePath;
};

AdaptiveConvergenceLimits LoadAdaptiveConvergenceConfig(
    const std::string &path);
void ValidateAdaptiveConvergenceLimits(
    const AdaptiveConvergenceLimits &limits,
    const std::string &sourceDescription);
std::string DescribeAdaptiveConvergenceLimits(
    const AdaptiveConvergenceLimits &limits);
double ResolveAdaptiveTolerance(double configured, double fallback);

