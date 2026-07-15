#pragma once

#include "AdaptiveConfig.h"

#include <string>
#include <vector>

class ArgPP;

enum class RunMethod
{
    PhysicalOptics,
    GeometricalOptics
};

enum class RunBackend
{
    Auto,
    Cpu,
    Cuda
};

enum class GeometryClassification
{
    Auto,
    Convex,
    Nonconvex
};

enum class ParticleSourceMode
{
    Builtin,
    File
};

enum class OrientationMode
{
    Fixed,
    EulerGrid,
    MonteCarlo,
    File,
    DiffractionGrid,
    Sobol,
    SO3Quaternion,
    SobolSeed,
    SobolRing,
    Hammersley,
    Lattice,
    LatticeGenerator,
    EulerQuadrature,
    EulerAdaptive,
    EulerConvergence,
    Adaptive,
    Auto,
    AutoFull,
    DiffractionAutoFull
};

enum class ThetaGridMode
{
    Default,
    Uniform,
    File,
    Adaptive
};

struct RunConfig
{
    RunConfig();

    static RunConfig FromCommandLine(const ArgPP &args,
                                     bool gpuDefaultBuild,
                                     bool cudaCompiled);

    RunMethod method;
    RunBackend backend;
    GeometryClassification geometry;
    ParticleSourceMode particleSource;
    OrientationMode orientation;
    ThetaGridMode thetaGrid;

    bool useGpu;
    bool useFft;
    int fftFactor;
    double fftTolerance;
    double refractiveReal;
    double refractiveImag;
    double wavelengthUm;
    int maxReflections;
    int threads;
    std::string cutoffProfile;
    AdaptiveConvergenceLimits adaptive;
    std::vector<std::string> warnings;
};

const char *RunMethodName(RunMethod method);
const char *RunBackendName(RunBackend backend);
const char *GeometryClassificationName(GeometryClassification geometry);
const char *OrientationModeName(OrientationMode orientation);
