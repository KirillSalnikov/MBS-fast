#include "TracerPOTotal.h"
#include "HandlerPOTotal.h"
#include "HandlerPO.h"
#include "BeamCache.h"
#include "IntegralCharacteristics.h"
#include "Sobol.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <future>
#include <algorithm>
#include <iterator>
#include <limits>
#include <atomic>
#include <exception>
#include <mutex>
#include <sys/stat.h>

static matrix MirrorGammaMuellerMatrix(const matrix &src, const matrix &mirror)
{
    static const double parity[4] = { 1.0, 1.0, -1.0, -1.0 };
    matrix out(4, 4);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            out[i][j] = 0.5 * (src[i][j] + parity[i] * parity[j] * mirror[i][j]);
    return out;
}

static long long ReadMeminfoKb(const char *key)
{
    std::ifstream f("/proc/meminfo");
    std::string name;
    long long value = 0;
    std::string unit;
    while (f >> name >> value >> unit)
    {
        if (name == key)
            return value;
    }
    return 0;
}

static long long ReadProcessRssKb()
{
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line))
    {
        if (line.find("VmRSS:") == 0)
        {
            long long kb = 0;
            std::sscanf(line.c_str(), "VmRSS: %lld", &kb);
            return kb;
        }
    }
    return 0;
}

static double EnvDouble(const char *name, double fallback)
{
    const char *value = std::getenv(name);
    if (!value || !*value)
        return fallback;
    char *end = nullptr;
    double parsed = std::strtod(value, &end);
    return (end && *end == '\0' && parsed > 0.0) ? parsed : fallback;
}

static int EnvInt(const char *name, int fallback)
{
    const char *value = std::getenv(name);
    if (!value || !*value)
        return fallback;
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    return (end && *end == '\0' && parsed > 0) ? (int)parsed : fallback;
}

static std::string SizeFileLabel(double value)
{
    std::ostringstream out;
    out << std::setprecision(12) << value;
    std::string label = out.str();
    for (char &c : label)
    {
        if (c == '.')
            c = 'p';
        else if (c == '-')
            c = 'm';
        else if (c == '+')
            c = 'p';
    }
    return label;
}

class ParallelExceptionState
{
public:
    ParallelExceptionState() : m_failed(false) {}

    bool Failed() const
    {
        return m_failed.load(std::memory_order_relaxed);
    }

    void Capture()
    {
        bool expected = false;
        if (m_failed.compare_exchange_strong(
                expected, true, std::memory_order_relaxed))
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_error = std::current_exception();
        }
    }

    void Rethrow()
    {
        std::exception_ptr error;
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            error = m_error;
        }
        if (error)
            std::rethrow_exception(error);
    }

private:
    std::atomic<bool> m_failed;
    mutable std::mutex m_mutex;
    std::exception_ptr m_error;
};

static long long HostMemoryBudgetOverrideKb()
{
    const char *value = std::getenv("MBS_HOST_MEM_BUDGET_MB");
    if (!value || !*value)
        return 0;
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    return (end && *end == '\0' && parsed > 0)
        ? (long long)parsed * 1024LL
        : 0;
}

static long long EffectiveMemAvailableMb()
{
    long long availableKb = ReadMeminfoKb("MemAvailable:");
    long long availableMb = availableKb > 0 ? availableKb / 1024 : 2048;
    long long overrideKb = HostMemoryBudgetOverrideKb();
    if (overrideKb > 0)
        availableMb = std::min(availableMb, overrideKb / 1024);
    return std::max(1LL, availableMb);
}

static long long Arr2DStorageKb(int nAz, int nZen, int rows, int cols,
                                bool complexValues)
{
    const long long nPhi = (long long)std::max(0, nAz) + 1LL;
    const long long nTheta = (long long)std::max(0, nZen) + 1LL;
    const long long cells = nPhi * nTheta;
    const long long scalarBytes = complexValues ? (long long)sizeof(::complex)
                                                : (long long)sizeof(double);
    const long long dataBytes = cells * (long long)rows * (long long)cols
        * scalarBytes;
    const long long pointerBytes = nPhi * (long long)sizeof(void **)
        + cells * (long long)sizeof(void *);
    return (dataBytes + pointerBytes + 1023LL) / 1024LL;
}

static long long EstimateOldautoTransientGridKb(int nAz, int nZen,
                                                bool computeNoShadow,
                                                bool coherent,
                                                bool gpuEnabled)
{
    const long long realGridKb = Arr2DStorageKb(nAz, nZen, 4, 4, false);
    const long long jonesGridKb = Arr2DStorageKb(nAz, nZen, 2, 2, true);
    long long transientKb = 0;

    if (gpuEnabled)
    {
        // betaM and localM are allocated after the chunk decision.  The global
        // M/M_noshadow arrays are already reflected in VmRSS at this point.
        transientKb += 2LL * realGridKb;
        if (computeNoShadow)
            transientKb += 2LL * realGridKb;
    }
    else
    {
        transientKb += realGridKb;
        if (computeNoShadow)
            transientKb += realGridKb;
        if (coherent)
        {
            transientKb += jonesGridKb;
            if (computeNoShadow)
                transientKb += jonesGridKb;
        }
    }

    const double safety = EnvDouble("MBS_OLDAUTO_GRID_MEM_SAFETY", 1.25);
    return (long long)std::ceil((double)transientKb
                                * std::max(1.0, std::min(4.0, safety)));
}

static int HostMemoryAwareGammaChunk(int requested,
                                     bool noBeamCutoff,
                                     int nAz,
                                     int nZen,
                                     bool computeNoShadow,
                                     bool coherent,
                                     bool gpuEnabled,
                                     long long *totalKbOut,
                                     long long *availableKbOut,
                                     long long *rssKbOut,
                                     long long *budgetKbOut,
                                     long long *gridKbOut)
{
    const long long totalKb = ReadMeminfoKb("MemTotal:");
    const long long availableKb = ReadMeminfoKb("MemAvailable:");
    const long long rssKb = ReadProcessRssKb();
    if (totalKbOut) *totalKbOut = totalKb;
    if (availableKbOut) *availableKbOut = availableKb;
    if (rssKbOut) *rssKbOut = rssKb;
    if (gridKbOut) *gridKbOut = 0;

    if (totalKb <= 0 || availableKb <= 0)
    {
        if (budgetKbOut) *budgetKbOut = 0;
        return std::max(1, requested);
    }

    const double fraction = EnvDouble("MBS_HOST_MEM_FRACTION",
                                      noBeamCutoff ? 0.90 : 0.55);
    const long long reserveKb =
        (long long)EnvInt("MBS_HOST_MEM_RESERVE_MB", 4096) * 1024LL;
    const long long freeBudgetKb = std::max(0LL, availableKb - reserveKb);
    const double expandableFraction = std::max(0.05, std::min(0.95, fraction));
    const long long expandableKb =
        (long long)(expandableFraction * (double)freeBudgetKb);
    const long long processCeilingKb = std::max(256LL * 1024LL,
        totalKb > reserveKb ? totalKb - reserveKb : totalKb / 2);
    long long budgetKb = std::max(256LL * 1024LL,
        std::min(processCeilingKb, rssKb + expandableKb));
    const long long overrideKb = HostMemoryBudgetOverrideKb();
    if (overrideKb > 0)
        budgetKb = std::min(budgetKb, overrideKb);
    if (budgetKbOut) *budgetKbOut = budgetKb;

    const int mbPerGamma = EnvInt("MBS_OLDAUTO_BYTES_PER_GAMMA_MB",
                                  noBeamCutoff ? 64 : 128);
    const long long perGammaKb = std::max(1, mbPerGamma) * 1024LL;
    const long long gridTransientKb = EstimateOldautoTransientGridKb(
        nAz, nZen, computeNoShadow, coherent, gpuEnabled);
    if (gridKbOut) *gridKbOut = gridTransientKb;
    const long long remainingKb = budgetKb - rssKb - gridTransientKb;
    int byMemory = (int)std::max(1LL, remainingKb / perGammaKb);
    return std::max(1, std::min(requested, byMemory));
}

static void ApplyMirrorGammaMueller(Arr2D &arr, int nAz, int nZen)
{
    Arr2D mirrored(nAz + 1, nZen + 1, 4, 4);
    mirrored.ClearArr();
    for (int p = 0; p < nAz; ++p)
    {
        int pm = (p == 0) ? 0 : (nAz - p);
        for (int t = 0; t <= nZen; ++t)
            mirrored.insert(p, t, MirrorGammaMuellerMatrix(arr(p, t), arr(pm, t)));
    }
    arr = mirrored;
}

static bool OldautoGammaStaggerEnabled()
{
    const char *env = std::getenv("MBS_OLDAUTO_GAMMA_STAGGER");
    if (!env || !*env)
        return false;
    return env[0] == '1' && env[1] == '\0';
}

static double OldautoGammaAngle(const AngleRange &gammaRange, int nGamma,
                                int gammaIndex, int betaIndex,
                                bool stagger)
{
    double unit = gammaIndex + 0.5;
    if (stagger && nGamma > 1)
    {
        // Deterministic Cranley-Patterson shift per beta ring.  This keeps the
        // midpoint rule in gamma but prevents every beta ring from sampling the
        // same periodic phase, which can alias narrow specular events.
        const double golden = 0.6180339887498948482;
        double shift = std::fmod((betaIndex + 0.5) * golden, 1.0);
        unit += shift;
        unit -= std::floor(unit / nGamma) * nGamma;
    }
    return gammaRange.min + unit * gammaRange.step;
}

static bool OldautoBetaMidpointEnabled()
{
    const char *env = std::getenv("MBS_OLDAUTO_BETA_MIDPOINT");
    if (!env || !*env)
        return true;
    return !(env[0] == '0' && env[1] == '\0');
}

static int OldautoBetaCount(const AngleRange &betaRange, bool midpoint)
{
    return midpoint ? betaRange.number : betaRange.number + 1;
}

static double OldautoBetaAngle(const AngleRange &betaRange, int betaIndex,
                               bool midpoint)
{
    return betaRange.min + (betaIndex + (midpoint ? 0.5 : 0.0))
        * betaRange.step;
}

static bool SharedFftGlobalRefineEnabled()
{
    const char *env = std::getenv("MBS_FFT_GLOBAL_REFINE");
    if (!env || !*env)
        return false;
    return env[0] == '1' && env[1] == '\0';
}

static double SharedFftGlobalRefineThreshold()
{
    const char *env = std::getenv("MBS_FFT_GLOBAL_REFINE_THRESHOLD");
    if (!env || !*env)
        return 0.04;
    char *end = nullptr;
    double value = std::strtod(env, &end);
    if (!end || *end != '\0' || value <= 0.0)
        return 0.04;
    return value;
}

static int ReadEnvInt(const char *name, int fallback, int minValue, int maxValue)
{
    const char *env = std::getenv(name);
    if (!env || !*env)
        return fallback;
    char *end = nullptr;
    long value = std::strtol(env, &end, 10);
    if (!end || *end != '\0')
        return fallback;
    if (value < minValue)
        return minValue;
    if (value > maxValue)
        return maxValue;
    return (int)value;
}

static double ReadEnvDouble(const char *name, double fallback,
                            double minValue, double maxValue)
{
    const char *env = std::getenv(name);
    if (!env || !*env)
        return fallback;
    char *end = nullptr;
    double value = std::strtod(env, &end);
    if (!end || *end != '\0' || !std::isfinite(value))
        return fallback;
    if (value < minValue)
        return minValue;
    if (value > maxValue)
        return maxValue;
    return value;
}

static int RoundUpToMultiple(int value, int step)
{
    if (step <= 1)
        return value;
    return ((value + step - 1) / step) * step;
}

static int RoundDownToMultiple(int value, int step)
{
    if (step <= 1)
        return value;
    return std::max(step, (value / step) * step);
}

static bool AutoFullThetaZonesEnabled()
{
    const char *env = std::getenv("MBS_AUTOFULL_THETA_ZONES");
    if (!env || !*env)
        return false;
    return !(env[0] == '0' && env[1] == '\0');
}

static bool AutoFullRowOrientEnabled()
{
    const char *env = std::getenv("MBS_AUTOFULL_ROW_ORIENT");
    if (!env || !*env)
        return false;
    return !(env[0] == '0' && env[1] == '\0');
}

static int AutoFullThetaZone(double thetaDeg)
{
    const double forwardMax =
        ReadEnvDouble("MBS_AUTOFULL_ZONE_FORWARD_DEG", 20.0, 0.0, 90.0);
    const double backwardMin =
        ReadEnvDouble("MBS_AUTOFULL_ZONE_BACKWARD_DEG", 150.0, 90.0, 180.0);
    if (thetaDeg <= forwardMax)
        return 0;
    if (thetaDeg >= backwardMin)
        return 2;
    return 1;
}

static const char *AutoFullThetaZoneName(int zone)
{
    if (zone < 0)
        return "mixed";
    if (zone == 0)
        return "forward";
    if (zone == 2)
        return "backward";
    return "middle";
}

static double AutoFullThetaZoneOrientFactor(int zone)
{
    if (zone == 0)
        return ReadEnvDouble("MBS_AUTOFULL_ZONE_FORWARD_ORIENT_FACTOR",
                             0.5, 0.015625, 1.0);
    if (zone == 2)
        return ReadEnvDouble("MBS_AUTOFULL_ZONE_BACKWARD_ORIENT_FACTOR",
                             1.0, 0.015625, 1.0);
    return ReadEnvDouble("MBS_AUTOFULL_ZONE_MIDDLE_ORIENT_FACTOR",
                         0.875, 0.015625, 1.0);
}

static double AutoFullRowOrientSafetyFactor(int zone)
{
    if (zone == 0)
        return ReadEnvDouble("MBS_AUTOFULL_ROW_ORIENT_FORWARD_SAFETY",
                             0.75, 0.0, 1.0);
    if (zone == 2)
        return ReadEnvDouble("MBS_AUTOFULL_ROW_ORIENT_BACKWARD_SAFETY",
                             0.0, 0.0, 1.0);
    return ReadEnvDouble("MBS_AUTOFULL_ROW_ORIENT_MIDDLE_SAFETY",
                         0.5, 0.0, 1.0);
}

static bool AutoFullBackscatterConeControlEnabled()
{
    const char *env = std::getenv("MBS_AUTOFULL_BACK_CONE_CONTROL");
    if (!env || !*env)
        return false;
    return !(env[0] == '0' && env[1] == '\0');
}

static double AutoFullBackscatterConeMinDeg()
{
    return ReadEnvDouble("MBS_AUTOFULL_BACK_CONE_DEG", 170.0, 120.0, 180.0);
}

static int AutoFullThetaZoneOrientLimit(int nOrient, int zone)
{
    if (!AutoFullThetaZonesEnabled())
        return nOrient;
    const double factor = AutoFullThetaZoneOrientFactor(zone);
    int raw = (int)std::floor(nOrient * factor + 0.5);
    raw = std::max(64, std::min(raw, nOrient));
    int limit = RoundDownToMultiple(raw, 16);
    return std::max(64, std::min(limit, nOrient));
}

struct AutoFullOldautoOrientEstimate
{
    int div = 1;
    int nBetaFull = 1;
    int nGammaFull = 1;
    int nBeta = 1;
    int nGamma = 1;
    int total = 1;
};

static int ClampOrientCount(long long value)
{
    if (value < 1)
        return 1;
    if (value > std::numeric_limits<int>::max())
        return std::numeric_limits<int>::max();
    return (int)value;
}

static int CeilDivInt(int value, int div)
{
    if (div < 1)
        div = 1;
    return (value + div - 1) / div;
}

static AutoFullOldautoOrientEstimate AutoFullOldautoEstimate(
    double betaSym, double gammaSym, double wave, double Dmax,
    int ringPoints, int div)
{
    AutoFullOldautoOrientEstimate est;
    est.div = std::max(1, div);
    const double deltaDeg = (Dmax > 0.0)
        ? 0.69 * wave / Dmax * (180.0 / M_PI)
        : 180.0;
    const double orientStep = deltaDeg / std::max(1, ringPoints);
    const double betaDeg = RadToDeg(betaSym);
    const double gammaDeg = RadToDeg(gammaSym);
    est.nBetaFull = std::max(1, (int)std::ceil(betaDeg / orientStep));
    est.nGammaFull = std::max(1, (int)std::ceil(gammaDeg / orientStep));
    est.nBeta = std::max(3, CeilDivInt(est.nBetaFull, est.div));
    est.nGamma = std::max(3, CeilDivInt(est.nGammaFull, est.div));
    est.total = ClampOrientCount((long long)est.nBeta * est.nGamma);
    return est;
}

static AutoFullOldautoOrientEstimate AutoFullOldautoEstimateForTarget(
    double betaSym, double gammaSym, double wave, double Dmax,
    int ringPoints, int targetOrient, bool evenGamma)
{
    AutoFullOldautoOrientEstimate best;
    bool haveBest = false;
    int bestActualTotal = 0;

    for (int div = 1; div <= 4096; ++div)
    {
        AutoFullOldautoOrientEstimate est = AutoFullOldautoEstimate(
            betaSym, gammaSym, wave, Dmax, ringPoints, div);
        if (evenGamma && (est.nGamma % 2 != 0))
            ++est.nGamma;
        est.total = ClampOrientCount((long long)est.nBeta * est.nGamma);
        const int actualTotal = ClampOrientCount(
            (long long)(est.nBeta + 1) * est.nGamma);

        if (actualTotal < targetOrient)
            continue;
        if (!haveBest || actualTotal < bestActualTotal)
        {
            best = est;
            bestActualTotal = actualTotal;
            haveBest = true;
        }
    }

    if (!haveBest)
    {
        best = AutoFullOldautoEstimate(betaSym, gammaSym, wave, Dmax,
                                       ringPoints, 1);
        if (evenGamma && (best.nGamma % 2 != 0))
            ++best.nGamma;
        best.total = ClampOrientCount((long long)best.nBeta * best.nGamma);
    }
    return best;
}

static int OldAutoFullFinalTargetOrient(int adaptiveTarget)
{
    const double factor = ReadEnvDouble("MBS_OLDAUTOFULL_ORIENT_FACTOR",
                                        8.0, 1.0, 64.0);
    const long long target = (long long)std::ceil(
        std::max(1, adaptiveTarget) * factor);
    return ClampOrientCount(target);
}

static int OldAutoFullMinFinalNphi()
{
    return ReadEnvInt("MBS_OLDAUTOFULL_MIN_NPHI", 600, 36, 2400);
}

static double AutoFullThetaZoneCutoffMultiplier(int zone)
{
    if (zone == 0)
        return ReadEnvDouble("MBS_AUTOFULL_ZONE_FORWARD_CUTOFF_MULT",
                             1.0, 1.0, 1000.0);
    if (zone == 2)
        return ReadEnvDouble("MBS_AUTOFULL_ZONE_BACKWARD_CUTOFF_MULT",
                             1.0, 1.0, 1000.0);
    return ReadEnvDouble("MBS_AUTOFULL_ZONE_MIDDLE_CUTOFF_MULT",
                         1.0, 1.0, 1000.0);
}

static double PreparedBeamImportance(const PreparedBeam &beam)
{
    if (beam.isExternal)
        return std::numeric_limits<double>::infinity();
    const double jn = beam.origBeam.J.Norm();
    const double importance = jn * beam.beam_area;
    return std::isfinite(importance) && importance > 0.0 ? importance : 0.0;
}

static std::vector<PreparedOrientation> BuildPreparedSliceForThetaZone(
    const std::vector<PreparedOrientation> &prepared,
    int start,
    int count,
    double weightScale,
    double extraImportanceRel)
{
    std::vector<PreparedOrientation> out;
    out.reserve(std::max(0, count));
    for (int i = 0; i < count; ++i)
    {
        PreparedOrientation po = prepared[start + i];
        po.sinZenith *= weightScale;
        if (extraImportanceRel > 0.0)
        {
            double maxImportance = 0.0;
            for (const PreparedBeam &beam : po.beams)
                if (!beam.isExternal)
                    maxImportance = std::max(maxImportance,
                                             PreparedBeamImportance(beam));
            if (maxImportance > 0.0)
            {
                std::vector<PreparedBeam> filtered;
                filtered.reserve(po.beams.size());
                const double threshold = maxImportance * extraImportanceRel;
                for (const PreparedBeam &beam : po.beams)
                {
                    if (beam.isExternal
                        || PreparedBeamImportance(beam) >= threshold)
                    {
                        filtered.push_back(beam);
                    }
                }
                po.beams.swap(filtered);
            }
        }
        out.push_back(po);
    }
    return out;
}

static long long CountPreparedBeams(
    const std::vector<PreparedOrientation> &prepared)
{
    long long count = 0;
    for (const PreparedOrientation &orientation : prepared)
        count += (long long)orientation.beams.size();
    return count;
}

struct AutoFullThetaWorkKey
{
    int nPhi;
    int nOrient;
    int zone;
    double beamCutoff;
    bool directFullPhi;

    bool operator<(const AutoFullThetaWorkKey &other) const
    {
        if (nOrient != other.nOrient)
            return nOrient < other.nOrient;
        if (zone != other.zone)
            return zone < other.zone;
        if (beamCutoff != other.beamCutoff)
            return beamCutoff < other.beamCutoff;
        if (directFullPhi != other.directFullPhi)
            return directFullPhi < other.directFullPhi;
        return nPhi < other.nPhi;
    }
};

struct MuellerRowAverage
{
    double v[16];
};

static std::vector<MuellerRowAverage> AverageMuellerRows(const Arr2D &arr,
                                                         int nAz,
                                                         int nZen)
{
    std::vector<MuellerRowAverage> rows(nZen + 1);
    for (int t = 0; t <= nZen; ++t)
    {
        for (int k = 0; k < 16; ++k)
            rows[t].v[k] = 0.0;
        for (int p = 0; p < nAz; ++p)
        {
            matrix m = arr(p, t);
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    rows[t].v[i * 4 + j] += m[i][j];
        }
        double inv = nAz > 0 ? 1.0 / nAz : 1.0;
        for (int k = 0; k < 16; ++k)
            rows[t].v[k] *= inv;
    }
    return rows;
}

static double MuellerRowRelativeError(const MuellerRowAverage &current,
                                      const MuellerRowAverage &previous)
{
    double maxErr = 0.0;
    double c11 = current.v[0];
    double p11 = previous.v[0];
    double m11Den = std::max(std::max(std::fabs(c11), std::fabs(p11)), 1e-30);
    maxErr = std::max(maxErr, std::fabs(c11 - p11) / m11Den);

    if (std::fabs(c11) < 1e-30 || std::fabs(p11) < 1e-30)
        return maxErr;

    for (int k = 1; k < 16; ++k)
    {
        double cn = current.v[k] / c11;
        double pn = previous.v[k] / p11;
        double den = std::max(1.0, std::max(std::fabs(cn), std::fabs(pn)));
        double err = std::fabs(cn - pn) / den;
        if (std::isfinite(err))
            maxErr = std::max(maxErr, err);
    }
    return maxErr;
}

static double MuellerRowRelativeErrorValues(const double *current,
                                            const double *previous,
                                            int *worstElement)
{
    double maxErr = 0.0;
    int worst = 0;
    double c11 = current[0];
    double p11 = previous[0];
    double m11Den = std::max(std::max(std::fabs(c11), std::fabs(p11)), 1e-30);
    maxErr = std::fabs(c11 - p11) / m11Den;

    if (std::fabs(c11) >= 1e-30 && std::fabs(p11) >= 1e-30)
    {
        for (int k = 1; k < 16; ++k)
        {
            double cn = current[k] / c11;
            double pn = previous[k] / p11;
            double den = std::max(1.0, std::max(std::fabs(cn), std::fabs(pn)));
            double err = std::fabs(cn - pn) / den;
            if (std::isfinite(err) && err > maxErr)
            {
                maxErr = err;
                worst = k;
            }
        }
    }

    if (worstElement)
        *worstElement = worst;
    return maxErr;
}

static double OwenMeanScalarRelativeError(const std::vector<double> &samples)
{
    const int n = (int)samples.size();
    if (n < 2)
        return 0.0;

    double sum = 0.0;
    for (double value : samples)
        sum += value;
    const double mean = sum / n;

    double sq = 0.0;
    for (double value : samples)
    {
        const double loo = (sum - value) / (n - 1);
        const double den = std::max(std::max(std::fabs(mean), std::fabs(loo)),
                                    1e-30);
        const double err = std::fabs(mean - loo) / den;
        if (std::isfinite(err))
            sq += err * err;
    }

    return std::sqrt(sq / n) / std::sqrt((double)n);
}

static double OwenMeanMuellerControlError(
    const std::vector<std::vector<double>> &seedCtrl,
    int ctrlRows,
    std::vector<double> *rowErrors,
    std::vector<int> *rowElements,
    int *worstRow,
    int *worstElement)
{
    const int seedCount = (int)seedCtrl.size();
    const int ctrlValues = ctrlRows * 16;
    if (seedCount < 2 || ctrlRows <= 0)
        return 0.0;

    std::vector<double> sum(ctrlValues, 0.0);
    for (const std::vector<double> &ctrl : seedCtrl)
    {
        if ((int)ctrl.size() < ctrlValues)
        {
            if (rowErrors)
                rowErrors->assign(ctrlRows,
                                  std::numeric_limits<double>::infinity());
            if (rowElements)
                rowElements->assign(ctrlRows, 0);
            if (worstRow)
                *worstRow = 0;
            if (worstElement)
                *worstElement = 0;
            return std::numeric_limits<double>::infinity();
        }
        for (int i = 0; i < ctrlValues; ++i)
            sum[i] += ctrl[i];
    }

    std::vector<double> mean(ctrlValues, 0.0);
    for (int i = 0; i < ctrlValues; ++i)
        mean[i] = sum[i] / seedCount;

    double maxEstimate = 0.0;
    int maxRow = 0;
    int maxElement = 0;
    if (rowErrors)
        rowErrors->assign(ctrlRows, 0.0);
    if (rowElements)
        rowElements->assign(ctrlRows, 0);

    for (int row = 0; row < ctrlRows; ++row)
    {
        double rowSq = 0.0;
        double rowWorstSeedErr = 0.0;
        int rowWorstElement = 0;
        const int base = row * 16;

        for (int s = 0; s < seedCount; ++s)
        {
            double loo[16];
            for (int e = 0; e < 16; ++e)
                loo[e] = (sum[base + e] - seedCtrl[s][base + e])
                    / (seedCount - 1);

            int element = 0;
            const double err = MuellerRowRelativeErrorValues(
                &mean[base], loo, &element);
            if (std::isfinite(err))
            {
                rowSq += err * err;
                if (err > rowWorstSeedErr)
                {
                    rowWorstSeedErr = err;
                    rowWorstElement = element;
                }
            }
        }

        const double estimate =
            std::sqrt(rowSq / seedCount) / std::sqrt((double)seedCount);
        if (rowErrors)
            (*rowErrors)[row] = estimate;
        if (rowElements)
            (*rowElements)[row] = rowWorstElement;
        if (std::isfinite(estimate) && estimate > maxEstimate)
        {
            maxEstimate = estimate;
            maxRow = row;
            maxElement = rowWorstElement;
        }
    }

    if (worstRow)
        *worstRow = maxRow;
    if (worstElement)
        *worstElement = maxElement;
    return maxEstimate;
}

static double OwenMeanMuellerSelectedRowsError(
    const std::vector<std::vector<double>> &seedCtrl,
    int ctrlRows,
    const std::vector<int> &activeRows,
    std::vector<double> *rowErrors,
    std::vector<int> *rowElements,
    int *worstRow,
    int *worstElement)
{
    const int seedCount = (int)seedCtrl.size();
    const int ctrlValues = ctrlRows * 16;
    if (seedCount < 2 || ctrlRows <= 0 || activeRows.empty())
        return 0.0;

    std::vector<double> sum(ctrlValues, 0.0);
    for (const std::vector<double> &ctrl : seedCtrl)
    {
        if ((int)ctrl.size() < ctrlValues)
        {
            if (rowErrors)
                rowErrors->assign(ctrlRows,
                                  std::numeric_limits<double>::infinity());
            if (rowElements)
                rowElements->assign(ctrlRows, 0);
            if (worstRow)
                *worstRow = activeRows.front();
            if (worstElement)
                *worstElement = 0;
            return std::numeric_limits<double>::infinity();
        }
        for (int i = 0; i < ctrlValues; ++i)
            sum[i] += ctrl[i];
    }

    std::vector<double> mean(ctrlValues, 0.0);
    for (int i = 0; i < ctrlValues; ++i)
        mean[i] = sum[i] / seedCount;

    if (rowErrors)
        rowErrors->assign(ctrlRows, 0.0);
    if (rowElements)
        rowElements->assign(ctrlRows, 0);

    double maxEstimate = 0.0;
    int maxRow = activeRows.front();
    int maxElement = 0;
    for (int row : activeRows)
    {
        if (row < 0 || row >= ctrlRows)
            continue;
        double rowSq = 0.0;
        double rowWorstSeedErr = 0.0;
        int rowWorstElement = 0;
        const int base = row * 16;
        for (int s = 0; s < seedCount; ++s)
        {
            double loo[16];
            for (int e = 0; e < 16; ++e)
                loo[e] = (sum[base + e] - seedCtrl[s][base + e])
                    / (seedCount - 1);

            int element = 0;
            const double err = MuellerRowRelativeErrorValues(
                &mean[base], loo, &element);
            if (std::isfinite(err))
            {
                rowSq += err * err;
                if (err > rowWorstSeedErr)
                {
                    rowWorstSeedErr = err;
                    rowWorstElement = element;
                }
            }
        }

        const double estimate =
            std::sqrt(rowSq / seedCount) / std::sqrt((double)seedCount);
        if (rowErrors)
            (*rowErrors)[row] = estimate;
        if (rowElements)
            (*rowElements)[row] = rowWorstElement;
        if (std::isfinite(estimate) && estimate > maxEstimate)
        {
            maxEstimate = estimate;
            maxRow = row;
            maxElement = rowWorstElement;
        }
    }

    if (worstRow)
        *worstRow = maxRow;
    if (worstElement)
        *worstElement = maxElement;
    return maxEstimate;
}

static double OwenMeanMuellerControlError(
    const std::vector<std::vector<double>> &seedCtrl,
    int ctrlRows,
    int *worstRow,
    int *worstElement)
{
    return OwenMeanMuellerControlError(seedCtrl, ctrlRows, nullptr, nullptr,
                                      worstRow, worstElement);
}

static std::string MuellerElementName(int element)
{
    element = std::max(0, std::min(15, element));
    std::ostringstream ss;
    ss << 'M' << (element / 4 + 1) << (element % 4 + 1);
    return ss.str();
}

static int AutoFullPilotOrientations(double eps)
{
    int fallback = 128;
    if (eps > 0.0 && eps <= 0.0025)
        fallback = 512;
    else if (eps > 0.0 && eps <= 0.01)
        fallback = 256;
    return ReadEnvInt("MBS_AUTOFULL_PILOT_ORIENT", fallback, 32, 4096);
}

static int AutoFullMaxNphi(int initialNphi)
{
    int fallback = std::max(600, RoundUpToMultiple(initialNphi, 6));
    return ReadEnvInt("MBS_AUTOFULL_MAX_NPHI", fallback, 36, 2400);
}

static std::vector<int> AutoFullPhiCandidates(int maxNphi)
{
    std::vector<int> candidates;
    const char *env = std::getenv("MBS_AUTOFULL_PHI_CANDIDATES");
    if (env && *env)
    {
        std::string text(env);
        for (char &ch : text)
            if (ch == ',' || ch == ';' || ch == ':')
                ch = ' ';
        std::stringstream ss(text);
        int value = 0;
        while (ss >> value)
        {
            if (value >= 2)
                candidates.push_back(value);
        }
    }
    if (candidates.empty())
        candidates = {36, 48, 60, 72, 90, 108, 132, 156, 180,
                      216, 252, 300, 360, 420, 480, 540, 600};

    candidates.push_back(maxNphi);
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()),
                     candidates.end());
    while (!candidates.empty() && candidates.back() > maxNphi)
        candidates.pop_back();
    if (candidates.empty() || candidates.front() != 36)
        candidates.insert(candidates.begin(), 36);
    return candidates;
}

static double AutoFullBoundedQualityEps(double eps)
{
    double value = eps > 0.0 ? eps : 0.01;
    return std::max(1e-4, std::min(value, 0.01));
}

static double AutoFullPhiSafetyFactor()
{
    return ReadEnvDouble("MBS_AUTOFULL_PHI_SAFETY", 0.85, 0.1, 1.0);
}

static bool AutoFullFullPhiValidationEnabled(double eps)
{
    const bool defaultEnabled = eps > 0.0 && eps <= 0.01;
    const char *env = std::getenv("MBS_AUTOFULL_FULL_PHI_VALIDATE");
    if (!env || !*env)
        return defaultEnabled;
    return !(env[0] == '0' && env[1] == '\0');
}

static bool AutoFullStrictFullPhiFinalEnabled(double eps)
{
    const bool defaultEnabled = eps > 0.0 && eps <= 0.01;
    const char *env = std::getenv("MBS_AUTOFULL_STRICT_FULL_PHI");
    if (!env || !*env)
        return defaultEnabled;
    return !(env[0] == '0' && env[1] == '\0');
}

static int AutoFullStrictMinFinalOrient(double eps)
{
    (void)eps;
    return ReadEnvInt("MBS_AUTOFULL_MIN_FINAL_ORIENT", 0, 0, 1048576);
}

static int AutoFullOldautoActualOrient(const AutoFullOldautoOrientEstimate &est)
{
    return ClampOrientCount((long long)(est.nBeta + 1) * est.nGamma);
}

static bool AutoFullBeamCutoffAutoEnabled()
{
    const char *env = std::getenv("MBS_AUTOFULL_BEAM_CUTOFF_AUTO");
    if (!env || !*env)
        return false;
    return !(env[0] == '0' && env[1] == '\0');
}

static double AutoFullBeamCutoffSafetyFactor(int zone)
{
    if (zone == 0)
        return ReadEnvDouble("MBS_AUTOFULL_BEAM_CUTOFF_FORWARD_SAFETY",
                             0.50, 0.0, 1.0);
    if (zone == 2)
        return ReadEnvDouble("MBS_AUTOFULL_BEAM_CUTOFF_BACKWARD_SAFETY",
                             0.20, 0.0, 1.0);
    return ReadEnvDouble("MBS_AUTOFULL_BEAM_CUTOFF_MIDDLE_SAFETY",
                         0.35, 0.0, 1.0);
}

static std::vector<double> AutoFullBeamCutoffCandidates(double baseCutoff)
{
    std::vector<double> candidates;
    const char *env = std::getenv("MBS_AUTOFULL_BEAM_CUTOFF_CANDIDATES");
    if (env && *env)
    {
        std::string text(env);
        for (char &ch : text)
            if (ch == ',' || ch == ';' || ch == ':')
                ch = ' ';
        std::stringstream ss(text);
        double value = 0.0;
        while (ss >> value)
        {
            if (value > 0.0 && std::isfinite(value))
                candidates.push_back(value);
        }
    }
    else if (baseCutoff > 0.0)
    {
        const double multipliers[] = {2.0, 4.0, 8.0, 16.0, 32.0};
        for (double mult : multipliers)
            candidates.push_back(baseCutoff * mult);
    }

    candidates.erase(
        std::remove_if(candidates.begin(), candidates.end(),
                       [baseCutoff](double value)
                       {
                           return !(value > baseCutoff * (1.0 + 1e-12))
                               || value > 0.25;
                       }),
        candidates.end());
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()),
                     candidates.end());
    return candidates;
}

static bool EnvFlagEnabled(const char *name, bool defaultValue)
{
    const char *value = std::getenv(name);
    if (!value || !*value)
        return defaultValue;
    return !(value[0] == '0' && value[1] == '\0');
}

static int FftPhiFactorOverrideValue()
{
    const char *value = std::getenv("MBS_FFT_PHI_FACTOR");
    if (!value || !*value)
        return 0;
    if ((value[0] == 'a' || value[0] == 'A')
        && (value[1] == 'u' || value[1] == 'U')
        && (value[2] == 't' || value[2] == 'T')
        && (value[3] == 'o' || value[3] == 'O')
        && value[4] == '\0')
        return 0;
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    if (!end || *end != '\0' || parsed < 1 || parsed > 64)
        return 0;
    return (int)parsed;
}

static int AutoFullFinalAveragePhi(int nPhi, double eps)
{
    int factor = FftPhiFactorOverrideValue();
    if (factor <= 0)
        factor = HandlerPO::SelectAutoFftPhiFactor(nPhi, eps);
    if (factor <= 1 || nPhi < 32)
        return nPhi;

    int nLow = nPhi / factor;
    if (nLow < 16)
        nLow = std::min(nPhi, 16);
    if (nLow >= nPhi)
        return nPhi;
    return nLow;
}

static bool AutoFullLowPhiAverageEnabled(const HandlerPO *handler, double eps)
{
    return handler
        && handler->IsGpuEnabled()
        && handler->IsFftEnabled()
        && eps > 0.0
        && EnvFlagEnabled("MBS_AUTOFULL_FINAL_LOWPHI_AVG", false);
}

static bool AutoFullDirectFullPhiEnabled(const HandlerPO *handler)
{
    return handler
        && handler->IsGpuEnabled()
        && handler->IsFftEnabled()
        && EnvFlagEnabled("MBS_AUTOFULL_DIRECT_FULL_NPHI", true);
}

static int AutoFullVariableCalcPhi(int nPhi, double eps, int minDirectPhi,
                                   bool canUseFftAverage,
                                   bool directFullPhi)
{
    if (directFullPhi)
        return nPhi;
    if (!canUseFftAverage)
        return nPhi;
    return std::max(AutoFullFinalAveragePhi(nPhi, eps), minDirectPhi);
}

static int AutoFullLowPhiDirectFloor(int outputNphi, double eps)
{
    return std::max(2, AutoFullFinalAveragePhi(outputNphi, eps));
}

static std::vector<int> BuildThetaControlIndices(const ScatteringRange &sphere,
                                                 double eps)
{
    std::vector<int> indices;
    auto addIndex = [&](int idx)
    {
        idx = std::max(0, std::min(idx, sphere.nZenith));
        indices.push_back(idx);
    };

    const int nTheta = sphere.nZenith + 1;
    int fallbackAllRowsMax = 64;
    if (eps > 0.0 && eps <= 0.01)
        fallbackAllRowsMax = 1024;
    if (eps > 0.0 && eps <= 0.005)
        fallbackAllRowsMax = std::max(fallbackAllRowsMax, 256);
    const int allRowsMax = ReadEnvInt("MBS_AUTOFULL_ORIENT_CONTROL_ALL_MAX",
                                      fallbackAllRowsMax, 0, 1000000);
    if (EnvFlagEnabled("MBS_AUTOFULL_ORIENT_CONTROL_ALL", false)
        || nTheta <= allRowsMax)
    {
        for (int j = 0; j <= sphere.nZenith; ++j)
            addIndex(j);
        return indices;
    }

    // Controls are selected from the actual grid, not from hard-coded halo
    // angles belonging to one particle family.
    addIndex(0);
    addIndex(sphere.nZenith);
    if (sphere.nZenith >= 1)
    {
        addIndex(1);
        addIndex(sphere.nZenith - 1);
    }

    int fallbackEven = 9;
    if (sphere.nZenith + 1 >= 512)
        fallbackEven = 65;
    else if (sphere.nZenith + 1 >= 128)
        fallbackEven = 33;
    else if (sphere.nZenith + 1 >= 64)
        fallbackEven = 17;
    if (eps > 0.0 && eps <= 0.005)
        fallbackEven = std::max(fallbackEven, 129);
    int nEven = std::min(ReadEnvInt("MBS_AUTOFULL_ORIENT_CONTROL_ROWS",
                                    fallbackEven, 9, 1025),
                         sphere.nZenith + 1);
    for (int k = 0; k < nEven; ++k)
    {
        int idx = nEven > 1
            ? (int)std::lround((double)k * sphere.nZenith / (nEven - 1))
            : 0;
        addIndex(idx);
    }

    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    return indices;
}

static std::vector<int> BuildBackscatterConeControlRows(
    const ScatteringRange &sphere,
    const std::vector<int> &ctrlIdx)
{
    std::vector<int> rows;
    if (!AutoFullBackscatterConeControlEnabled())
        return rows;
    const double minDeg = AutoFullBackscatterConeMinDeg();
    for (int row = 0; row < (int)ctrlIdx.size(); ++row)
    {
        int idx = std::max(0, std::min(ctrlIdx[row], sphere.nZenith));
        if (RadToDeg(sphere.GetZenith(idx)) >= minDeg)
            rows.push_back(row);
    }
    if ((int)rows.size() < 2)
        rows.clear();
    return rows;
}

static std::vector<double> AggregateControlRows(
    const std::vector<double> &values,
    const std::vector<int> &rows)
{
    std::vector<double> out(16, 0.0);
    if (rows.empty())
        return out;
    for (int row : rows)
    {
        const int base = row * 16;
        if (base < 0 || base + 15 >= (int)values.size())
            continue;
        for (int e = 0; e < 16; ++e)
            out[e] += values[base + e];
    }
    const double inv = 1.0 / (double)rows.size();
    for (double &value : out)
        value *= inv;
    return out;
}

static std::vector<std::vector<double>> AggregateSeedControlRows(
    const std::vector<std::vector<double>> &seedCtrl,
    const std::vector<int> &rows)
{
    std::vector<std::vector<double>> out;
    out.reserve(seedCtrl.size());
    for (const std::vector<double> &ctrl : seedCtrl)
        out.push_back(AggregateControlRows(ctrl, rows));
    return out;
}

static void ApplyForwardPoleSymmetryLocal(matrix &m)
{
    const double M00 = m[0][0];
    const double M11 = 0.5 * (m[1][1] + m[2][2]);
    const double M33 = m[3][3];
    const double M03 = m[0][3];

    m.Fill(0.0);
    m[0][0] = M00;
    m[1][1] = M11;
    m[2][2] = M11;
    m[3][3] = M33;
    m[0][3] = M03;
    m[3][0] = M03;
}

static void ApplyBackwardPoleSymmetryLocal(matrix &m)
{
    const double M00 = m[0][0];
    const double M11 = 0.5 * (m[1][1] - m[2][2]);
    const double M33 = m[3][3];
    const double M03 = m[0][3];

    m.Fill(0.0);
    m[0][0] = M00;
    m[1][1] = M11;
    m[2][2] = -M11;
    m[3][3] = M33;
    m[0][3] = M03;
    m[3][0] = M03;
}

static matrix AzimuthAverageOutputMatrix(const Arr2D &src,
                                         const ScatteringRange &sphere,
                                         int nAz,
                                         int localRow)
{
    matrix sum(4, 4);
    sum.Fill(0.0);
    matrix Lp(4, 4);
    Lp.Fill(0.0);
    for (int i = 0; i < 4; ++i)
        Lp[i][i] = 1.0;

    const double theta = sphere.GetZenith(localRow);
    const bool isForwardPole = theta < __FLT_EPSILON__;
    const bool isBackwardPole = theta > M_PI - __FLT_EPSILON__;

    for (int p = 0; p < nAz; ++p)
    {
        const double radAz = -p * sphere.azinuthStep;
        matrix m = src(p, localRow);
        Lp[1][1] = cos(2 * radAz);
        Lp[1][2] = sin(2 * radAz);
        Lp[2][1] = -Lp[1][2];
        Lp[2][2] = Lp[1][1];

        if (isForwardPole || isBackwardPole)
            sum += m;
        else
            sum += m * Lp;
    }
    sum /= nAz;
    if (isBackwardPole)
        ApplyBackwardPoleSymmetryLocal(sum);
    else if (isForwardPole)
        ApplyForwardPoleSymmetryLocal(sum);
    return sum;
}

struct PreparedProbeResult
{
    std::vector<MuellerRowAverage> rows;
    double scatteringIntegral = 0.0;
    bool usedGpu = false;
};

static PreparedProbeResult EvaluatePreparedProbe(
    HandlerPO *handler,
    const std::vector<PreparedOrientation> &prepared,
    ScatteringRange sphere,
    const Light &incidentLight,
    bool mirrorGamma)
{
    PreparedProbeResult result;
    if (!handler || prepared.empty() || sphere.nAzimuth < 1
        || sphere.nZenith < 0)
        return result;

    sphere.ComputeSphereDirections(incidentLight);
    const ScatteringRange savedSphere = handler->m_sphere;
    handler->SetScatteringSphere(sphere);

    const int nAz = sphere.nAzimuth;
    const int nZen = sphere.nZenith;
    Arr2D localM(nAz + 1, nZen + 1, 4, 4);
    localM.ClearArr();
    Arr2D localMNoShadow(nAz + 1, nZen + 1, 4, 4);
    localMNoShadow.ClearArr();

    bool usedGpu = false;
    if (handler->IsGpuEnabled())
    {
        usedGpu = true;
        for (int start = 0; start < (int)prepared.size(); )
        {
            const int batch = handler->SelectGpuOrientationBatchSize(
                prepared, start, (int)prepared.size() - start);
            const int end = std::min(start + std::max(1, batch),
                                     (int)prepared.size());
            if (!handler->HandleOrientationsToLocalGpu(
                    prepared, start, end - start, localM, localMNoShadow))
            {
                usedGpu = false;
                localM.ClearArr();
                localMNoShadow.ClearArr();
                break;
            }
            start = end;
        }
    }

    if (!usedGpu)
    {
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            Arr2D threadM(nAz + 1, nZen + 1, 4, 4);
            threadM.ClearArr();
            Arr2D threadMNoShadow(nAz + 1, nZen + 1, 4, 4);
            threadMNoShadow.ClearArr();
            std::vector<Arr2DC> localJ;
            std::vector<Arr2DC> localJNoShadow;
            if (handler->isCoh)
            {
                Arr2DC tmp(nAz + 1, nZen + 1, 2, 2);
                tmp.ClearArr();
                localJ.push_back(tmp);
                Arr2DC tmpNoShadow(nAz + 1, nZen + 1, 2, 2);
                tmpNoShadow.ClearArr();
                localJNoShadow.push_back(tmpNoShadow);
            }

#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 1)
#endif
            for (int i = 0; i < (int)prepared.size(); ++i)
            {
                if (!prepared[i].beams.empty())
                    handler->HandleBeamsToLocal(
                        prepared[i], threadM, localJ,
                        handler->isCoh ? &localJNoShadow : nullptr);
                if (handler->isCoh && !localJ.empty())
                {
                    const double weight = prepared[i].sinZenith;
                    HandlerPO::AddToMuellerLocal(localJ, weight, threadM,
                                                 nAz, nZen);
                    HandlerPO::AddToMuellerLocal(localJNoShadow, weight,
                                                 threadMNoShadow, nAz, nZen);
                    localJ[0].ClearArr();
                    localJNoShadow[0].ClearArr();
                }
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                for (int p = 0; p < nAz; ++p)
                    for (int t = 0; t <= nZen; ++t)
                    {
                        localM.insert(p, t, threadM(p, t));
                        localMNoShadow.insert(p, t,
                                              threadMNoShadow(p, t));
                    }
            }
        }
    }

    if (mirrorGamma)
        ApplyMirrorGammaMueller(localM, nAz, nZen);

    result.rows.resize(nZen + 1);
    for (int t = 0; t <= nZen; ++t)
    {
        const matrix averaged = AzimuthAverageOutputMatrix(
            localM, sphere, nAz, t);
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                result.rows[t].v[r * 4 + c] = averaged[r][c];
        result.scatteringIntegral += averaged[0][0]
            * sphere.Compute2PiDcos(t);
    }
    result.usedGpu = usedGpu;

    handler->SetScatteringSphere(savedSphere);
    return result;
}

static double ProbeRowsRelativeError(
    const std::vector<MuellerRowAverage> &current,
    const std::vector<MuellerRowAverage> &previous,
    const ScatteringRange &sphere,
    double currentIntegral,
    double previousIntegral,
    double *worstThetaDeg,
    int *worstElement)
{
    double maxError = 0.0;
    int maxRow = 0;
    int maxElement = 0;
    const int rows = std::min((int)current.size(), (int)previous.size());
    for (int row = 0; row < rows; ++row)
    {
        int element = 0;
        const double error = MuellerRowRelativeErrorValues(
            current[row].v, previous[row].v, &element);
        if (std::isfinite(error) && error > maxError)
        {
            maxError = error;
            maxRow = row;
            maxElement = element;
        }
    }

    const double integralDen = std::max(
        std::max(std::fabs(currentIntegral), std::fabs(previousIntegral)),
        1e-30);
    const double integralError = std::fabs(currentIntegral - previousIntegral)
        / integralDen;
    if (std::isfinite(integralError) && integralError > maxError)
    {
        maxError = integralError;
        maxElement = 16; // Integrated M11 / Qsca criterion.
    }

    if (worstThetaDeg)
        *worstThetaDeg = rows > 0
            ? RadToDeg(sphere.GetZenith(std::min(maxRow, rows - 1))) : 0.0;
    if (worstElement)
        *worstElement = maxElement;
    return maxError;
}

static double AdaptiveThetaIntegral(
    const std::map<double, MuellerRowAverage> &values)
{
    if (values.size() < 2)
        return 0.0;
    ScatteringRange sphere(values.begin()->first,
                           values.rbegin()->first, 1,
                           (int)values.size() - 1);
    sphere.isNonUniform = true;
    sphere.thetaValues.reserve(values.size());
    for (const auto &entry : values)
        sphere.thetaValues.push_back(entry.first);
    sphere.nZenith = (int)sphere.thetaValues.size() - 1;
    sphere.zenithStart = sphere.thetaValues.front();
    sphere.zenithEnd = sphere.thetaValues.back();

    double integral = 0.0;
    int row = 0;
    for (const auto &entry : values)
    {
        integral += entry.second.v[0] * sphere.Compute2PiDcos(row);
        ++row;
    }
    return integral;
}

static double MuellerMatrixRelativeError(const matrix &current,
                                         const matrix &previous,
                                         int *worstElement)
{
    double currentValues[16];
    double previousValues[16];
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
        {
            currentValues[r * 4 + c] = current[r][c];
            previousValues[r * 4 + c] = previous[r][c];
        }
    return MuellerRowRelativeErrorValues(currentValues, previousValues,
                                         worstElement);
}

static double WriteOwenSpreadByThetaFile(
    const std::string &fileName,
    const ScatteringRange &sphere,
    const std::vector<std::vector<matrix>> &seedRows)
{
    if (seedRows.size() < 2)
        return 0.0;

    std::ofstream out(fileName.c_str(), std::ios::out);
    out << std::setprecision(10);
    out << "# theta_deg m11_pair_spread_pct max_mueller_pair_spread_pct"
        << " pair_worst_element m11_mean_error_est_pct"
        << " max_mueller_mean_error_est_pct mean_worst_element\n";

    double globalMueller = 0.0;
    double globalM11 = 0.0;
    int globalTheta = 0;
    int globalElement = 0;
    double globalMeanMueller = 0.0;
    double globalMeanM11 = 0.0;
    int globalMeanTheta = 0;
    int globalMeanElement = 0;
    std::vector<int> backConeRows;
    std::vector<char> isBackConeRow(sphere.nZenith + 1, 0);
    if (AutoFullBackscatterConeControlEnabled())
    {
        const double minDeg = AutoFullBackscatterConeMinDeg();
        for (int t = 0; t <= sphere.nZenith; ++t)
        {
            if (RadToDeg(sphere.GetZenith(t)) >= minDeg)
            {
                backConeRows.push_back(t);
                isBackConeRow[t] = 1;
            }
        }
        if ((int)backConeRows.size() < 2)
        {
            backConeRows.clear();
            std::fill(isBackConeRow.begin(), isBackConeRow.end(), 0);
        }
    }

    for (int t = 0; t <= sphere.nZenith; ++t)
    {
        double rowM11 = 0.0;
        double rowMueller = 0.0;
        int rowElement = 0;
        double rowMeanM11 = 0.0;
        double rowMeanMueller = 0.0;
        int rowMeanElement = 0;

        for (size_t a = 0; a < seedRows.size(); ++a)
            for (size_t b = a + 1; b < seedRows.size(); ++b)
            {
                const matrix &ma = seedRows[a][t];
                const matrix &mb = seedRows[b][t];
                const double den = std::max(std::max(std::fabs(ma[0][0]),
                                                     std::fabs(mb[0][0])),
                                            1e-30);
                const double m11 = std::fabs(ma[0][0] - mb[0][0]) / den;
                rowM11 = std::max(rowM11, m11);

                int element = 0;
                const double mueller =
                    MuellerMatrixRelativeError(ma, mb, &element);
                if (std::isfinite(mueller) && mueller > rowMueller)
                {
                    rowMueller = mueller;
                    rowElement = element;
                }
            }

        if (seedRows.size() >= 3)
        {
            const int seedCount = (int)seedRows.size();
            double sum[16] = {0.0};
            double mean[16] = {0.0};
            for (size_t s = 0; s < seedRows.size(); ++s)
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c)
                        sum[r * 4 + c] += seedRows[s][t][r][c];
            for (int e = 0; e < 16; ++e)
                mean[e] = sum[e] / seedCount;

            double m11Sq = 0.0;
            double muellerSq = 0.0;
            double rowWorstSeedErr = 0.0;
            for (size_t s = 0; s < seedRows.size(); ++s)
            {
                double loo[16];
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c)
                    {
                        const int e = r * 4 + c;
                        loo[e] = (sum[e] - seedRows[s][t][r][c])
                            / (seedCount - 1);
                    }

                const double denM11 = std::max(std::max(std::fabs(mean[0]),
                                                        std::fabs(loo[0])),
                                               1e-30);
                const double m11Err = std::fabs(mean[0] - loo[0]) / denM11;
                if (std::isfinite(m11Err))
                    m11Sq += m11Err * m11Err;

                int element = 0;
                const double muellerErr =
                    MuellerRowRelativeErrorValues(mean, loo, &element);
                if (std::isfinite(muellerErr))
                {
                    muellerSq += muellerErr * muellerErr;
                    if (muellerErr > rowWorstSeedErr)
                    {
                        rowWorstSeedErr = muellerErr;
                        rowMeanElement = element;
                    }
                }
            }

            rowMeanM11 = std::sqrt(m11Sq / seedCount)
                / std::sqrt((double)seedCount);
            rowMeanMueller = std::sqrt(muellerSq / seedCount)
                / std::sqrt((double)seedCount);
        }

        const bool exactRowControlsGlobal =
            t >= (int)isBackConeRow.size() || !isBackConeRow[t];
        if (exactRowControlsGlobal && rowMueller > globalMueller)
        {
            globalMueller = rowMueller;
            globalM11 = rowM11;
            globalTheta = t;
            globalElement = rowElement;
        }
        if (exactRowControlsGlobal && rowMeanMueller > globalMeanMueller)
        {
            globalMeanMueller = rowMeanMueller;
            globalMeanM11 = rowMeanM11;
            globalMeanTheta = t;
            globalMeanElement = rowMeanElement;
        }

        out << RadToDeg(sphere.GetZenith(t)) << ' '
            << rowM11 * 100.0 << ' '
            << rowMueller * 100.0 << ' '
            << MuellerElementName(rowElement) << ' '
            << rowMeanM11 * 100.0 << ' '
            << rowMeanMueller * 100.0 << ' '
            << MuellerElementName(rowMeanElement) << '\n';
    }

    if (!backConeRows.empty())
    {
        std::vector<matrix> coneSeed;
        coneSeed.reserve(seedRows.size());
        for (size_t s = 0; s < seedRows.size(); ++s)
            coneSeed.push_back(matrix(4, 4));
        for (matrix &m : coneSeed)
            m.Fill(0.0);
        for (size_t s = 0; s < seedRows.size(); ++s)
        {
            for (int t : backConeRows)
                coneSeed[s] += seedRows[s][t];
            coneSeed[s] *= (1.0 / (double)backConeRows.size());
        }

        double coneM11 = 0.0;
        double coneMueller = 0.0;
        int coneElement = 0;
        for (size_t a = 0; a < coneSeed.size(); ++a)
            for (size_t b = a + 1; b < coneSeed.size(); ++b)
            {
                const matrix &ma = coneSeed[a];
                const matrix &mb = coneSeed[b];
                const double den = std::max(std::max(std::fabs(ma[0][0]),
                                                     std::fabs(mb[0][0])),
                                            1e-30);
                coneM11 = std::max(coneM11,
                    std::fabs(ma[0][0] - mb[0][0]) / den);
                int element = 0;
                const double err = MuellerMatrixRelativeError(ma, mb,
                                                              &element);
                if (std::isfinite(err) && err > coneMueller)
                {
                    coneMueller = err;
                    coneElement = element;
                }
            }

        double coneMeanM11 = 0.0;
        double coneMeanMueller = 0.0;
        int coneMeanElement = 0;
        if (coneSeed.size() >= 3)
        {
            const int seedCount = (int)coneSeed.size();
            double sum[16] = {0.0};
            double mean[16] = {0.0};
            for (const matrix &m : coneSeed)
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c)
                        sum[r * 4 + c] += m[r][c];
            for (int e = 0; e < 16; ++e)
                mean[e] = sum[e] / seedCount;

            double m11Sq = 0.0;
            double muellerSq = 0.0;
            double worstSeedErr = 0.0;
            for (const matrix &m : coneSeed)
            {
                double loo[16];
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c)
                    {
                        const int e = r * 4 + c;
                        loo[e] = (sum[e] - m[r][c]) / (seedCount - 1);
                    }
                const double denM11 = std::max(std::max(std::fabs(mean[0]),
                                                        std::fabs(loo[0])),
                                               1e-30);
                const double m11Err = std::fabs(mean[0] - loo[0]) / denM11;
                if (std::isfinite(m11Err))
                    m11Sq += m11Err * m11Err;

                int element = 0;
                const double muellerErr =
                    MuellerRowRelativeErrorValues(mean, loo, &element);
                if (std::isfinite(muellerErr))
                {
                    muellerSq += muellerErr * muellerErr;
                    if (muellerErr > worstSeedErr)
                    {
                        worstSeedErr = muellerErr;
                        coneMeanElement = element;
                    }
                }
            }
            coneMeanM11 = std::sqrt(m11Sq / seedCount)
                / std::sqrt((double)seedCount);
            coneMeanMueller = std::sqrt(muellerSq / seedCount)
                / std::sqrt((double)seedCount);
        }

        if (coneMueller > globalMueller)
        {
            globalMueller = coneMueller;
            globalM11 = coneM11;
            globalTheta = backConeRows.back();
            globalElement = coneElement;
        }
        if (coneMeanMueller > globalMeanMueller)
        {
            globalMeanMueller = coneMeanMueller;
            globalMeanM11 = coneMeanM11;
            globalMeanTheta = backConeRows.back();
            globalMeanElement = coneMeanElement;
        }
        out << "# backscatter_cone theta_min_deg "
            << AutoFullBackscatterConeMinDeg()
            << " rows " << backConeRows.size()
            << " m11_pair_spread_pct " << coneM11 * 100.0
            << " max_mueller_pair_spread_pct " << coneMueller * 100.0
            << " pair_worst_element " << MuellerElementName(coneElement)
            << " m11_mean_error_est_pct " << coneMeanM11 * 100.0
            << " max_mueller_mean_error_est_pct "
            << coneMeanMueller * 100.0
            << " mean_worst_element "
            << MuellerElementName(coneMeanElement) << '\n';
    }
    out.close();

    std::cout << "Autofull final Owen spread by theta written: "
              << fileName << std::endl;
    std::cout << "Worst final Owen inter-seed theta spread: "
              << std::fixed << std::setprecision(2)
              << globalMueller * 100.0 << "% @ "
              << RadToDeg(sphere.GetZenith(globalTheta)) << "deg/"
              << MuellerElementName(globalElement)
              << " (M11 spread " << globalM11 * 100.0 << "%)"
              << std::endl;
    if (seedRows.size() >= 3)
    {
        std::cout << "Worst final Owen mean-error estimate: "
                  << std::fixed << std::setprecision(2)
                  << globalMeanMueller * 100.0 << "% @ "
                  << RadToDeg(sphere.GetZenith(globalMeanTheta)) << "deg/"
                  << MuellerElementName(globalMeanElement)
                  << " (M11 estimate " << globalMeanM11 * 100.0 << "%)"
                  << std::endl;
    }
    return seedRows.size() >= 3 ? globalMeanMueller : 0.0;
}

static void WriteAveragedRowsFile(const std::string &destName,
                                  const ScatteringRange &sphere,
                                  const std::vector<matrix> &rows,
                                  double incomingEnergy,
                                  double cExtOt,
                                  bool hasExtinctionOt,
                                  bool hasAbsorption,
                                  std::string &integralSummary)
{
    std::ofstream outFile(destName + ".dat", std::ios::out);
    if (!outFile.is_open())
        throw std::runtime_error(
            "cannot open Mueller output '" + destName
            + ".dat'.\n  Fix: verify output permissions and free disk space.");
    outFile << std::setprecision(10);
    outFile << "ScAngle 2pi*dcos "
            << "M11 M12 M13 M14 "
            << "M21 M22 M23 M24 "
            << "M31 M32 M33 M34 "
            << "M41 M42 M43 M44";

    double cscaIntegral = 0.0;
    std::vector<double> thetaRows;
    std::vector<double> m11Rows;
    thetaRows.reserve(sphere.nZenith + 1);
    m11Rows.reserve(sphere.nZenith + 1);
    for (int t = 0; t <= sphere.nZenith; ++t)
    {
        double dcos = sphere.Compute2PiDcos(t);
        cscaIntegral += rows[t][0][0] * dcos;
        thetaRows.push_back(sphere.GetZenith(t));
        m11Rows.push_back(rows[t][0][0]);
        outFile << '\n' << RadToDeg(sphere.GetZenith(t)) << ' ' << dcos << ' ';
        outFile << rows[t];
    }
    outFile.close();
    if (!outFile)
        throw std::runtime_error(
            "failed while writing Mueller output '" + destName
            + ".dat'.\n  Fix: verify free disk space and filesystem health.");

    const std::string label =
        (destName.find("noshadow") != std::string::npos) ? "no-shadow" : "full";
    const double cScaError = EstimateAngularIntegralRelativeError(
        thetaRows, m11Rows, cscaIntegral);
    const IntegralCharacteristics characteristics = ComputeIntegralCharacteristics(
        IntegralMethod::PhysicalOptics, incomingEnergy, cscaIntegral, cExtOt,
        hasExtinctionOt, hasAbsorption, cScaError);
    const std::string log = FormatIntegralCharacteristicsLog(
        characteristics, label);
    WriteIntegralCharacteristicsTsv(destName, label, characteristics);
    AppendIntegralCharacteristicsLog(destName, log);
    std::cerr << log;
    if (label == "full")
        integralSummary = log;
}

static bool AutoFullSymmetryErrorProjectionEnabled(double betaSym,
                                                   double gammaSym)
{
    const char *env = std::getenv("MBS_AUTOFULL_SYMMETRY_ERROR_PROJECT");
    if (env && *env)
        return !(env[0] == '0' && env[1] == '\0');
    (void)betaSym;
    (void)gammaSym;
    return false;
}

static bool AutoFullSymmetryOutputProjectionEnabled(double betaSym,
                                                    double gammaSym)
{
    const char *env = std::getenv("MBS_AUTOFULL_SYMMETRY_OUTPUT_PROJECT");
    if (!env || !*env)
        env = std::getenv("MBS_AUTOFULL_SYMMETRY_PROJECT");
    if (env && *env)
        return !(env[0] == '0' && env[1] == '\0');
    (void)betaSym;
    (void)gammaSym;
    return false;
}

static void ProjectRandomOrientationMuellerValues(double *v)
{
    static const int zeroElems[] = {2, 3, 6, 7, 8, 9, 12, 13};
    for (int idx : zeroElems)
        v[idx] = 0.0;
}

static void ProjectMuellerControlValues(std::vector<double> &values)
{
    const int rows = (int)values.size() / 16;
    for (int row = 0; row < rows; ++row)
        ProjectRandomOrientationMuellerValues(&values[row * 16]);
}

static void ProjectRandomOrientationMuellerSymmetry(matrix &m)
{
    // Full random-orientation hex-column averaging makes these components odd
    // under mirror/symmetry operations.  Oldauto's paired grid cancels them
    // exactly; Sobol/Owen only cancels statistically unless we project them.
    static const int zeroElems[][2] = {
        {0, 2}, {0, 3},
        {1, 2}, {1, 3},
        {2, 0}, {2, 1},
        {3, 0}, {3, 1}
    };
    for (const auto &rc : zeroElems)
        m[rc[0]][rc[1]] = 0.0;
}

static void ProjectRowsRandomOrientationMuellerSymmetry(
    std::vector<matrix> &rows)
{
    for (matrix &m : rows)
        ProjectRandomOrientationMuellerSymmetry(m);
}

static void ProjectSeedRowsRandomOrientationMuellerSymmetry(
    std::vector<std::vector<matrix>> &seedRows)
{
    for (std::vector<matrix> &rows : seedRows)
        ProjectRowsRandomOrientationMuellerSymmetry(rows);
}

static bool AutoFullRobustOwenOutputEnabled()
{
    return EnvFlagEnabled("MBS_AUTOFULL_ROBUST_OWEN_OUTPUT", false);
}

static double MedianValue(std::vector<double> values)
{
    if (values.empty())
        return 0.0;
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    if (n & 1u)
        return values[n / 2];
    return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

static void ReplaceRowsByRobustOwenMedian(
    std::vector<matrix> &rows,
    const std::vector<std::vector<matrix>> &seedRows)
{
    const int seedCount = (int)seedRows.size();
    if (seedCount < 3 || rows.empty())
        return;

    std::vector<double> values(seedCount);
    for (size_t t = 0; t < rows.size(); ++t)
    {
        matrix robust(4, 4);
        robust.Fill(0.0);
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
            {
                for (int s = 0; s < seedCount; ++s)
                    values[s] = seedRows[s][t][r][c] * seedCount;
                robust[r][c] = MedianValue(values);
            }
        rows[t] = robust;
    }
}

static int SharedFftGlobalRefineMaxRows()
{
    const char *env = std::getenv("MBS_FFT_GLOBAL_REFINE_MAX_ROWS");
    if (!env || !*env)
        return 8;
    char *end = nullptr;
    long value = std::strtol(env, &end, 10);
    if (!end || *end != '\0' || value < 0 || value > 128)
        return 8;
    return (int)value;
}

static bool SharedFftSpikeDebugEnabled()
{
    const char *env = std::getenv("MBS_FFT_SPIKE_DEBUG");
    return env && std::atoi(env) != 0;
}

static std::string SharedFftSpikeDebugPath()
{
    const char *env = std::getenv("MBS_FFT_SPIKE_DEBUG_FILE");
    if (env && *env)
        return std::string(env);
    return std::string("fft_spike_groups.csv");
}

static void AzimuthAverageOutputCell(const Arr2D &src, int nAz, int nZen,
                                     int t, double out[16])
{
    for (int k = 0; k < 16; ++k)
        out[k] = 0.0;

    const bool isPole = (t == 0 || t == nZen);
    for (int p = 0; p < nAz; ++p)
    {
        matrix m = src(p, t);
        if (isPole)
        {
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    out[r * 4 + c] += m[r][c];
            continue;
        }

        const double radAz = -p * (M_2PI / (double)nAz);
        matrix Lp(4, 4);
        Lp.Fill(0.0);
        for (int i = 0; i < 4; ++i)
            Lp[i][i] = 1.0;
        Lp[1][1] = cos(2 * radAz);
        Lp[1][2] = sin(2 * radAz);
        Lp[2][1] = -Lp[1][2];
        Lp[2][2] = Lp[1][1];
        matrix rotated = m * Lp;
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                out[r * 4 + c] += rotated[r][c];
    }

    const double inv = 1.0 / (double)nAz;
    for (int k = 0; k < 16; ++k)
        out[k] *= inv;
}

static std::vector<int> DetectGlobalFftSpikeRows(const Arr2D &m, int nAz,
                                                 int nZen, double threshold,
                                                 int maxRows)
{
    std::vector<int> rows;
    if (maxRows <= 0 || nAz <= 0 || nZen < 2)
        return rows;

    const int elems[] = {
        1 * 4 + 2, 1 * 4 + 3,
        2 * 4 + 1, 2 * 4 + 3,
        3 * 4 + 1, 3 * 4 + 2,
        0 * 4 + 2, 0 * 4 + 3
    };
    const int nElems = (int)(sizeof(elems) / sizeof(elems[0]));
    std::vector<std::pair<double, int>> scored;
    double prev[16], cur[16], next[16];

    for (int t = 1; t < nZen; ++t)
    {
        AzimuthAverageOutputCell(m, nAz, nZen, t - 1, prev);
        AzimuthAverageOutputCell(m, nAz, nZen, t, cur);
        AzimuthAverageOutputCell(m, nAz, nZen, t + 1, next);
        const double c11 = std::max(fabs(cur[0]), 1e-30);
        const double p11 = std::max(fabs(prev[0]), 1e-30);
        const double n11 = std::max(fabs(next[0]), 1e-30);
        double score = 0.0;
        for (int i = 0; i < nElems; ++i)
        {
            int e = elems[i];
            double yc = cur[e] / c11;
            double yn = 0.5 * (prev[e] / p11 + next[e] / n11);
            score = std::max(score, fabs(yc - yn));
        }
        if (score >= threshold)
            scored.push_back(std::make_pair(score, t));
    }

    std::sort(scored.begin(), scored.end(),
              [](const std::pair<double, int> &a,
                 const std::pair<double, int> &b) {
                  return a.first > b.first;
              });
    for (const auto &entry : scored)
    {
        rows.push_back(entry.second);
        if ((int)rows.size() >= maxRows)
            break;
    }
    std::sort(rows.begin(), rows.end());
    return rows;
}

static void WriteSpikeDebugHeader(std::ofstream &f)
{
    f << "size_index,label,theta_index,theta_deg,orient_start,orient_end,"
      << "beta_start_deg,beta_end_deg,gamma_start_deg,gamma_end_deg,"
      << "M11,M13,M14,M23,M24,M32,M42,"
      << "R13,R14,R23,R24,R32,R42\n";
}

static void WriteSpikeDebugRow(std::ofstream &f, size_t sizeIndex,
                               const std::string &label, int thetaIndex,
                               double thetaDeg, int orientStart, int orientEnd,
                               int nGamma, const AngleRange &betaRange,
                               const AngleRange &gammaRange,
                               const double avg[16])
{
    const int firstBeta = orientStart / nGamma;
    const int lastBeta = (orientEnd - 1) / nGamma;
    const int firstGamma = orientStart - firstBeta * nGamma;
    const int lastGamma = (orientEnd - 1) - lastBeta * nGamma;
    const double betaStart = betaRange.min + firstBeta * betaRange.step;
    const double betaEnd = betaRange.min + lastBeta * betaRange.step;
    const double gammaStart = gammaRange.min + (firstGamma + 0.5) * gammaRange.step;
    const double gammaEnd = gammaRange.min + (lastGamma + 0.5) * gammaRange.step;
    const double m11 = avg[0];
    const double denom = std::max(std::fabs(m11), 1e-300);
    f << sizeIndex << ','
      << label << ','
      << thetaIndex << ','
      << thetaDeg << ','
      << orientStart << ','
      << orientEnd << ','
      << RadToDeg(betaStart) << ','
      << RadToDeg(betaEnd) << ','
      << RadToDeg(gammaStart) << ','
      << RadToDeg(gammaEnd) << ','
      << avg[0] << ','
      << avg[2] << ','
      << avg[3] << ','
      << avg[6] << ','
      << avg[7] << ','
      << avg[9] << ','
      << avg[13] << ','
      << (avg[2] / denom) << ','
      << (avg[3] / denom) << ','
      << (avg[6] / denom) << ','
      << (avg[7] / denom) << ','
      << (avg[9] / denom) << ','
      << (avg[13] / denom) << '\n';
}

// ---- Checkpoint save/load for resume after crash ----
static void SaveCheckpoint(const std::string &path, const Arr2D &M, const Arr2D &M_ns,
                            double energy, double outputEnergy,
                            int completedOrient, int totalOrient,
                            unsigned long paramHash, int nAz, int nZen,
                            bool saveNoShadow)
{
    std::ofstream f(path, std::ios::binary);
    if (!f.is_open()) return;
    uint32_t magic = 0x4D425345; // "MBSE", checkpoint format v2.
    f.write((char*)&magic, 4);
    f.write((char*)&paramHash, sizeof(paramHash));
    f.write((char*)&completedOrient, sizeof(completedOrient));
    f.write((char*)&totalOrient, sizeof(totalOrient));
    f.write((char*)&energy, sizeof(energy));
    f.write((char*)&outputEnergy, sizeof(outputEnergy));
    f.write((char*)&nAz, 4);
    f.write((char*)&nZen, 4);
    const uint8_t hasNoShadow = saveNoShadow ? 1 : 0;
    f.write((char*)&hasNoShadow, 1);
    // Dump M data
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t < nZen; ++t) {
            matrix m = M(p, t);
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c) {
                    double v = m[r][c];
                    f.write((char*)&v, 8);
                }
        }
    if (saveNoShadow)
        for (int p = 0; p < nAz; ++p)
            for (int t = 0; t < nZen; ++t) {
                matrix m = M_ns(p, t);
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c) {
                        double v = m[r][c];
                        f.write((char*)&v, 8);
                    }
            }
    f.close();
}

static bool LoadCheckpoint(const std::string &path, Arr2D &M, Arr2D &M_ns,
                            double &energy, double &outputEnergy,
                            int &completedOrient, int totalOrient,
                            unsigned long paramHash, bool loadNoShadow)
{
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) return false;
    uint32_t magic; f.read((char*)&magic, 4);
    if (magic != 0x4D425345) return false;
    unsigned long storedHash; f.read((char*)&storedHash, sizeof(storedHash));
    if (storedHash != paramHash) {
        std::cerr << "Checkpoint param mismatch, ignoring" << std::endl;
        return false;
    }
    int storedCompleted, storedTotal;
    f.read((char*)&storedCompleted, sizeof(storedCompleted));
    f.read((char*)&storedTotal, sizeof(storedTotal));
    if (storedTotal != totalOrient) return false;
    f.read((char*)&energy, sizeof(energy));
    f.read((char*)&outputEnergy, sizeof(outputEnergy));
    int nAz, nZen; f.read((char*)&nAz, 4); f.read((char*)&nZen, 4);
    uint8_t hasNoShadow = 0;
    f.read((char*)&hasNoShadow, 1);
    if ((hasNoShadow != 0) != loadNoShadow)
        return false;
    // Load M
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t < nZen; ++t)
        {
            matrix m(4, 4);
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c) {
                    double v; f.read((char*)&v, 8);
                    m[r][c] = v;
                }
            M.replace(p, t, m);
        }
    if (loadNoShadow)
        for (int p = 0; p < nAz; ++p)
            for (int t = 0; t < nZen; ++t)
            {
                matrix m(4, 4);
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c) {
                        double v; f.read((char*)&v, 8);
                        m[r][c] = v;
                    }
                M_ns.replace(p, t, m);
            }
    completedOrient = storedCompleted;
    f.close();
    return true;
}

static unsigned long HashParams(double wave, double ri_re, double ri_im, int nActs, int nOrient,
                                 double L, double D, int nAz, int nZen)
{
    unsigned long h = 0;
    auto mix = [&](double v) { h ^= std::hash<double>{}(v) + 0x9e3779b9 + (h<<6) + (h>>2); };
    auto mixi = [&](int v) { h ^= std::hash<int>{}(v) + 0x9e3779b9 + (h<<6) + (h>>2); };
    mix(wave); mix(ri_re); mix(ri_im); mixi(nActs); mixi(nOrient); mix(L); mix(D); mixi(nAz); mixi(nZen);
    return h;
}

static void MixHashDouble(unsigned long &h, double v)
{
    h ^= std::hash<double>{}(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
}

static void MixHashInt(unsigned long &h, int v)
{
    h ^= std::hash<int>{}(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
}

static void SaveOldautoCheckpoint(const std::string &path,
                                  const Arr2D &M, const Arr2D &M_ns,
                                  double energy, double outputEnergy,
                                  double extinctionOt, bool hasExtinctionOt,
                                  int completedBeta, int nBeta, int nGamma,
                                  unsigned long paramHash, int nAz, int nZen,
                                  bool saveNoShadow)
{
    std::string tmp = path + ".tmp";
    std::ofstream f(tmp, std::ios::binary);
    if (!f.is_open()) return;
    uint32_t magic = 0x4D425341; // "MBSA" oldauto checkpoint
    uint32_t version = 2;
    int hasExt = hasExtinctionOt ? 1 : 0;
    int hasNoShadow = saveNoShadow ? 1 : 0;
    int completedOrient = completedBeta * nGamma;
    int totalOrient = nBeta * nGamma;
    f.write((char*)&magic, 4);
    f.write((char*)&version, 4);
    f.write((char*)&paramHash, sizeof(paramHash));
    f.write((char*)&completedBeta, sizeof(completedBeta));
    f.write((char*)&nBeta, sizeof(nBeta));
    f.write((char*)&nGamma, sizeof(nGamma));
    f.write((char*)&completedOrient, sizeof(completedOrient));
    f.write((char*)&totalOrient, sizeof(totalOrient));
    f.write((char*)&energy, sizeof(energy));
    f.write((char*)&outputEnergy, sizeof(outputEnergy));
    f.write((char*)&extinctionOt, sizeof(extinctionOt));
    f.write((char*)&hasExt, sizeof(hasExt));
    f.write((char*)&hasNoShadow, sizeof(hasNoShadow));
    f.write((char*)&nAz, 4);
    f.write((char*)&nZen, 4);
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t < nZen; ++t) {
            matrix m = M(p, t);
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c) {
                    double v = m[r][c];
                    f.write((char*)&v, 8);
                }
        }
    if (saveNoShadow)
        for (int p = 0; p < nAz; ++p)
            for (int t = 0; t < nZen; ++t) {
                matrix m = M_ns(p, t);
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c) {
                        double v = m[r][c];
                        f.write((char*)&v, 8);
                    }
            }
    f.close();
    if (f.good())
        std::rename(tmp.c_str(), path.c_str());
    else
        std::remove(tmp.c_str());
}

static bool LoadOldautoCheckpoint(const std::string &path,
                                  Arr2D &M, Arr2D &M_ns,
                                  double &energy, double &outputEnergy,
                                  double &extinctionOt, bool &hasExtinctionOt,
                                  int &completedBeta, int nBeta, int nGamma,
                                  unsigned long paramHash)
{
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) return false;
    uint32_t magic = 0, version = 0;
    f.read((char*)&magic, 4);
    f.read((char*)&version, 4);
    if (magic != 0x4D425341 || version != 2) return false;
    unsigned long storedHash = 0;
    f.read((char*)&storedHash, sizeof(storedHash));
    if (storedHash != paramHash) {
        std::cerr << "Oldauto checkpoint param mismatch, ignoring" << std::endl;
        return false;
    }
    int storedCompletedBeta = 0, storedBeta = 0, storedGamma = 0;
    int completedOrient = 0, totalOrient = 0;
    f.read((char*)&storedCompletedBeta, sizeof(storedCompletedBeta));
    f.read((char*)&storedBeta, sizeof(storedBeta));
    f.read((char*)&storedGamma, sizeof(storedGamma));
    f.read((char*)&completedOrient, sizeof(completedOrient));
    f.read((char*)&totalOrient, sizeof(totalOrient));
    if (storedBeta != nBeta || storedGamma != nGamma || totalOrient != nBeta * nGamma)
        return false;
    f.read((char*)&energy, sizeof(energy));
    f.read((char*)&outputEnergy, sizeof(outputEnergy));
    f.read((char*)&extinctionOt, sizeof(extinctionOt));
    int hasExt = 0;
    f.read((char*)&hasExt, sizeof(hasExt));
    int hasNoShadow = 0;
    f.read((char*)&hasNoShadow, sizeof(hasNoShadow));
    int nAz = 0, nZen = 0;
    f.read((char*)&nAz, 4);
    f.read((char*)&nZen, 4);
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t < nZen; ++t)
        {
            matrix m(4, 4);
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c) {
                    double v = 0.0;
                    f.read((char*)&v, 8);
                    m[r][c] = v;
                }
            M.replace(p, t, m);
        }
    const bool loadNoShadow = hasNoShadow != 0 && StrArr(M_ns) > 0 && ColArr(M_ns) > 0;
    if (hasNoShadow)
        for (int p = 0; p < nAz; ++p)
            for (int t = 0; t < nZen; ++t)
            {
                matrix m(4, 4);
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c) {
                        double v = 0.0;
                        f.read((char*)&v, 8);
                        m[r][c] = v;
                    }
                if (loadNoShadow)
                    M_ns.replace(p, t, m);
            }
    if (!f.good()) return false;
    completedBeta = storedCompletedBeta;
    hasExtinctionOt = hasExt != 0;
    return completedBeta >= 0 && completedBeta <= nBeta
        && completedOrient == completedBeta * nGamma;
}
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace std;
using ::complex;

// Helper: pack Arr2D (N x M grid of n x m matrices) into flat double array
static void Arr2DToFlat(const Arr2D &arr, int N, int M, double *buf)
{
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j) {
            matrix m = arr(i, j);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    buf[((i*M + j)*4 + a)*4 + b] = m[a][b];
        }
}

// Helper: unpack flat double array into Arr2D (replace mode)
static void FlatToArr2D(const double *buf, int N, int M, Arr2D &arr)
{
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j) {
            matrix mt(4, 4);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    mt[a][b] = buf[((i*M + j)*4 + a)*4 + b];
            arr.replace(i, j, mt);
        }
}

// Helper: MPI reduce Arr2D + energy counters, rank 0 gets result
static void MPI_ReduceMueller(HandlerPO *hp, int nAz, int nZen,
                               double &incomingEnergy, int mpi_rank)
{
#ifdef USE_MPI
    int totalDoubles = (nAz+1) * (nZen+1) * 16;
    std::vector<double> sendbuf(totalDoubles), recvbuf(totalDoubles);

    // Reduce M
    Arr2DToFlat(hp->M, nAz+1, nZen+1, sendbuf.data());
    MPI_Reduce(sendbuf.data(), recvbuf.data(), totalDoubles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) FlatToArr2D(recvbuf.data(), nAz+1, nZen+1, hp->M);

    if (hp->ComputeNoShadow())
    {
        Arr2DToFlat(hp->M_noshadow, nAz+1, nZen+1, sendbuf.data());
        MPI_Reduce(sendbuf.data(), recvbuf.data(), totalDoubles, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (mpi_rank == 0) FlatToArr2D(recvbuf.data(), nAz+1, nZen+1, hp->M_noshadow);
    }

    // Reduce incomingEnergy
    double totalEnergy = 0;
    MPI_Reduce(&incomingEnergy, &totalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) incomingEnergy = totalEnergy;

    // Reduce output energy used for absorption/extinction accounting.
    double totalOutputEnergy = 0;
    MPI_Reduce(&hp->m_outputEnergy, &totalOutputEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) hp->m_outputEnergy = totalOutputEnergy;

    double totalExtOt = 0;
    MPI_Reduce(&hp->m_extinctionCrossSectionOt, &totalExtOt, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0)
    {
        hp->m_extinctionCrossSectionOt = totalExtOt;
        hp->m_hasExtinctionOt = true;
    }
#endif
}

static PreparedOrientation ScalePreparedOrientation(const PreparedOrientation &src,
                                                    double scale,
                                                    double waveIndex,
                                                    double cAbs)
{
    PreparedOrientation dst = src;
    const double scale2 = scale * scale;
    for (PreparedBeam &pb : dst.beams)
    {
        for (int e = 0; e < pb.edgeData.nVertices; ++e)
        {
            pb.edgeData.x[e] *= scale;
            pb.edgeData.y[e] *= scale;
            pb.edgeData.intercept_x[e] *= scale;
            pb.edgeData.intercept_y[e] *= scale;
        }

        pb.info.area *= scale2;
        pb.info.projLenght *= scale;
        pb.info.center.x *= scale;
        pb.info.center.y *= scale;
        pb.info.center.z *= scale;
        pb.info.projectedCenter.x *= scale;
        pb.info.projectedCenter.y *= scale;
        pb.info.projectedCenter.z *= scale;
        for (double &v : pb.info.opticalLengths)
            v *= scale;

        pb.cenx *= scale;
        pb.ceny *= scale;
        pb.cenz *= scale;
        pb.beam_area *= scale2;
        pb.origBeam.opticalPath *= scale;

        ::complex absorption(1.0, 0.0);
        if (!pb.absorptionPaths.empty())
        {
            double sum = 0.0;
            int count = 0;
            for (double path : pb.absorptionPaths)
            {
                sum += (path > DBL_EPSILON) ? std::exp(cAbs * path * scale) : 1.0;
                ++count;
            }
            if (count > 0)
                absorption = sum / count;
        }

        matrixC J_base = pb.origBeam.J * absorption;
        matrixC J_phased = pb.isExternal
            ? pb.origBeam.J * exp_im(waveIndex * pb.origBeam.opticalPath)
            : J_base * exp_im(waveIndex * pb.info.projLenght);
        if (pb.isExternal)
            J_phased *= -1.0;
        if (!pb.isExternal && (pb.origBeam.nActs & 1))
            J_phased *= -1.0;

        ::complex jp00 = J_phased[0][0], jp01 = J_phased[0][1];
        ::complex jp10 = J_phased[1][0], jp11 = J_phased[1][1];
        pb.jp00r = real(jp00); pb.jp00i = imag(jp00);
        pb.jp01r = real(jp01); pb.jp01i = imag(jp01);
        pb.jp10r = real(jp10); pb.jp10i = imag(jp10);
        pb.jp11r = real(jp11); pb.jp11i = imag(jp11);
    }
    return dst;
}

static double PreparedAbsorptionFactor(const PreparedBeam &pb, double scale,
                                       double cAbs)
{
    if (pb.absorptionPaths.empty())
        return 1.0;
    double sum = 0.0;
    int count = 0;
    for (double path : pb.absorptionPaths)
    {
        sum += (path > DBL_EPSILON) ? std::exp(cAbs * path * scale) : 1.0;
        ++count;
    }
    return count > 0 ? sum / count : 1.0;
}

static double PreparedOutputEnergy(const PreparedOrientation &po, double scale,
                                   double cAbs)
{
    const double scale2 = scale * scale;
    double energy = 0.0;
    for (const PreparedBeam &pb : po.beams)
    {
        if (pb.outputCrossSection <= 0.0 || pb.outputMueller00 == 0.0)
            continue;
        const double absorption = PreparedAbsorptionFactor(pb, scale, cAbs);
        energy += pb.outputCrossSection * scale2 * pb.outputMueller00
            * absorption * absorption * po.sinZenith;
    }
    return energy;
}

TracerPOTotal::TracerPOTotal(Particle *particle, int nActs,
                             const string &resultFileName)
    : TracerPO(particle, nActs, resultFileName)
{
}

void TracerPOTotal::ResetConvergenceReport(const std::string &mode, double eps)
{
    m_convergenceTarget = eps;
    if (m_mpiRank != 0)
        return;
    const std::string path = m_resultDirName + "_convergence.tsv";
    std::ofstream out(path.c_str(), std::ios::out);
    out << "# adaptive_mode\t" << mode << "\n";
    if (eps > 0.0)
        out << "# requested_relative_error\t" << std::setprecision(12)
            << eps << "\n";
    else
        out << "# requested_relative_error\tmixed; see target_relative_error column\n";
    out << "parameter\tsweep\tprimary\tsecondary\ttarget_relative_error"
        << "\tmax_relative_error"
        << "\tworst_theta_deg\tworst_element\tstatus\n";
}

void TracerPOTotal::RecordConvergenceStep(
    const std::string &parameter, int sweep, int primary, int secondary,
    double error, double thetaDeg, int element, const std::string &status)
{
    if (m_mpiRank != 0)
        return;
    const std::string path = m_resultDirName + "_convergence.tsv";
    std::ofstream out(path.c_str(), std::ios::out | std::ios::app);
    out << parameter << '\t' << sweep << '\t' << primary << '\t'
        << secondary << '\t' << std::setprecision(12)
        << m_convergenceTarget << '\t' << error << '\t'
        << thetaDeg << '\t' << element << '\t' << status << '\n';
}

void TracerPOTotal::TraceRandom(const AngleRange &betaRange,
                                const AngleRange &gammaRange)
{
    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) {
        std::cerr << "Error: handler is not HandlerPO in TraceRandom" << std::endl;
        return;
    }

    m_incomingEnergy = 0;
    handlerPO->m_outputEnergy = 0;
    handlerPO->m_extinctionCrossSectionOt = 0;
    handlerPO->m_hasExtinctionOt = false;
    handlerPO->M.ClearArr();
    handlerPO->M_noshadow.ClearArr();
    handlerPO->CleanJ();

    int nGamma = gammaRange.number; // MBS-raw: same as gammaRange.number
    const bool betaMidpoint = OldautoBetaMidpointEnabled() && !m_fastPoleGamma;
    int nBeta = OldautoBetaCount(betaRange, betaMidpoint);
    int nOrientations = nBeta * nGamma;
    int betaNorm = (m_symmetry.beta < M_PI_2+FLT_EPSILON && m_symmetry.beta > M_PI_2-FLT_EPSILON) ? 1 : 2;
    double normGamma = gammaRange.number * betaNorm; // MBS-raw: normalize by gammaRange.number
    const bool gammaStagger = OldautoGammaStaggerEnabled();

    // --coh_orient: legacy mode — HandleBeams accumulates Jones coherently
    // across ALL orientations, then single AddToMueller at end.
    // No chunking, no OpenMP (matches old MBS-raw exactly).
    if (m_cohOrient) {
        m_handler->SetNormIndex(normGamma);
        vector<Beam> outBeams;
        for (int i = 0; i < nBeta; ++i) {
            double beta = OldautoBetaAngle(betaRange, i, betaMidpoint);
            double dcos;
            CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normGamma, dcos);
            const bool pole = (fabs(beta) <= FLT_EPSILON
                               || fabs(beta - M_PI) <= FLT_EPSILON);
            const bool fastPole = pole && m_fastPoleGamma;
            const int gammaCount = fastPole ? 1 : gammaRange.number;
            const double gammaWeight = fastPole ? dcos * gammaRange.number : dcos;
            m_handler->SetSinZenith(gammaWeight);
            for (int j = 0; j < gammaCount; ++j) {
                double gamma = OldautoGammaAngle(gammaRange, nGamma, j, i,
                                                 gammaStagger && !fastPole);
                m_particle->Rotate(beta, gamma, 0);
                if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
                m_scattering->ScatterLight(0, 0, outBeams);
                m_handler->HandleBeams(outBeams, gammaWeight);
                m_incomingEnergy += m_scattering->GetIncedentEnergy() * gammaWeight;
                outBeams.clear();
            }
        }
        m_handler->WriteTotalMatricesToFile(m_resultDirName);
        m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
        CalcTimer timer; timer.Start();
        OutputStatisticsPO(timer, nOrientations, m_resultDirName);
        return;
    }

    // Default: incoherent per-orientation (physically correct)
    // Loop over beta explicitly — enables per-beta Mueller saving

    if (m_mpiRank == 0)
        std::cout << "Random grid: " << nBeta << " x " << nGamma
                  << " = " << nOrientations << " orientations" << std::endl;

    CalcTimer timer;
    timer.Start();
    if (m_mpiRank == 0) OutputStartTime(timer);

    m_handler->SetNormIndex(normGamma);

    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;
    long long count = 0;
    bool oldautoCheckpoint = m_enableCheckpoint;
    if (oldautoCheckpoint && m_mpiSize > 1)
    {
        if (m_mpiRank == 0)
            std::cerr << "WARNING: --checkpoint for oldauto/random is disabled under MPI; "
                      << "use single-process runs or --save_betas for per-rank diagnostics."
                      << std::endl;
        oldautoCheckpoint = false;
    }
    std::string oldautoCkptPath = m_resultDirName + "_oldauto_checkpoint.bin";
    const ::complex ri = m_particle->GetRefractiveIndex();
    unsigned long oldautoParamHash = HashParams(m_scattering->m_wave,
        real(ri), imag(ri), m_scattering->GetMaxReflections(), nOrientations,
        m_particle->MaximalDimention(), 0, nAz + 1, nZen + 1);
    MixHashDouble(oldautoParamHash, betaRange.min);
    MixHashDouble(oldautoParamHash, betaRange.max);
    MixHashDouble(oldautoParamHash, betaRange.step);
    MixHashInt(oldautoParamHash, betaRange.number);
    MixHashDouble(oldautoParamHash, gammaRange.min);
    MixHashDouble(oldautoParamHash, gammaRange.max);
    MixHashDouble(oldautoParamHash, gammaRange.step);
    MixHashInt(oldautoParamHash, gammaRange.number);
    MixHashInt(oldautoParamHash, betaMidpoint ? 1 : 0);
    MixHashInt(oldautoParamHash, gammaStagger ? 1 : 0);
    MixHashInt(oldautoParamHash, m_mirrorGamma ? 1 : 0);
    MixHashInt(oldautoParamHash, shadowOff ? 1 : 0);
    MixHashInt(oldautoParamHash, handlerPO->IsGpuEnabled() ? 1 : 0);
    MixHashInt(oldautoParamHash, handlerPO->IsFftEnabled() ? 1 : 0);
    int resumeBeta = 0;
    if (oldautoCheckpoint && m_mpiRank == 0)
    {
        if (LoadOldautoCheckpoint(oldautoCkptPath,
                                  handlerPO->M, handlerPO->M_noshadow,
                                  m_incomingEnergy, handlerPO->m_outputEnergy,
                                  handlerPO->m_extinctionCrossSectionOt,
                                  handlerPO->m_hasExtinctionOt,
                                  resumeBeta, nBeta, nGamma,
                                  oldautoParamHash))
        {
            count = (long long)resumeBeta * nGamma;
            std::ostringstream line;
            line << "*** RESUMED oldauto checkpoint: beta " << resumeBeta
                 << "/" << nBeta << ", orientations " << count
                 << "/" << nOrientations;
            std::cout << line.str() << std::endl;
            AppendTextLog(line.str() + "\n");
        }
    }
    int nThreads = 1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        nThreads = omp_get_num_threads();
    }
#endif
    int parallelTraceMinGamma = 2;
    const char *traceMinEnv = std::getenv("MBS_TRACE_MIN_GAMMA");
    if (traceMinEnv && *traceMinEnv)
    {
        char *end = nullptr;
        long value = std::strtol(traceMinEnv, &end, 10);
        if (end && *end == '\0' && value > 0)
            parallelTraceMinGamma = (int)value;
    }
    if (m_mpiRank == 0)
    {
        std::ostringstream line;
        line << "Oldauto/random Phase 1: ";
        if (nThreads > 1)
        {
            line << "gamma blocks < " << parallelTraceMinGamma
                 << " trace sequentially; larger blocks trace in parallel ("
                 << nThreads << " threads)";
        }
        else
        {
            line << "sequential tracing (1 thread)";
        }
        if (gammaStagger)
            line << "; gamma stagger enabled";
        if (betaMidpoint)
            line << "; beta midpoint enabled";
        else if (m_fastPoleGamma)
            line << "; beta endpoints enabled for --pole";
        std::cout << line.str() << std::endl;
        AppendTextLog(line.str() + "\n");
    }

    // --save_betas: create output directory
    std::string betaDir;
    if (m_saveBetas && m_mpiRank == 0) {
        betaDir = m_resultDirName + "_betas";
        mkdir(betaDir.c_str(), 0755);
    }

    double phase1_total = 0, phase2_total = 0;

    struct BetaBlock
    {
        int ib = 0;
        int gammaStart = 0;
        int gammaCount = 0;
        int gammaFullCount = 0;
        double beta = 0.0;
        double traceBeta = 0.0;
        double dcos = 0.0;
        double phase1 = 0.0;
        std::vector<PreparedOrientation> prepared;
        std::vector<double> energies;
        std::vector<double> outputEnergies;
        std::vector<int> beamCounts;
    };

    int gammaChunk = nGamma;
    const double effectiveJCutoff =
        (handlerPO->m_beamCutoffJRel >= 0) ? handlerPO->m_beamCutoffJRel
                                           : handlerPO->m_targetEps;
    const double effectiveAreaCutoff =
        (handlerPO->m_beamCutoffAreaRel >= 0) ? handlerPO->m_beamCutoffAreaRel
                                             : handlerPO->m_targetEps;
    const bool noBeamCutoff =
        effectiveJCutoff <= 0.0
        && effectiveAreaCutoff <= 0.0
        && handlerPO->m_beamCutoffImportanceRel <= 0.0
        && handlerPO->m_beamCutoff <= 0.0;
    bool autoGammaChunk = false;
    bool hostMemGuarded = false;
    bool hostMemClamped = false;
    bool hostMemConservative = false;
    long long hostMemTotalKb = 0;
    long long hostMemAvailableKb = 0;
    long long hostMemRssKb = 0;
    long long hostMemBudgetKb = 0;
    long long hostMemGridKb = 0;
    if (m_sobolChunkSize > 0)
    {
        gammaChunk = std::max(1, std::min(nGamma, m_sobolChunkSize));
    }
    else if (handlerPO->IsGpuEnabled())
    {
        int defaultChunk = nGamma;
        const char *chunkEnv = std::getenv("MBS_OLDAUTO_GAMMA_CHUNK");
        if (chunkEnv && *chunkEnv)
        {
            char *end = nullptr;
            long value = std::strtol(chunkEnv, &end, 10);
            if (end && *end == '\0' && value > 0)
                defaultChunk = (int)value;
        }
        gammaChunk = std::max(1, std::min(nGamma, defaultChunk));
        autoGammaChunk = gammaChunk < nGamma;
    }

    if (handlerPO->IsGpuEnabled())
    {
        const int requestedGammaChunk = std::max(1, std::min(nGamma, gammaChunk));
        gammaChunk = HostMemoryAwareGammaChunk(
            requestedGammaChunk,
            noBeamCutoff,
            nAz,
            nZen,
            handlerPO->ComputeNoShadow(),
            handlerPO->isCoh,
            handlerPO->IsGpuEnabled(),
            &hostMemTotalKb,
            &hostMemAvailableKb,
            &hostMemRssKb,
            &hostMemBudgetKb,
            &hostMemGridKb);
        gammaChunk = std::max(1, std::min(nGamma, gammaChunk));
        hostMemGuarded = true;
        hostMemClamped = gammaChunk < requestedGammaChunk;
        hostMemConservative = noBeamCutoff;
    }

    if (m_mpiRank == 0 && gammaChunk < nGamma)
    {
        std::ostringstream line;
        line << "Oldauto/random memory: gamma chunk=" << gammaChunk
             << "/" << nGamma << " per beta";
        if (hostMemGuarded)
        {
            line << " (GPU host-RAM guard";
            if (hostMemTotalKb > 0)
            {
                line << ", MemAvailable=" << (hostMemAvailableKb / 1024)
                     << " MB, VmRSS=" << (hostMemRssKb / 1024)
                     << " MB, grid-transient=" << (hostMemGridKb / 1024)
                     << " MB, budget=" << (hostMemBudgetKb / 1024)
                     << " MB";
                if (HostMemoryBudgetOverrideKb() > 0)
                    line << ", shared-budget=" << (HostMemoryBudgetOverrideKb() / 1024)
                         << " MB";
                if (hostMemRssKb + hostMemGridKb > hostMemBudgetKb)
                    line << ", grid dominates budget";
            }
            if (hostMemConservative)
                line << ", no-cutoff conservative mode";
            if (hostMemClamped)
                line << ", clamped requested chunk";
            else if (autoGammaChunk)
                line << ", auto chunk";
            line << "; tune with MBS_HOST_MEM_FRACTION/MBS_HOST_MEM_RESERVE_MB)";
        }
        else
            line << " (--chunk)";
        std::cout << line.str() << std::endl;
        AppendTextLog(line.str() + "\n");
    }
    m_scattering->PrepareForParallelTrace();

    auto prepareBetaBlock = [&](int ib, int gammaStart, int gammaLimit) {
        BetaBlock block;
        block.ib = ib;
        block.gammaStart = gammaStart;
        block.beta = OldautoBetaAngle(betaRange, ib, betaMidpoint);
        CalcCsBeta(betaNorm, block.beta, betaRange, gammaRange, normGamma, block.dcos);
        const bool pole = (fabs(block.beta) <= FLT_EPSILON
                           || fabs(block.beta - M_PI) <= FLT_EPSILON);
        const bool fastPole = pole && m_fastPoleGamma;
        block.gammaFullCount = fastPole ? 1 : nGamma;
        block.gammaCount = fastPole ? 1 : std::min(gammaLimit, nGamma - gammaStart);
        const double gammaWeight = fastPole ? block.dcos * nGamma : block.dcos;
        block.traceBeta = block.beta;
        if (pole && !fastPole && betaRange.step > 0.0)
        {
            if (fabs(block.beta) <= FLT_EPSILON)
                block.traceBeta = block.beta + 0.5 * betaRange.step;
            else
                block.traceBeta = block.beta - 0.5 * betaRange.step;
        }

        auto tp1 = std::chrono::high_resolution_clock::now();
        block.prepared.assign(block.gammaCount, PreparedOrientation());
        block.energies.assign(block.gammaCount, 0.0);
        block.outputEnergies.assign(block.gammaCount, 0.0);
        block.beamCounts.assign(block.gammaCount, 0);

        if (nThreads > 1 && block.gammaCount >= parallelTraceMinGamma)
        {
            ParallelExceptionState parallelError;
            #pragma omp parallel
            {
                Particle localParticle = *m_particle;
                Scattering *localScatter = m_scattering->CloneFor(&localParticle, &m_incidentLight);

                HandlerPO localHandler(&localParticle, &m_incidentLight,
                                       handlerPO->nTheta, m_scattering->m_wave);
                localHandler.ConfigureForThreadLocalPrepare(*handlerPO, localScatter);

                std::vector<Beam> localBeams;

                #pragma omp for schedule(dynamic, 1)
                for (int jj = 0; jj < block.gammaCount; ++jj)
                {
                    if (parallelError.Failed())
                        continue;
                    try
                    {
                    int globalGamma = fastPole ? 0 : (block.gammaStart + jj);
                    double gamma = OldautoGammaAngle(gammaRange, nGamma,
                                                     globalGamma, ib,
                                                     gammaStagger && !fastPole);
                    localParticle.Rotate(block.traceBeta, gamma, 0);
                    if (!shadowOff) localScatter->FormShadowBeam(localBeams);
                    bool ok = localScatter->ScatterLight(0, 0, localBeams);
                    block.beamCounts[jj] = (int)localBeams.size();
                    if (ok)
                    {
                        double beforeOutput = localHandler.m_outputEnergy;
                        localHandler.PrepareBeams(localBeams, gammaWeight, block.prepared[jj]);
                        block.outputEnergies[jj] = localHandler.m_outputEnergy - beforeOutput;
                    }
                    else
                    {
                        block.prepared[jj].sinZenith = gammaWeight;
                    }
                    block.energies[jj] = localScatter->GetIncedentEnergy() * gammaWeight;
                    localBeams.clear();
                    }
                    catch (...)
                    {
                        parallelError.Capture();
                    }
                }

                delete localScatter;
            }
            parallelError.Rethrow();
        }
        else
        {
            Particle localParticle = *m_particle;
            Scattering *localScatter = m_scattering->CloneFor(&localParticle, &m_incidentLight);

            HandlerPO localHandler(&localParticle, &m_incidentLight,
                                   handlerPO->nTheta, m_scattering->m_wave);
            localHandler.ConfigureForThreadLocalPrepare(*handlerPO, localScatter);

            std::vector<Beam> localBeams;
            for (int jj = 0; jj < block.gammaCount; ++jj)
            {
                int globalGamma = fastPole ? 0 : (block.gammaStart + jj);
                double gamma = OldautoGammaAngle(gammaRange, nGamma,
                                                 globalGamma, ib,
                                                 gammaStagger && !fastPole);
                localParticle.Rotate(block.traceBeta, gamma, 0);
                if (!shadowOff) localScatter->FormShadowBeam(localBeams);
                bool ok = localScatter->ScatterLight(0, 0, localBeams);
                block.beamCounts[jj] = (int)localBeams.size();
                if (ok)
                {
                    double beforeOutput = localHandler.m_outputEnergy;
                    localHandler.PrepareBeams(localBeams, gammaWeight, block.prepared[jj]);
                    block.outputEnergies[jj] = localHandler.m_outputEnergy - beforeOutput;
                }
                else
                {
                    block.prepared[jj].sinZenith = gammaWeight;
                }
                block.energies[jj] = localScatter->GetIncedentEnergy() * gammaWeight;
                localBeams.clear();
            }

            delete localScatter;
        }

        block.phase1 = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - tp1).count();
        return block;
    };

    int myBetaStart = m_mpiRank * nBeta / m_mpiSize;
    int myBetaEnd = (m_mpiRank + 1) * nBeta / m_mpiSize;
    if (m_mpiRank == 0 && m_mpiSize > 1)
    {
        std::ostringstream line;
        line << "MPI oldauto/random: beta split across " << m_mpiSize
             << " ranks, rank0 beta range [" << myBetaStart << ", "
             << myBetaEnd << ")";
        std::cout << line.str() << std::endl;
        AppendTextLog(line.str() + "\n");
    }

    // Process beta-by-beta.  Under MPI each rank owns a disjoint beta range.
    // Inside a beta, gamma is streamed in chunks to bound per-rank memory.
    for (int ib = 0; ib < nBeta; ++ib)
    {
        if (ib < resumeBeta)
            continue;
        if (ib < myBetaStart || ib >= myBetaEnd)
            continue;

        double beta = OldautoBetaAngle(betaRange, ib, betaMidpoint);
        double dcos = 0.0;
        CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normGamma, dcos);
        const bool pole = (fabs(beta) <= FLT_EPSILON
                           || fabs(beta - M_PI) <= FLT_EPSILON);
        const bool fastPole = pole && m_fastPoleGamma;
        const int gammaFullCount = fastPole ? 1 : nGamma;

        const bool computeNoShadow = handlerPO->ComputeNoShadow();
        Arr2D betaM(nAz+1, nZen+1, 4, 4); betaM.ClearArr();
        Arr2D betaM_ns(computeNoShadow ? nAz+1 : 0,
                       computeNoShadow ? nZen+1 : 0,
                       computeNoShadow ? 4 : 0,
                       computeNoShadow ? 4 : 0);
        if (computeNoShadow)
            betaM_ns.ClearArr();

        auto processBlock = [&](BetaBlock &block)
        {
            const int gammaCount = block.gammaCount;
            std::vector<PreparedOrientation> &chunkPrepared = block.prepared;

            for (int jj = 0; jj < gammaCount; ++jj)
            {
                int globalGamma = fastPole ? 0 : (block.gammaStart + jj);
                m_incomingEnergy += block.energies[jj];
                handlerPO->m_outputEnergy += block.outputEnergies[jj];
                handlerPO->m_extinctionCrossSectionOt +=
                    chunkPrepared[jj].extinctionOt;
                handlerPO->m_hasExtinctionOt = true;
                ++count;
                if (m_mpiRank == 0)
                    OutputProgress(nOrientations, count, ib*nGamma+globalGamma, 0,
                                   timer, block.beamCounts[jj]);
            }
            phase1_total += block.phase1;

            // Phase 2: diffraction (OpenMP parallel over this gamma chunk)
            auto tp2 = std::chrono::high_resolution_clock::now();

            if (handlerPO->IsGpuEnabled())
            {
                Arr2D localM(nAz+1, nZen+1, 4, 4); localM.ClearArr();
                Arr2D localM_ns(computeNoShadow ? nAz+1 : 0,
                                computeNoShadow ? nZen+1 : 0,
                                computeNoShadow ? 4 : 0,
                                computeNoShadow ? 4 : 0);
                if (computeNoShadow)
                    localM_ns.ClearArr();
                for (int gpuStart = 0; gpuStart < gammaCount; )
                {
                    int gpuBatchSize = handlerPO->SelectGpuOrientationBatchSize(
                        chunkPrepared, gpuStart, gammaCount - gpuStart);
                    int gpuEnd = std::min(gpuStart + gpuBatchSize, gammaCount);
                    bool ok = handlerPO->IsFftEnabled()
                        ? handlerPO->HandleOrientationsToLocalGpuFftPhi(
                            chunkPrepared, gpuStart, gpuEnd - gpuStart,
                            localM, localM_ns)
                        : handlerPO->HandleOrientationsToLocalGpu(
                            chunkPrepared, gpuStart, gpuEnd - gpuStart,
                            localM, localM_ns);
                    if (!ok)
                    {
                        std::cerr << "ERROR: --gpu requested but GPU diffraction backend "
                                  << "could not process oldauto/random beta block." << std::endl;
                        throw std::runtime_error("GPU diffraction backend failed");
                    }
                    gpuStart = gpuEnd;
                }
                for (int p=0;p<nAz;++p) for (int t=0;t<=nZen;++t) {
                    betaM.insert(p,t,localM(p,t));
                    if (computeNoShadow)
                        betaM_ns.insert(p,t,localM_ns(p,t));
                }
            }
            else
            {
                #pragma omp parallel
                {
                    Arr2D localM(nAz+1, nZen+1, 4, 4); localM.ClearArr();
                    Arr2D localM_ns(computeNoShadow ? nAz+1 : 0,
                                    computeNoShadow ? nZen+1 : 0,
                                    computeNoShadow ? 4 : 0,
                                    computeNoShadow ? 4 : 0);
                    if (computeNoShadow)
                        localM_ns.ClearArr();
                    std::vector<Arr2DC> localJ, localJ_ns;
                    if (handlerPO->isCoh) {
                        Arr2DC t1(nAz+1,nZen+1,2,2); t1.ClearArr(); localJ.push_back(t1);
                        if (computeNoShadow) {
                            Arr2DC t2(nAz+1,nZen+1,2,2); t2.ClearArr(); localJ_ns.push_back(t2);
                        }
                    }
                    #pragma omp for schedule(dynamic, 1)
                    for (int i = 0; i < gammaCount; ++i) {
                        if (!chunkPrepared[i].beams.empty())
                            handlerPO->HandleBeamsToLocal(chunkPrepared[i], localM, localJ,
                                                           (handlerPO->isCoh && computeNoShadow) ? &localJ_ns : nullptr);
                        if (handlerPO->isCoh && !localJ.empty()) {
                            double w = chunkPrepared[i].sinZenith;
                            HandlerPO::AddToMuellerLocal(localJ, w, localM, nAz, nZen);
                            if (computeNoShadow)
                                HandlerPO::AddToMuellerLocal(localJ_ns, w, localM_ns, nAz, nZen);
                            localJ[0].ClearArr();
                            if (computeNoShadow)
                                localJ_ns[0].ClearArr();
                        }
                    }
                    #pragma omp critical
                    { for (int p=0;p<nAz;++p) for (int t=0;t<=nZen;++t) {
                        betaM.insert(p,t,localM(p,t));
                        if (computeNoShadow)
                            betaM_ns.insert(p,t,localM_ns(p,t));
                    } }
                }
            }

            phase2_total += std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - tp2).count();
            chunkPrepared.clear(); chunkPrepared.shrink_to_fit();
        };

        const bool pipelineTraceGpu = handlerPO->IsGpuEnabled()
            && !handlerPO->IsFftEnabled()
            && gammaChunk < gammaFullCount
            && !hostMemConservative;
        if (pipelineTraceGpu)
        {
            auto launchPrepare = [&](int start) {
                return std::async(std::launch::async, prepareBetaBlock, ib, start, gammaChunk);
            };

            int gammaStart = 0;
            std::future<BetaBlock> futureBlock = launchPrepare(gammaStart);
            gammaStart += gammaChunk;

            while (true)
            {
                BetaBlock block = futureBlock.get();
                const bool haveNext = gammaStart < gammaFullCount;
                if (haveNext)
                {
                    futureBlock = launchPrepare(gammaStart);
                    gammaStart += gammaChunk;
                }

                processBlock(block);

                if (!haveNext)
                    break;
            }
        }
        else
        {
            for (int gammaStart = 0; gammaStart < gammaFullCount; gammaStart += gammaChunk)
            {
                BetaBlock block = prepareBetaBlock(ib, gammaStart, gammaChunk);
                processBlock(block);
            }
        }

        if (m_mirrorGamma)
        {
            ApplyMirrorGammaMueller(betaM, nAz, nZen);
            if (computeNoShadow)
                ApplyMirrorGammaMueller(betaM_ns, nAz, nZen);
        }

        // Accumulate into global Mueller
        for (int p=0;p<nAz;++p) for (int t=0;t<=nZen;++t) {
            handlerPO->M.insert(p,t,betaM(p,t));
            if (computeNoShadow)
                handlerPO->M_noshadow.insert(p,t,betaM_ns(p,t));
        }

        // --save_betas: write per-beta Mueller (phi-averaged)
        if (m_saveBetas && m_mpiRank == 0)
        {
            auto &sphere = handlerPO->m_sphere;
            matrix *Lp = handlerPO->m_Lp;

            // Per-beta contribution
            std::string fname = betaDir + "/beta_" + std::to_string(ib)
                + "_" + std::to_string((int)RadToDeg(beta)) + "deg.dat";
            std::ofstream bf(fname, std::ios::out);
            bf << std::setprecision(10);
            bf << "# Per-beta Mueller: beta=" << RadToDeg(beta) << " deg, dcos=" << dcos
               << ", " << nGamma << " gamma orientations" << std::endl;
            bf << "ScAngle 2pi*dcos M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44";
            for (int iZen = 0; iZen <= nZen; ++iZen) {
                matrix Msum(4,4); Msum.Fill(0.0);
                double radZen = sphere.GetZenith(iZen);
                for (int iAz = 0; iAz < nAz; ++iAz) {
                    double radAz = -iAz * sphere.azinuthStep;
                    matrix m = betaM(iAz, iZen);
                    (*Lp)[1][1] = cos(2*radAz); (*Lp)[1][2] = sin(2*radAz);
                    (*Lp)[2][1] = -(*Lp)[1][2]; (*Lp)[2][2] = (*Lp)[1][1];
                    Msum += m * (*Lp);
                }
                Msum /= nAz;
                double _2PiDcos = sphere.Compute2PiDcos(iZen);
                bf << std::endl << RadToDeg(radZen) << ' ' << _2PiDcos << ' ';
                bf << Msum;
            }
            bf.close();

            // Also write cumulative
            std::string cfname = betaDir + "/cumul_" + std::to_string(ib)
                + "_" + std::to_string((int)RadToDeg(beta)) + "deg.dat";
            std::ofstream cf(cfname, std::ios::out);
            cf << std::setprecision(10);
            cf << "# Cumulative Mueller up to beta=" << RadToDeg(beta) << " deg"
               << " (" << (ib+1)*nGamma << "/" << nOrientations << " orientations)" << std::endl;
            cf << "ScAngle 2pi*dcos M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44";
            for (int iZen = 0; iZen <= nZen; ++iZen) {
                matrix Msum(4,4); Msum.Fill(0.0);
                double radZen = sphere.GetZenith(iZen);
                for (int iAz = 0; iAz < nAz; ++iAz) {
                    double radAz = -iAz * sphere.azinuthStep;
                    matrix m = handlerPO->M(iAz, iZen);
                    (*Lp)[1][1] = cos(2*radAz); (*Lp)[1][2] = sin(2*radAz);
                    (*Lp)[2][1] = -(*Lp)[1][2]; (*Lp)[2][2] = (*Lp)[1][1];
                    Msum += m * (*Lp);
                }
                Msum /= nAz;
                double _2PiDcos = sphere.Compute2PiDcos(iZen);
                cf << std::endl << RadToDeg(radZen) << ' ' << _2PiDcos << ' ';
                cf << Msum;
            }
            cf.close();
        }

        if (oldautoCheckpoint && m_mpiRank == 0)
        {
            SaveOldautoCheckpoint(oldautoCkptPath,
                                  handlerPO->M, handlerPO->M_noshadow,
                                  m_incomingEnergy, handlerPO->m_outputEnergy,
                                  handlerPO->m_extinctionCrossSectionOt,
                                  handlerPO->m_hasExtinctionOt,
                                  ib + 1, nBeta, nGamma,
                                  oldautoParamHash, nAz + 1, nZen + 1,
                                  handlerPO->ComputeNoShadow());
            std::ostringstream line;
            struct stat st;
            const bool haveStat = stat(oldautoCkptPath.c_str(), &st) == 0;
            line << "Oldauto checkpoint v2 saved: beta " << (ib + 1)
                 << "/" << nBeta
                 << ", orient " << ((long long)(ib + 1) * nGamma)
                 << "/" << nOrientations
                 << ", noShadow=" << (handlerPO->ComputeNoShadow() ? 1 : 0);
            if (haveStat)
                line << ", size=" << (st.st_size / (1024 * 1024)) << " MB";
            line << " -> " << oldautoCkptPath;
            std::cout << line.str() << std::endl;
            AppendTextLog(line.str() + "\n");
        }
    }

    MPI_ReduceMueller(handlerPO, nAz, nZen, m_incomingEnergy, m_mpiRank);

    EraseConsoleLine(60);
    if (m_mpiRank == 0) {
        std::ostringstream line;
        line << "Sequential tracing completed: " << count << "/" << nOrientations
             << " trace calls, phase1=" << std::fixed << std::setprecision(2)
             << phase1_total << " s";
        std::cout << line.str() << std::endl;
        AppendTextLog(line.str() + "\n");
        std::cout << "Phase 1 (tracing): " << std::fixed << std::setprecision(2) << phase1_total << " s" << std::endl;
        std::cout << "Phase 2 (diffraction, "
                  << (handlerPO->IsGpuEnabled() ? "CUDA" : "OpenMP")
                  << "): " << phase2_total << " s" << std::endl;
        std::cout << "Total: " << phase1_total + phase2_total << " s" << std::endl;
        if (m_saveBetas)
            std::cout << "Saved " << nBeta << " beta files to " << betaDir << "/" << std::endl;

        m_handler->WriteTotalMatricesToFile(m_resultDirName);
        m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
        if (handlerPO->ComputeNoShadow())
        {
            std::swap(handlerPO->M, handlerPO->M_noshadow);
            std::string nsName = m_resultDirName + "_noshadow";
            handlerPO->WriteMatricesToFile(nsName, m_incomingEnergy);
            std::swap(handlerPO->M, handlerPO->M_noshadow);
        }
        OutputStatisticsPO(timer, nOrientations, m_resultDirName);
        if (oldautoCheckpoint)
            std::remove(oldautoCkptPath.c_str());
    }
}


void TracerPOTotal::TraceMonteCarlo(const AngleRange &betaRange,
                                    const AngleRange &gammaRange,
                                    int nOrientations)
{
    // Generate random orientations, then use same chunked+OpenMP pipeline
    std::vector<std::pair<double,double>> orientations;
    std::vector<double> weights;

    unsigned int lo, hi;
    asm volatile("rdtsc" : "=a"(lo), "=d"(hi));
    srand(lo ^ hi);

    const double betaSpan = betaRange.max - betaRange.min;
    const double betaIntegral = cos(betaRange.min) - cos(betaRange.max);
    const double betaWeightNorm = (fabs(betaIntegral) > DBL_EPSILON)
        ? betaSpan / (betaIntegral * nOrientations)
        : 0.0;

    for (int i = 0; i < nOrientations; ++i)
    {
        double beta = betaRange.min + RandomDouble(0, 1) * betaSpan;
        double gamma = gammaRange.min + RandomDouble(0, 1) * (gammaRange.max - gammaRange.min);
        orientations.push_back({beta, gamma});
        weights.push_back(sin(beta) * betaWeightNorm);
    }

    std::cout << "Monte Carlo: " << nOrientations << " orientations" << std::endl;

    CalcTimer timer;
    timer.Start();
    if (m_mpiRank == 0) OutputStartTime(timer);

    m_handler->SetNormIndex(1);

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) { std::cerr << "Error: not HandlerPO" << std::endl; return; }
    const bool computeNoShadow = handlerPO->ComputeNoShadow();

    m_incomingEnergy = 0;
    handlerPO->m_outputEnergy = 0;
    handlerPO->m_extinctionCrossSectionOt = 0;
    handlerPO->m_hasExtinctionOt = false;
    handlerPO->M.ClearArr();
    handlerPO->M_noshadow.ClearArr();
    handlerPO->CleanJ();

    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    // MPI split
    int myStart = m_mpiRank * nOrientations / m_mpiSize;
    int myEnd = (m_mpiRank + 1) * nOrientations / m_mpiSize;
    int myCount = myEnd - myStart;

    // Chunked + OpenMP (same as TraceRandom)
    std::vector<Beam> outBeams;
    for (int i = myStart; i < myEnd; ++i)
    {
        m_particle->Rotate(orientations[i].first, orientations[i].second, 0);
        if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
        bool ok = m_scattering->ScatterLight(0, 0, outBeams);

        if (ok)
        {
            PreparedOrientation prepared;
            handlerPO->PrepareBeams(outBeams, weights[i], prepared);

            std::vector<Arr2DC> localJ, localJ_ns;
            if (handlerPO->isCoh) {
                Arr2DC t1(nAz+1,nZen+1,2,2); t1.ClearArr(); localJ.push_back(t1);
                Arr2DC t2(nAz+1,nZen+1,2,2); t2.ClearArr(); localJ_ns.push_back(t2);
            }
            handlerPO->HandleBeamsToLocal(prepared, handlerPO->M, localJ,
                                           handlerPO->isCoh ? &localJ_ns : nullptr);
            if (handlerPO->isCoh && !localJ.empty()) {
                HandlerPO::AddToMuellerLocal(localJ, prepared.sinZenith,
                                              handlerPO->M, nAz, nZen);
                if (computeNoShadow)
                    HandlerPO::AddToMuellerLocal(localJ_ns, prepared.sinZenith,
                                                  handlerPO->M_noshadow, nAz, nZen);
            }
        }
        m_incomingEnergy += m_scattering->GetIncedentEnergy() * weights[i];
        outBeams.clear();
    }

    MPI_ReduceMueller(handlerPO, nAz, nZen, m_incomingEnergy, m_mpiRank);

    if (m_mpiRank == 0) {
        m_handler->WriteTotalMatricesToFile(m_resultDirName);
        m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
        if (handlerPO->ComputeNoShadow())
        {
            std::swap(handlerPO->M, handlerPO->M_noshadow);
            std::string nsName = m_resultDirName + "_noshadow";
            handlerPO->WriteMatricesToFile(nsName, m_incomingEnergy);
            std::swap(handlerPO->M, handlerPO->M_noshadow);
        }
        OutputStatisticsPO(timer, nOrientations, m_resultDirName);
    }
}

void TracerPOTotal::TraceFromFile(const std::string &orientFile)
{
    // Read orientations from file
    std::ifstream inFile(orientFile);
    if (!inFile.is_open())
    {
        std::cerr << "Error! Cannot open orientation file: " << orientFile << std::endl;
        throw std::exception();
    }

    std::vector<std::pair<double,double>> orientations;
    std::string line;
    while (std::getline(inFile, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        double b, g;
        if (iss >> b >> g)
            orientations.push_back({b, g});
    }
    inFile.close();

    int nOrientations = orientations.size();
    if (nOrientations == 0)
    {
        std::cerr << "Error! No orientations in file: " << orientFile << std::endl;
        throw std::exception();
    }

    CalcTimer timer;
    timer.Start();
    OutputStartTime(timer);

    double weight = 1.0 / nOrientations;
    m_handler->SetNormIndex(1);

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO)
    {
        std::cerr << "Error! Handler is not HandlerPO in TraceFromFile" << std::endl;
        throw std::exception();
    }
    const bool computeNoShadow = handlerPO->ComputeNoShadow();

    m_incomingEnergy = 0;
    handlerPO->m_outputEnergy = 0;
    handlerPO->m_extinctionCrossSectionOt = 0;
    handlerPO->m_hasExtinctionOt = false;
    handlerPO->M.ClearArr();
    handlerPO->M_noshadow.ClearArr();
    handlerPO->CleanJ();

    // =========================================================================
    // Phase 1 (sequential): Trace beams for all orientations, preprocess them.
    //
    // WHY SEQUENTIAL: Particle::Rotate() and ScatterLight() modify shared
    // state (particle geometry, scattering internals). Parallelizing would
    // require N_THREADS copies of Particle and Scattering (no Clone() method
    // exists, and the virtual hierarchy makes it complex to implement).
    //
    // PERFORMANCE: Phase 1 is typically <5% of total time. The dominant
    // Phase 2 (diffraction integrals) IS parallelized with OpenMP.
    // =========================================================================
    // =========================================================================
    // Chunked streaming: process orientations in chunks to limit memory.
    // Each chunk: Phase 1 (sequential trace) → Phase 2 (parallel diffraction).
    // Memory: O(chunkSize * beams_per_orient * sizeof(PreparedBeam)).
    //
    // Chunk size: auto-select based on available RAM.
    // PreparedBeam ~ 3.5 KB, ~100 beams/orient → ~350 KB/orient.
    // Target: use at most 50% of available RAM for beam storage.
    // =========================================================================
    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    // Estimate available memory. In --multigrid_parallel the parent can pass a
    // per-child cap via MBS_HOST_MEM_BUDGET_MB to avoid N children overbooking RAM.
    long long availMemBytes = EffectiveMemAvailableMb() * 1024LL * 1024LL;
    // Reserve memory for thread-local Mueller arrays: nThreads * nAz * nZen * 16 * 8 bytes
    int nThreads = 1;
    #pragma omp parallel
    {
        #pragma omp single
        nThreads = omp_get_num_threads();
    }
    long long muellerMem = (long long)nThreads * (nAz+1) * (nZen+1) * 16 * 8 * 2; // Mueller + Jones
    long long beamBudget = (availMemBytes / 2) - muellerMem; // 50% of RAM for beams
    if (beamBudget < 100LL * 1024 * 1024) beamBudget = 100LL * 1024 * 1024; // min 100 MB

    long long bytesPerOrient = 350LL * 1024; // ~350 KB estimated (100 beams × 3.5 KB)
    int chunkSize = std::max(32, std::min(4096, std::min(nOrientations, (int)(beamBudget / bytesPerOrient))));
    // Round up to nice number for Sobol (power of 2 or at least multiple of nThreads)
    if (chunkSize >= nOrientations) chunkSize = nOrientations;

    int nChunks = (nOrientations + chunkSize - 1) / chunkSize;

    std::cerr << "Memory: " << availMemBytes / (1024*1024) << " MB available, "
              << "chunk=" << chunkSize << " orientations (" << nChunks << " chunks), "
              << nThreads << " threads" << std::endl;

    auto t_total_start = std::chrono::high_resolution_clock::now();
    double phase1_total = 0, phase2_total = 0;

    std::vector<Beam> outBeams;
    long long count = 0;

    // Checkpoint is opt-in because frequent binary dumps add avoidable I/O
    // overhead to ordinary orientation-file runs.
    std::string ckptPath = m_resultDirName + "_checkpoint.bin";
    const ::complex ri = m_particle->GetRefractiveIndex();
    unsigned long paramHash = HashParams(m_scattering->m_wave,
        real(ri), imag(ri), m_scattering->GetMaxReflections(),
        nOrientations, m_particle->MaximalDimention(), 0, nAz, nZen);
    MixHashInt(paramHash, computeNoShadow ? 1 : 0);
    for (const auto &bg : orientations)
    {
        paramHash ^= std::hash<double>{}(bg.first)
            + 0x9e3779b9 + (paramHash << 6) + (paramHash >> 2);
        paramHash ^= std::hash<double>{}(bg.second)
            + 0x9e3779b9 + (paramHash << 6) + (paramHash >> 2);
    }
    int resumeChunk = 0;
    if (m_enableCheckpoint)
    {
        int completedOrient = 0;
        if (LoadCheckpoint(ckptPath, handlerPO->M, handlerPO->M_noshadow,
                           m_incomingEnergy, handlerPO->m_outputEnergy,
                           completedOrient, nOrientations, paramHash,
                           computeNoShadow))
        {
            resumeChunk = completedOrient / chunkSize;
            count = completedOrient;
            if (m_mpiRank == 0)
                std::cout << "*** RESUMED from checkpoint: " << completedOrient
                          << "/" << nOrientations << " orientations (chunk " << resumeChunk
                          << "/" << nChunks << ")" << std::endl;
        }
    }

    for (int chunk = resumeChunk; chunk < nChunks; ++chunk)
    {
        int iStart = chunk * chunkSize;
        int iEnd = std::min(iStart + chunkSize, nOrientations);
        int thisChunkSize = iEnd - iStart;

        // Phase 1: trace and preprocess this chunk
        auto t_p1 = std::chrono::high_resolution_clock::now();

        std::vector<PreparedOrientation> chunkPrepared(thisChunkSize);

        for (int i = 0; i < thisChunkSize; ++i)
        {
            int idx = iStart + i;
            m_particle->Rotate(orientations[idx].first, orientations[idx].second, 0);

            if (!shadowOff)
                m_scattering->FormShadowBeam(outBeams);

            bool ok = m_scattering->ScatterLight(0, 0, outBeams);

            if (ok)
                handlerPO->PrepareBeams(outBeams, weight, chunkPrepared[i]);
            else
                chunkPrepared[i].sinZenith = weight;

            m_incomingEnergy += m_scattering->GetIncedentEnergy() * weight;
            if (m_mpiRank == 0) OutputProgress(nOrientations, count + 1, iStart + i, 0, timer, outBeams.size());
            outBeams.clear();
            ++count;
        }

        auto t_p1_end = std::chrono::high_resolution_clock::now();
        phase1_total += std::chrono::duration<double>(t_p1_end - t_p1).count();

        // Phase 2: parallel diffraction for this chunk
        auto t_p2 = std::chrono::high_resolution_clock::now();

        #pragma omp parallel
        {
            Arr2D localM(nAz + 1, nZen + 1, 4, 4);
            localM.ClearArr();
            Arr2D localM_ns(nAz + 1, nZen + 1, 4, 4); // no-shadow
            localM_ns.ClearArr();

            std::vector<Arr2DC> localJ, localJ_ns;
            if (handlerPO->isCoh)
            {
                Arr2DC tmp(nAz + 1, nZen + 1, 2, 2);
                tmp.ClearArr();
                localJ.push_back(tmp);
                Arr2DC tmp2(nAz + 1, nZen + 1, 2, 2);
                tmp2.ClearArr();
                localJ_ns.push_back(tmp2);
            }

            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < thisChunkSize; ++i)
            {
                if (!chunkPrepared[i].beams.empty())
                {
                    handlerPO->HandleBeamsToLocal(chunkPrepared[i], localM, localJ,
                                                   handlerPO->isCoh ? &localJ_ns : nullptr);
                }

                if (handlerPO->isCoh && !localJ.empty())
                {
                    double w = chunkPrepared[i].sinZenith;
                    HandlerPO::AddToMuellerLocal(localJ, w, localM, nAz, nZen);
                    if (computeNoShadow)
                        HandlerPO::AddToMuellerLocal(
                            localJ_ns, w, localM_ns, nAz, nZen);
                    localJ[0].ClearArr();
                    localJ_ns[0].ClearArr();
                }
            }

            #pragma omp critical
            {
                for (int p = 0; p < nAz; ++p)
                    for (int t = 0; t <= nZen; ++t)
                    {
                        handlerPO->M.insert(p, t, localM(p, t));
                        if (computeNoShadow)
                            handlerPO->M_noshadow.insert(
                                p, t, localM_ns(p, t));
                    }
            }
        } // end omp parallel

        auto t_p2_end = std::chrono::high_resolution_clock::now();
        phase2_total += std::chrono::duration<double>(t_p2_end - t_p2).count();

        // Free chunk memory immediately
        chunkPrepared.clear();
        chunkPrepared.shrink_to_fit();

        // Save checkpoint after each chunk.
        if (m_enableCheckpoint && m_mpiRank == 0) {
            SaveCheckpoint(ckptPath, handlerPO->M, handlerPO->M_noshadow,
                           m_incomingEnergy, handlerPO->m_outputEnergy,
                           iEnd, nOrientations, paramHash, nAz+1, nZen+1,
                           computeNoShadow);
        }
    }

    EraseConsoleLine(60);
    if (m_mpiRank == 0)
    {
        std::ostringstream line;
        line << "Sequential tracing completed: " << count << "/" << nOrientations
             << " trace calls, phase1=" << std::fixed << std::setprecision(2)
             << phase1_total << " s";
        std::cout << line.str() << std::endl;
        AppendTextLog(line.str() + "\n");
    }
    std::cout << "Phase 1 (tracing + preprocessing): " << std::fixed
              << std::setprecision(2) << phase1_total << " s" << std::endl;

    // (t_phase2_end - t_phase2_start) compatibility: set from totals
    auto t_phase2_start = t_total_start; // dummy for downstream code
    auto t_phase2_end = std::chrono::high_resolution_clock::now();
    double phase2_sec = std::chrono::duration<double>(t_phase2_end - t_phase2_start).count();

    std::cout << "Phase 2 (diffraction, OpenMP): " << std::fixed
              << std::setprecision(2) << phase2_total << " s" << std::endl;
    std::cout << "Total: " << phase1_total + phase2_total << " s" << std::endl;

    m_handler->WriteTotalMatricesToFile(m_resultDirName);
    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);

    // Write no-shadow Mueller matrix (swap M <-> M_noshadow temporarily)
    {
        HandlerPO *hp = dynamic_cast<HandlerPO*>(m_handler);
        if (hp && hp->ComputeNoShadow()) {
            std::swap(hp->M, hp->M_noshadow);
            std::string nsName = m_resultDirName + "_noshadow";
            hp->WriteMatricesToFile(nsName, m_incomingEnergy);
            std::swap(hp->M, hp->M_noshadow); // swap back
        }
    }

    OutputStatisticsPO(timer, nOrientations, m_resultDirName);

    // Remove checkpoint on successful completion.
    if (m_enableCheckpoint)
        std::remove(ckptPath.c_str());
}

void TracerPOTotal::TraceFromFileMultiSize(const std::string &orientFile,
                                            const std::vector<double> &x_sizes,
                                            double x_ref)
{
    // Phase 1: Read orientations
    std::ifstream inFile(orientFile);
    if (!inFile.is_open())
    {
        std::cerr << "Error! Cannot open orientation file: " << orientFile << std::endl;
        throw std::exception();
    }

    std::vector<std::pair<double,double>> orientations;
    std::string line;
    while (std::getline(inFile, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        double b, g;
        if (iss >> b >> g)
            orientations.push_back({b, g});
    }
    inFile.close();

    int nOrientations = orientations.size();
    if (nOrientations == 0)
    {
        std::cerr << "Error! No orientations in file: " << orientFile << std::endl;
        throw std::exception();
    }

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO)
    {
        std::cerr << "Error! Handler is not HandlerPO in TraceFromFileMultiSize" << std::endl;
        throw std::exception();
    }
    // D_ref = maximal dimension of the particle as traced.
    // The particle is already set to x_ref size by the caller.
    double D_ref = m_particle->MaximalDimention();

    CalcTimer timer;
    long long count = 0;
    std::vector<Beam> outBeams;
    timer.Start();
    OutputStartTime(timer);

    double weight = 1.0 / nOrientations;

    // Phase 1: Trace and cache
    BeamCache cache;
    cache.D_ref = D_ref;
    cache.orientations.resize(nOrientations);

    auto t_cache_start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < nOrientations; ++i)
    {
        double beta  = orientations[i].first;
        double gamma = orientations[i].second;

        m_particle->Rotate(beta, gamma, 0);

        if (!shadowOff)
        {
            m_scattering->FormShadowBeam(outBeams);
        }

        bool ok = m_scattering->ScatterLight(0, 0, outBeams);

        double incomingE = m_scattering->GetIncedentEnergy() * weight;

        if (ok)
        {
            handlerPO->CacheBeams(outBeams, weight, D_ref, incomingE,
                                   cache.orientations[i]);
        }
        else
        {
            cache.orientations[i].weight = weight;
            cache.orientations[i].incomingEnergy = incomingE;
            if (m_mpiRank == 0) std::cout << std::endl << "Orientation " << i
                      << " (beta=" << beta << ", gamma=" << gamma
                      << ") has been skipped!!!" << std::endl;
        }

        OutputProgress(nOrientations, count + 1, i, 0, timer, outBeams.size());
        outBeams.clear();
        ++count;
    }

    auto t_cache_end = std::chrono::high_resolution_clock::now();
    double cache_sec = std::chrono::duration<double>(t_cache_end - t_cache_start).count();

    EraseConsoleLine(60);
    if (m_mpiRank == 0) std::cout << "Phase 1 (tracing): 100%" << std::endl;
    if (m_mpiRank == 0) std::cout << "Cached " << cache.totalBeams() << " beams from "
              << nOrientations << " orientations in " << cache_sec << " s" << std::endl;

    // Phase 2: Compute diffraction for all sizes
    auto t_diff_start = std::chrono::high_resolution_clock::now();

    std::vector<Arr2D> results_M;
    std::vector<double> results_energy;
    handlerPO->ComputeFromCache(cache, x_sizes, results_M, results_energy);

    auto t_diff_end = std::chrono::high_resolution_clock::now();
    double diff_sec = std::chrono::duration<double>(t_diff_end - t_diff_start).count();

    std::cout << "Phase 2 (diffraction for " << x_sizes.size() << " sizes): "
              << diff_sec << " s" << std::endl;
    std::cout << "Total: " << cache_sec + diff_sec << " s" << std::endl;

    // Write results for each size
    std::string baseName = m_resultDirName;

    // Save original M, write each size's result
    Arr2D &origM = handlerPO->M;

    for (size_t s = 0; s < x_sizes.size(); ++s)
    {
        // Swap in this size's Mueller matrix
        Arr2D savedM = origM;
        origM = results_M[s];
        m_incomingEnergy = results_energy[s];

        // Write with size suffix
        std::string sizeName = baseName + "_x" + SizeFileLabel(x_sizes[s]);
        handlerPO->WriteMatricesToFile(sizeName, m_incomingEnergy);

        // Also write in HandlerPOTotal format
        // Write the azimuth-averaged format
        {
            std::ofstream outFile(sizeName + ".dat", std::ios::out);
            outFile << std::setprecision(10);
            outFile << "ScAngle 2pi*dcos "
                    "M11 M12 M13 M14 "
                    "M21 M22 M23 M24 "
                    "M31 M32 M33 M34 "
                    "M41 M42 M43 M44";

            auto &sphere = handlerPO->m_sphere;
            matrix *Lp = handlerPO->m_Lp;
            int nZen = sphere.nZenith;
            int nAz = sphere.nAzimuth;

            for (int iZen = 0; iZen <= nZen; ++iZen)
            {
                matrix Msum(4, 4);
                Msum.Fill(0.0);
                double radZen = sphere.GetZenith(iZen);

                for (int iAz = 0; iAz < nAz; ++iAz)
                {
                    double radAz = -iAz * sphere.azinuthStep;
                    matrix m = results_M[s](iAz, iZen);

                    (*Lp)[1][1] = cos(2*radAz);
                    (*Lp)[1][2] = sin(2*radAz);
                    (*Lp)[2][1] = -(*Lp)[1][2];
                    (*Lp)[2][2] = (*Lp)[1][1];

                    matrix Ln = *Lp;
                    Ln[1][2] *= -1;
                    Ln[2][1] *= -1;

                    if (radZen > M_PI - __FLT_EPSILON__)
                        Msum += (*Lp) * m * (*Lp);
                    else if (radZen < __FLT_EPSILON__)
                        Msum += Ln * m * (*Lp);
                    else
                        Msum += m * (*Lp);
                }

                double _2Pi_dcos = sphere.Compute2PiDcos(iZen);

                Msum /= nAz;
                outFile << std::endl << RadToDeg(radZen) << ' ' << _2Pi_dcos << ' ';
                outFile << Msum;
            }
            outFile.close();
        }

        // Restore
        origM = savedM;
    }

    std::cout << "Results written to " << baseName << "_x*.dat" << std::endl;
}

void TracerPOTotal::TraceRandomMultiSize(const AngleRange &betaRange,
                                          const AngleRange &gammaRange,
                                          const std::vector<double> &x_sizes,
                                          const std::vector<std::string> &labels)
{
    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO)
    {
        std::cerr << "Error! Handler is not HandlerPO in TraceRandomMultiSize" << std::endl;
        throw std::exception();
    }
    if (x_sizes.empty())
        return;
    if (x_sizes.size() > 32)
    {
        std::cerr << "ERROR: shared multikeq currently supports at most 32 sizes per run "
                  << "(got " << x_sizes.size() << "). Split the range." << std::endl;
        throw std::runtime_error("too many multikeq sizes");
    }
    if (m_mpiSize > 1)
    {
        std::cerr << "ERROR: shared multikeq cache mode is single-process only for now." << std::endl;
        throw std::runtime_error("shared multikeq does not support MPI");
    }

    int nGamma = gammaRange.number;
    const bool betaMidpoint = OldautoBetaMidpointEnabled() && !m_fastPoleGamma;
    int nBeta = OldautoBetaCount(betaRange, betaMidpoint);
    int nOrientations = nBeta * nGamma;
    int betaNorm = (m_symmetry.beta < M_PI_2 + FLT_EPSILON
                    && m_symmetry.beta > M_PI_2 - FLT_EPSILON) ? 1 : 2;
    double normGamma = gammaRange.number * betaNorm;
    const bool gammaStagger = OldautoGammaStaggerEnabled();
    double D_ref = m_particle->MaximalDimention();
    double x_ref = M_PI * D_ref / m_scattering->m_wave;
    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    std::cout << "Shared multikeq: trace once on " << nBeta << " x " << nGamma
              << " = " << nOrientations << " orientations, then diffract "
              << x_sizes.size() << " sizes" << std::endl;
    std::cout << "  D_ref=" << D_ref << ", x_ref=" << x_ref
              << ", backend=" << (handlerPO->IsGpuEnabled()
                    ? (handlerPO->IsFftEnabled() ? "CUDA FFT" : "CUDA")
                    : "OpenMP") << std::endl;

    CalcTimer timer;
    timer.Start();
    OutputStartTime(timer);

    int nThreads = 1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        nThreads = omp_get_num_threads();
    }
#endif
    std::cout << "Shared multikeq Phase 1: trace/prepare with "
              << nThreads << " OpenMP threads"
              << (gammaStagger ? "; gamma stagger enabled" : "")
              << (betaMidpoint ? "; beta midpoint enabled" : "")
              << (!betaMidpoint && m_fastPoleGamma ? "; beta endpoints enabled for --pole" : "")
              << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    long long count = 0;
    double phase1 = 0.0;
    double phase2 = 0.0;
    const bool computeNoShadow = handlerPO->ComputeNoShadow();
    std::vector<Arr2D> results_M;
    std::vector<Arr2D> results_M_ns;
    std::vector<double> results_energy(x_sizes.size(), 0.0);
    std::vector<double> results_output_energy(x_sizes.size(), 0.0);
    std::vector<double> results_ext_ot(x_sizes.size(), 0.0);
    for (size_t s = 0; s < x_sizes.size(); ++s)
    {
        results_M.push_back(Arr2D(nAz + 1, nZen + 1, 4, 4));
        results_M.back().ClearArr();
        if (computeNoShadow)
        {
            results_M_ns.push_back(Arr2D(nAz + 1, nZen + 1, 4, 4));
            results_M_ns.back().ClearArr();
        }
    }

    int sharedBetaGroup = 1;
    if (const char *env = std::getenv("MBS_SHARED_BETA_GROUP"))
    {
        char *end = nullptr;
        long value = std::strtol(env, &end, 10);
        if (end && *end == '\0' && value > 0)
            sharedBetaGroup = (int)std::min<long>(value, nBeta);
    }
    else if (nThreads > 1 && nBeta > 1)
    {
        // A single beta block has only nGamma independent trace tasks. For
        // non-convex particles dynamic scheduling still leaves a long serial
        // tail when a few hard orientations split many beams. Group a few beta
        // blocks by default so OpenMP has enough independent orientations while
        // keeping prepared-beam memory bounded.
        int autoGroup = std::max(1, nThreads / 4);
        autoGroup = std::min(autoGroup, 4);
        autoGroup = std::min(autoGroup, nBeta);

        const int maxGroupOrient = std::max(nGamma, nThreads * 96);
        while (autoGroup > 1 && autoGroup * nGamma > maxGroupOrient)
            --autoGroup;

        sharedBetaGroup = autoGroup;
    }
    if (m_mpiRank == 0)
    {
        std::cout << "Shared multikeq beta grouping: " << sharedBetaGroup
                  << " beta block"
                  << (sharedBetaGroup == 1 ? "" : "s")
                  << " per diffraction pass";
        if (!std::getenv("MBS_SHARED_BETA_GROUP"))
            std::cout << " (auto; set MBS_SHARED_BETA_GROUP to override)";
        std::cout << std::endl;
    }

    int sharedOrientChunk = sharedBetaGroup * nGamma;
    bool orientChunkFromEnv = false;
    if (const char *env = std::getenv("MBS_SHARED_ORIENT_CHUNK"))
    {
        char *end = nullptr;
        long value = std::strtol(env, &end, 10);
        if (end && *end == '\0' && value > 0)
        {
            sharedOrientChunk = (int)std::min<long>(value, nOrientations);
            orientChunkFromEnv = true;
        }
    }
    else if (handlerPO->IsGpuEnabled())
    {
        // Keep GPU diffraction passes large enough to amortize launch,
        // transfer, and FFT interpolation overhead. The GPU backend still
        // splits this chunk internally by available VRAM.
        int autoChunk = handlerPO->ComputeNoShadow() ? 1024 : 4096;
        sharedOrientChunk = std::min(autoChunk, nOrientations);
    }
    sharedOrientChunk = std::max(1, std::min(sharedOrientChunk, nOrientations));
    if (m_mpiRank == 0)
    {
        std::cout << "Shared multikeq global scheduler: chunks of "
                  << sharedOrientChunk << " orientations";
        if (!orientChunkFromEnv)
            std::cout << " (auto; set MBS_SHARED_ORIENT_CHUNK to override)";
        std::cout << std::endl;
    }
    int debugOrientBegin = 0;
    int debugOrientEnd = nOrientations;
    if (const char *env = std::getenv("MBS_DEBUG_ORIENT_BEGIN"))
    {
        char *end = nullptr;
        long value = std::strtol(env, &end, 10);
        if (end && *end == '\0' && value >= 0 && value < nOrientations)
            debugOrientBegin = (int)value;
    }
    if (const char *env = std::getenv("MBS_DEBUG_ORIENT_END"))
    {
        char *end = nullptr;
        long value = std::strtol(env, &end, 10);
        if (end && *end == '\0' && value > debugOrientBegin
            && value <= nOrientations)
            debugOrientEnd = (int)value;
    }
    if (m_mpiRank == 0
        && (debugOrientBegin != 0 || debugOrientEnd != nOrientations))
    {
        std::cout << "DEBUG: limiting shared multikeq orientations to ["
                  << debugOrientBegin << ", " << debugOrientEnd
                  << ") via MBS_DEBUG_ORIENT_BEGIN/END" << std::endl;
    }

    bool sharedPipeline = false;
    if (const char *env = std::getenv("MBS_SHARED_PIPELINE"))
        sharedPipeline = (std::atoi(env) != 0);
    if (m_mpiRank == 0 && sharedPipeline)
        std::cout << "Shared multikeq pipeline: CPU tracing of the next beta group "
                  << "overlaps current GPU diffraction" << std::endl;
    const bool sharedGpuExplicitScale =
        handlerPO->IsGpuEnabled()
        && std::getenv("MBS_SHARED_GPU_EXPLICIT_SCALE")
        && std::atoi(std::getenv("MBS_SHARED_GPU_EXPLICIT_SCALE")) != 0;
    if (m_mpiRank == 0 && sharedGpuExplicitScale)
        std::cout << "Shared multikeq GPU: explicit CPU scaling enabled "
                  << "(MBS_SHARED_GPU_EXPLICIT_SCALE=1)" << std::endl;
    const bool sharedGpuMultiKFull =
        handlerPO->IsGpuEnabled()
        && !computeNoShadow
        && !sharedGpuExplicitScale
        && std::getenv("MBS_GPU_MULTI_K_FULL")
        && std::atoi(std::getenv("MBS_GPU_MULTI_K_FULL")) != 0;
    if (m_mpiRank == 0 && sharedGpuMultiKFull)
        std::cout << "Shared multikeq GPU: experimental fused multi-k kernel "
                  << "enabled" << std::endl;
    HandlerPO prepareTemplate(m_particle, &m_incidentLight,
                              handlerPO->nTheta, m_scattering->m_wave);
    prepareTemplate.ConfigureForThreadLocalPrepare(*handlerPO, m_scattering);
    m_scattering->PrepareForParallelTrace();

    struct SharedGroupData
    {
        int orientStart = 0;
        int orientEnd = 0;
        int groupOrient = 0;
        std::vector<PreparedOrientation> prepared;
        std::vector<double> energies;
        std::vector<double> outputEnergies;
        std::vector<int> beamCounts;
        double traceSeconds = 0.0;
    };

    auto traceGroup = [&](int orientStart) -> SharedGroupData
    {
        SharedGroupData group;
        group.orientStart = orientStart;
        group.orientEnd = std::min(debugOrientEnd, orientStart + sharedOrientChunk);
        group.groupOrient = group.orientEnd - group.orientStart;
        group.prepared.resize(group.groupOrient);
        group.energies.assign(group.groupOrient, 0.0);
        group.outputEnergies.assign(group.groupOrient, 0.0);
        group.beamCounts.assign(group.groupOrient, 0);

        auto groupTraceStart = std::chrono::high_resolution_clock::now();
        ParallelExceptionState parallelError;

        #pragma omp parallel
        {
            Particle localParticle = *m_particle;
            Scattering *localScatter = m_scattering->CloneFor(&localParticle, &m_incidentLight);
            HandlerPO localHandler(&localParticle, &m_incidentLight,
                                   handlerPO->nTheta, m_scattering->m_wave);
            localHandler.ConfigureForThreadLocalPrepare(prepareTemplate, localScatter);
            std::vector<Beam> localBeams;

            #pragma omp for schedule(dynamic, 1)
            for (int localIndex = 0; localIndex < group.groupOrient; ++localIndex)
            {
                if (parallelError.Failed())
                    continue;
                try
                {
                const int globalIndex = group.orientStart + localIndex;
                const int ib = globalIndex / nGamma;
                const int jg = globalIndex - ib * nGamma;
                const int gi = localIndex;

                double beta = OldautoBetaAngle(betaRange, ib, betaMidpoint);
                double dcos = 0.0;
                CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normGamma, dcos);
                const bool pole = (fabs(beta) <= FLT_EPSILON
                                   || fabs(beta - M_PI) <= FLT_EPSILON);
                double traceBeta = beta;
                if (pole && !m_fastPoleGamma && betaRange.step > 0.0)
                {
                    if (fabs(beta) <= FLT_EPSILON)
                        traceBeta = beta + 0.5 * betaRange.step;
                    else
                        traceBeta = beta - 0.5 * betaRange.step;
                }

                double gamma = OldautoGammaAngle(gammaRange, nGamma, jg, ib,
                                                 gammaStagger);
                localParticle.Rotate(traceBeta, gamma, 0);
                if (!shadowOff)
                    localScatter->FormShadowBeam(localBeams);
                bool ok = localScatter->ScatterLight(0, 0, localBeams);
                group.beamCounts[gi] = (int)localBeams.size();
                group.energies[gi] = localScatter->GetIncedentEnergy() * dcos;
                if (ok)
                {
                    double beforeOutput = localHandler.m_outputEnergy;
                    localHandler.PrepareBeams(localBeams, dcos, group.prepared[gi]);
                    group.outputEnergies[gi] = localHandler.m_outputEnergy - beforeOutput;
                }
                else
                {
                    group.prepared[gi].sinZenith = dcos;
                }
                localBeams.clear();
                }
                catch (...)
                {
                    parallelError.Capture();
                }
            }

            delete localScatter;
        }
        parallelError.Rethrow();
        group.traceSeconds = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - groupTraceStart).count();
        return group;
    };

    auto processGroup = [&](const SharedGroupData &group)
    {
        phase1 += group.traceSeconds;
        for (int i = 0; i < group.groupOrient; ++i)
        {
            int idx = group.orientStart + i;
            ++count;
            OutputProgress(nOrientations, count, idx, 0, timer,
                           group.beamCounts[i]);
        }

        auto betaDiffStart = std::chrono::high_resolution_clock::now();
        if (sharedGpuMultiKFull)
        {
            std::vector<double> scales(x_sizes.size(), 1.0);
            std::vector<double> localExtOt(x_sizes.size(), 0.0);
            for (size_t s = 0; s < x_sizes.size(); ++s)
            {
                const double scale = x_sizes[s] / x_ref;
                scales[s] = scale;
                const double waveIndex = 2.0 * M_PI / m_scattering->m_wave;
                for (const PreparedOrientation &po : group.prepared)
                {
                    localExtOt[s] +=
                        handlerPO->ComputeForwardExtinctionOtScaled(
                            po, scale, waveIndex,
                            handlerPO->AbsorptionCoefficient());
                }
            }

            std::vector<Arr2D> localMs;
            localMs.reserve(x_sizes.size());
            for (size_t s = 0; s < x_sizes.size(); ++s)
            {
                localMs.push_back(Arr2D(nAz + 1, nZen + 1, 4, 4));
                localMs.back().ClearArr();
            }

            bool fusedOk = true;
            for (int gpuStart = 0; gpuStart < group.groupOrient; )
            {
                int gpuBatchSize = handlerPO->SelectGpuOrientationBatchSize(
                    group.prepared, gpuStart, group.groupOrient - gpuStart);
                int gpuEnd = std::min(gpuStart + gpuBatchSize, group.groupOrient);
                fusedOk = handlerPO->IsFftEnabled()
                    ? handlerPO->HandleOrientationsToLocalGpuFftPhiMultiK(
                        group.prepared, gpuStart, gpuEnd - gpuStart, scales,
                        2.0 * M_PI / m_scattering->m_wave, localMs)
                    : handlerPO->HandleOrientationsToLocalGpuMultiK(
                        group.prepared, gpuStart, gpuEnd - gpuStart, scales,
                        2.0 * M_PI / m_scattering->m_wave, localMs);
                if (!fusedOk)
                    break;
                gpuStart = gpuEnd;
            }

            if (fusedOk)
            {
                for (size_t s = 0; s < x_sizes.size(); ++s)
                {
                    const double scale = scales[s];
                    const double scale2 = scale * scale;
                    results_ext_ot[s] += localExtOt[s];
                    if (m_mirrorGamma)
                        ApplyMirrorGammaMueller(localMs[s], nAz, nZen);
                    for (int p=0; p<nAz; ++p)
                        for (int t=0; t<=nZen; ++t)
                            results_M[s].insert(p, t, localMs[s](p, t));
                    for (int i = 0; i < group.groupOrient; ++i)
                    {
                        results_energy[s] += group.energies[i] * scale2;
                        results_output_energy[s] += PreparedOutputEnergy(
                            group.prepared[i], scale,
                            handlerPO->AbsorptionCoefficient());
                    }
                }
                phase2 += std::chrono::duration<double>(
                    std::chrono::high_resolution_clock::now() - betaDiffStart).count();
                return;
            }
            if (m_mpiRank == 0)
                std::cerr << "Shared multikeq GPU: fused multi-k unavailable "
                          << "for this chunk; falling back to per-size GPU path."
                          << std::endl;
        }
        for (size_t s = 0; s < x_sizes.size(); ++s)
        {
            double scale = x_sizes[s] / x_ref;
            double scale2 = scale * scale;
            const double waveIndex = 2.0 * M_PI / m_scattering->m_wave;

            for (const PreparedOrientation &po : group.prepared)
            {
                results_ext_ot[s] +=
                    handlerPO->ComputeForwardExtinctionOtScaled(
                        po, scale, waveIndex, handlerPO->AbsorptionCoefficient());
            }

            Arr2D localM(nAz + 1, nZen + 1, 4, 4); localM.ClearArr();
            Arr2D localM_ns(nAz + 1, nZen + 1, 4, 4); localM_ns.ClearArr();

            if (handlerPO->IsGpuEnabled())
            {
                const std::vector<PreparedOrientation> *gpuPrepared = &group.prepared;
                std::vector<PreparedOrientation> scaled;
                double gpuScale = scale;
                double gpuWaveIndex = waveIndex;
                if (sharedGpuExplicitScale)
                {
                    scaled.reserve(group.prepared.size());
                    for (const PreparedOrientation &po : group.prepared)
                        scaled.push_back(ScalePreparedOrientation(
                            po, scale, waveIndex, handlerPO->AbsorptionCoefficient()));
                    gpuPrepared = &scaled;
                    gpuScale = 1.0;
                    gpuWaveIndex = 0.0;
                }

                for (int gpuStart = 0; gpuStart < group.groupOrient; )
                {
                    int gpuBatchSize = handlerPO->SelectGpuOrientationBatchSize(
                        *gpuPrepared, gpuStart, group.groupOrient - gpuStart);
                    int gpuEnd = std::min(gpuStart + gpuBatchSize, group.groupOrient);
                    bool ok = handlerPO->IsFftEnabled()
                        ? handlerPO->HandleOrientationsToLocalGpuFftPhi(
                            *gpuPrepared, gpuStart, gpuEnd - gpuStart, localM, localM_ns,
                            gpuScale, gpuWaveIndex)
                        : handlerPO->HandleOrientationsToLocalGpu(
                            *gpuPrepared, gpuStart, gpuEnd - gpuStart, localM, localM_ns,
                            gpuScale, gpuWaveIndex);
                    if (!ok)
                    {
                        std::cerr << "ERROR: shared multikeq GPU backend failed." << std::endl;
                        throw std::runtime_error("shared multikeq GPU backend failed");
                    }
                    gpuStart = gpuEnd;
                }
            }
            else
            {
                std::vector<PreparedOrientation> scaled;
                scaled.reserve(group.prepared.size());
                for (const PreparedOrientation &po : group.prepared)
                    scaled.push_back(ScalePreparedOrientation(
                        po, scale, 2.0 * M_PI / m_scattering->m_wave,
                        handlerPO->AbsorptionCoefficient()));

                #pragma omp parallel
                {
                    Arr2D threadM(nAz + 1, nZen + 1, 4, 4); threadM.ClearArr();
                    Arr2D threadMns(nAz + 1, nZen + 1, 4, 4); threadMns.ClearArr();
                    std::vector<Arr2DC> localJ, localJns;
                    if (handlerPO->isCoh) {
                        Arr2DC t1(nAz+1,nZen+1,2,2); t1.ClearArr(); localJ.push_back(t1);
                        Arr2DC t2(nAz+1,nZen+1,2,2); t2.ClearArr(); localJns.push_back(t2);
                    }
                    #pragma omp for schedule(dynamic, 1)
                    for (int i = 0; i < group.groupOrient; ++i) {
                        if (!scaled[i].beams.empty())
                            handlerPO->HandleBeamsToLocal(scaled[i], threadM, localJ,
                                                           (handlerPO->isCoh && computeNoShadow) ? &localJns : nullptr);
                        if (handlerPO->isCoh && !localJ.empty()) {
                            double w = scaled[i].sinZenith;
                            HandlerPO::AddToMuellerLocal(localJ, w, threadM, nAz, nZen);
                            if (computeNoShadow)
                                HandlerPO::AddToMuellerLocal(localJns, w, threadMns, nAz, nZen);
                            localJ[0].ClearArr();
                            if (computeNoShadow)
                                localJns[0].ClearArr();
                        }
                    }
                    #pragma omp critical
                    {
                        for (int p=0; p<nAz; ++p) for (int t=0; t<=nZen; ++t) {
                            localM.insert(p, t, threadM(p, t));
                            if (computeNoShadow)
                                localM_ns.insert(p, t, threadMns(p, t));
                        }
                    }
                }
            }

            if (m_mirrorGamma)
            {
                ApplyMirrorGammaMueller(localM, nAz, nZen);
                if (computeNoShadow)
                    ApplyMirrorGammaMueller(localM_ns, nAz, nZen);
            }

            for (int p=0; p<nAz; ++p) for (int t=0; t<=nZen; ++t) {
                results_M[s].insert(p, t, localM(p, t));
                if (computeNoShadow)
                    results_M_ns[s].insert(p, t, localM_ns(p, t));
            }
            for (int i = 0; i < group.groupOrient; ++i) {
                results_energy[s] += group.energies[i] * scale2;
                results_output_energy[s] += PreparedOutputEnergy(
                    group.prepared[i], scale, handlerPO->AbsorptionCoefficient());
            }
        }
        phase2 += std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - betaDiffStart).count();
    };

    if (sharedPipeline && nOrientations > sharedOrientChunk)
    {
        int orientStart = debugOrientBegin;
        std::future<SharedGroupData> current =
            std::async(std::launch::async, traceGroup, orientStart);
        orientStart += sharedOrientChunk;
        while (orientStart < debugOrientEnd)
        {
            SharedGroupData group = current.get();
            current = std::async(std::launch::async, traceGroup, orientStart);
            orientStart += sharedOrientChunk;
            processGroup(group);
        }
        processGroup(current.get());
    }
    else
    {
        for (int orientStart = debugOrientBegin; orientStart < debugOrientEnd;
             orientStart += sharedOrientChunk)
            processGroup(traceGroup(orientStart));
    }

    double refineSeconds = 0.0;
    if (handlerPO->IsGpuEnabled() && handlerPO->IsFftEnabled()
        && SharedFftGlobalRefineEnabled())
    {
        const double refineThreshold = SharedFftGlobalRefineThreshold();
        const int refineMaxRows = SharedFftGlobalRefineMaxRows();
        const bool spikeDebug = SharedFftSpikeDebugEnabled();
        std::ofstream spikeDebugFile;
        if (spikeDebug && m_mpiRank == 0)
        {
            spikeDebugFile.open(SharedFftSpikeDebugPath().c_str(), std::ios::out);
            if (spikeDebugFile.is_open())
            {
                spikeDebugFile << std::setprecision(17);
                WriteSpikeDebugHeader(spikeDebugFile);
            }
            else
                std::cerr << "WARNING: cannot open MBS_FFT_SPIKE_DEBUG_FILE="
                          << SharedFftSpikeDebugPath() << std::endl;
        }
        ScatteringRange savedSphere = handlerPO->m_sphere;
        std::vector<std::vector<int>> refineRows(x_sizes.size());
        for (size_t s = 0; s < x_sizes.size(); ++s)
        {
            refineRows[s] = DetectGlobalFftSpikeRows(results_M[s], nAz, nZen,
                                                     refineThreshold,
                                                     refineMaxRows);
            if (!refineRows[s].empty() && m_mpiRank == 0)
            {
                std::cout << "GPU FFT global phi refine for "
                          << ((s < labels.size() && !labels[s].empty())
                                  ? labels[s] : std::to_string(s))
                          << ": " << refineRows[s].size()
                          << " theta rows at full N_phi=" << nAz
                          << " (threshold=" << refineThreshold << "):";
                for (int t : refineRows[s])
                    std::cout << ' ' << RadToDeg(savedSphere.GetZenith(t));
                std::cout << std::endl;
            }
        }

        auto refineStart = std::chrono::high_resolution_clock::now();
        std::vector<ScatteringRange> rowSpheres;
        rowSpheres.reserve(x_sizes.size());
        std::vector<int> nRefineRows(x_sizes.size(), 0);
        std::vector<Arr2D> refinedM(x_sizes.size());
        std::vector<Arr2D> refinedMns(x_sizes.size());
        bool anyRefineRows = false;
        for (size_t s = 0; s < x_sizes.size(); ++s)
        {
            rowSpheres.push_back(ScatteringRange(0.0, 0.0, nAz, 1));
            const std::vector<int> &rows = refineRows[s];
            if (rows.empty())
                continue;

            anyRefineRows = true;
            nRefineRows[s] = (int)rows.size();
            ScatteringRange &rowSphere = rowSpheres[s];
            rowSphere.isNonUniform = true;
            rowSphere.thetaValues.clear();
            rowSphere.thetaValues.reserve(rows.size());
            for (int row : rows)
                rowSphere.thetaValues.push_back(savedSphere.GetZenith(row));
            rowSphere.nAzimuth = nAz;
            rowSphere.azinuthStep = M_2PI / nAz;
            rowSphere.nZenith = (int)rowSphere.thetaValues.size() - 1;
            rowSphere.zenithStart = rowSphere.thetaValues.front();
            rowSphere.zenithEnd = rowSphere.thetaValues.back();
            rowSphere.zenithStep = rowSphere.nZenith > 0
                ? (rowSphere.zenithEnd - rowSphere.zenithStart) / rowSphere.nZenith
                : 0.0;
            rowSphere.ComputeSphereDirections(m_incidentLight);

            refinedM[s] = Arr2D(nAz + 1, nRefineRows[s], 4, 4);
            refinedM[s].ClearArr();
            if (computeNoShadow)
            {
                refinedMns[s] = Arr2D(nAz + 1, nRefineRows[s], 4, 4);
                refinedMns[s].ClearArr();
            }
        }

        const double waveIndex = 2.0 * M_PI / m_scattering->m_wave;
        if (anyRefineRows)
        {
            for (int orientStart = debugOrientBegin; orientStart < debugOrientEnd;
                 orientStart += sharedOrientChunk)
            {
                SharedGroupData group = traceGroup(orientStart);

                for (size_t s = 0; s < x_sizes.size(); ++s)
                {
                    const int rowsForSize = nRefineRows[s];
                    if (rowsForSize <= 0)
                        continue;

                    const double scale = x_sizes[s] / x_ref;
                    const std::vector<PreparedOrientation> *gpuPrepared = &group.prepared;
                    std::vector<PreparedOrientation> scaled;
                    double gpuScale = scale;
                    double gpuWaveIndex = waveIndex;
                    if (sharedGpuExplicitScale)
                    {
                        scaled.reserve(group.prepared.size());
                        for (const PreparedOrientation &po : group.prepared)
                            scaled.push_back(ScalePreparedOrientation(
                                po, scale, waveIndex, handlerPO->AbsorptionCoefficient()));
                        gpuPrepared = &scaled;
                        gpuScale = 1.0;
                        gpuWaveIndex = 0.0;
                    }

                    Arr2D localM(nAz + 1, rowsForSize, 4, 4); localM.ClearArr();
                    Arr2D localMns(nAz + 1, rowsForSize, 4, 4); localMns.ClearArr();
                    handlerPO->m_sphere = rowSpheres[s];
                    for (int gpuStart = 0; gpuStart < group.groupOrient; )
                    {
                        int gpuBatchSize = handlerPO->SelectGpuOrientationBatchSize(
                            *gpuPrepared, gpuStart, group.groupOrient - gpuStart);
                        int gpuEnd = std::min(gpuStart + gpuBatchSize, group.groupOrient);
                        bool ok = handlerPO->HandleOrientationsToLocalGpu(
                            *gpuPrepared, gpuStart, gpuEnd - gpuStart,
                            localM, localMns, gpuScale, gpuWaveIndex);
                        if (!ok)
                        {
                            handlerPO->m_sphere = savedSphere;
                            std::cerr << "ERROR: shared multikeq GPU global refine failed." << std::endl;
                            throw std::runtime_error("shared multikeq GPU global refine failed");
                        }
                        gpuStart = gpuEnd;
                    }
                    handlerPO->m_sphere = savedSphere;

                    if (m_mirrorGamma)
                    {
                        ApplyMirrorGammaMueller(localM, nAz, rowsForSize - 1);
                        if (computeNoShadow)
                            ApplyMirrorGammaMueller(localMns, nAz, rowsForSize - 1);
                    }
                    if (spikeDebug && m_mpiRank == 0 && spikeDebugFile.is_open())
                    {
                        const std::string label = (s < labels.size() && !labels[s].empty())
                            ? labels[s] : std::to_string(s);
                        const std::vector<int> &rows = refineRows[s];
                        for (int localRow = 0; localRow < rowsForSize; ++localRow)
                        {
                            double avg[16];
                            AzimuthAverageOutputCell(localM, nAz, rowsForSize - 1,
                                                     localRow, avg);
                            WriteSpikeDebugRow(spikeDebugFile, s, label,
                                               rows[localRow],
                                               RadToDeg(savedSphere.GetZenith(rows[localRow])),
                                               group.orientStart, group.orientEnd,
                                               nGamma, betaRange, gammaRange, avg);
                        }
                    }

                    for (int p = 0; p < nAz; ++p)
                        for (int localRow = 0; localRow < rowsForSize; ++localRow)
                        {
                            refinedM[s].insert(p, localRow, localM(p, localRow));
                            if (computeNoShadow)
                                refinedMns[s].insert(p, localRow, localMns(p, localRow));
                        }
                }
            }

            for (size_t s = 0; s < x_sizes.size(); ++s)
            {
                const int rowsForSize = nRefineRows[s];
                if (rowsForSize <= 0)
                    continue;
                const std::vector<int> &rows = refineRows[s];
                for (int p = 0; p < nAz; ++p)
                    for (int localRow = 0; localRow < rowsForSize; ++localRow)
                    {
                        const int row = rows[localRow];
                        results_M[s].replace(p, row, refinedM[s](p, localRow));
                        if (computeNoShadow)
                            results_M_ns[s].replace(p, row, refinedMns[s](p, localRow));
                    }
            }
        }
        handlerPO->m_sphere = savedSphere;
        refineSeconds = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - refineStart).count();
    }

    EraseConsoleLine(60);
    (void)t1;
    std::cout << "Shared multikeq Phase 1 (tracing/prepare once): " << std::fixed
              << std::setprecision(2) << phase1 << " s" << std::endl;
    std::cout << "Shared multikeq Phase 2 (diffraction for " << x_sizes.size()
              << " sizes, " << (handlerPO->IsGpuEnabled() ? "CUDA" : "OpenMP")
              << "): " << std::fixed << std::setprecision(2)
              << phase2 << " s" << std::endl;
    std::cout << "Shared multikeq total: " << phase1 + phase2 << " s" << std::endl;
    if (refineSeconds > 0.0)
        std::cout << "Shared multikeq FFT global refine: " << std::fixed
                  << std::setprecision(2) << refineSeconds << " s" << std::endl;

    std::string baseName = m_resultDirName;
    Arr2D savedM = handlerPO->M;
    Arr2D savedMns = handlerPO->M_noshadow;
    double savedOutEnergy = handlerPO->m_outputEnergy;
    double savedExtOt = handlerPO->m_extinctionCrossSectionOt;
    bool savedHasExtOt = handlerPO->m_hasExtinctionOt;
    for (size_t s = 0; s < x_sizes.size(); ++s)
    {
        handlerPO->M = results_M[s];
        m_incomingEnergy = results_energy[s];
        handlerPO->m_outputEnergy = results_output_energy[s];
        handlerPO->m_extinctionCrossSectionOt = results_ext_ot[s];
        handlerPO->m_hasExtinctionOt = true;
        std::string suffix = (s < labels.size() && !labels[s].empty())
            ? labels[s]
            : ("x" + std::to_string((int)x_sizes[s]));
        std::string outName = baseName + "_" + suffix;
        handlerPO->WriteMatricesToFile(outName, m_incomingEnergy);
        if (computeNoShadow)
        {
            handlerPO->M_noshadow = results_M_ns[s];
            std::swap(handlerPO->M, handlerPO->M_noshadow);
            std::string nsName = outName + "_noshadow";
            handlerPO->WriteMatricesToFile(nsName, m_incomingEnergy);
            std::swap(handlerPO->M, handlerPO->M_noshadow);
        }
    }
    handlerPO->M = savedM;
    handlerPO->M_noshadow = savedMns;
    handlerPO->m_outputEnergy = savedOutEnergy;
    handlerPO->m_extinctionCrossSectionOt = savedExtOt;
    handlerPO->m_hasExtinctionOt = savedHasExtOt;
}

void TracerPOTotal::TraceRandomMultiSizeIndependent(
    const AngleRange &betaRange,
    const AngleRange &gammaRange,
    const std::vector<double> &x_sizes,
    const std::vector<std::string> &labels,
    double x_ref)
{
    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO)
    {
        std::cerr << "Error! Handler is not HandlerPO in "
                  << "TraceRandomMultiSizeIndependent" << std::endl;
        throw std::exception();
    }
    if (x_sizes.empty())
        return;
    if (x_ref <= 0.0)
    {
        std::cerr << "ERROR: independent multikeq requires positive "
                  << "reference size." << std::endl;
        throw std::runtime_error("invalid independent multikeq reference");
    }

    const double D_ref = m_particle->MaximalDimention();
    const int nAz = handlerPO->m_sphere.nAzimuth;
    const int nZen = handlerPO->m_sphere.nZenith;
    const std::string baseName = m_resultDirName;

    if (m_mpiRank == 0)
    {
        std::cout << "Conservative multikeq: retrace each size independently "
                  << "on " << x_sizes.size() << " sizes"
                  << " (set MBS_OLDAUTOFULL_SHARED_MULTI=1 to use the "
                  << "experimental shared reference cache)" << std::endl;
    }

    auto tTotalStart = std::chrono::high_resolution_clock::now();
    for (size_t s = 0; s < x_sizes.size(); ++s)
    {
        const double D_target = D_ref * (x_sizes[s] / x_ref);
        m_particle->Resize(D_target);

        handlerPO->M = Arr2D(nAz + 1, nZen + 1, 4, 4);
        handlerPO->M_noshadow = Arr2D(nAz + 1, nZen + 1, 4, 4);
        handlerPO->M.ClearArr();
        handlerPO->M_noshadow.ClearArr();
        handlerPO->CleanJ();
        m_incomingEnergy = 0.0;
        handlerPO->m_outputEnergy = 0.0;
        handlerPO->m_extinctionCrossSectionOt = 0.0;
        handlerPO->m_hasExtinctionOt = false;

        const std::string suffix =
            (s < labels.size() && !labels[s].empty())
                ? labels[s]
                : ("x" + std::to_string((int)x_sizes[s]));
        const std::string savedName = m_resultDirName;
        m_resultDirName = baseName + "_" + suffix;

        if (m_mpiRank == 0)
        {
            std::cout << "  Independent size " << (s + 1) << "/"
                      << x_sizes.size() << ": "
                      << suffix << ", Dmax="
                      << m_particle->MaximalDimention() << std::endl;
        }
        TraceRandom(betaRange, gammaRange);
        m_resultDirName = savedName;
    }

    m_particle->Resize(D_ref);

    auto tTotalEnd = std::chrono::high_resolution_clock::now();
    const double totalSeconds =
        std::chrono::duration<double>(tTotalEnd - tTotalStart).count();
    if (m_mpiRank == 0)
        std::cout << "Conservative multikeq total: " << std::fixed
                  << std::setprecision(2) << totalSeconds << " s for "
                  << x_sizes.size() << " sizes" << std::endl;
}

void TracerPOTotal::TraceSobolMultiSize(int nOrient, double betaSym, double gammaSym,
                                         const std::vector<double> &x_sizes, double x_ref)
{
    // Simple approach: for each size, resize particle and run TraceFromSobol.
    // Uses the fully optimized HandleBeamsToLocal path (batched sincos, OpenMP).
    // Tracing is fast (<1s), dominated by diffraction which is already optimal.

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) {
        std::cerr << "Error! Handler is not HandlerPO in TraceSobolMultiSize" << std::endl;
        throw std::exception();
    }

    double D_ref = m_particle->MaximalDimention();
    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    std::string baseName = m_resultDirName;

    if (m_mpiRank == 0)
        std::cout << "MultiSize: " << x_sizes.size() << " sizes, "
                  << nOrient << " Sobol orientations each" << std::endl;

    auto t_total_start = std::chrono::high_resolution_clock::now();

    for (size_t s = 0; s < x_sizes.size(); ++s)
    {
        // Scale particle to this size
        double D_target = x_sizes[s] / x_ref * D_ref;
        m_particle->Resize(D_target);

        if (m_mpiRank == 0)
            std::cout << "  Size " << (s+1) << "/" << x_sizes.size()
                      << ": x=" << std::setprecision(12) << x_sizes[s]
                      << " (D=" << m_particle->MaximalDimention() << ")" << std::endl;

        // Reset Mueller arrays
        handlerPO->M = Arr2D(nAz+1, nZen+1, 4, 4);
        handlerPO->M_noshadow = Arr2D(nAz+1, nZen+1, 4, 4);
        m_incomingEnergy = 0;

        // Use fully optimized TraceFromSobol (batched sincos, parallel Phase 1+2)
        // Save/restore m_resultDirName since TraceFromSobol writes output
        std::string savedName = m_resultDirName;
        std::string sizeName = baseName + "_x" + SizeFileLabel(x_sizes[s]);
        m_resultDirName = sizeName;
        TraceFromSobol(nOrient, betaSym, gammaSym);
        m_resultDirName = savedName;
    }

    // Restore particle to original size
    m_particle->Resize(D_ref);

    auto t_total_end = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(t_total_end - t_total_start).count();

    if (m_mpiRank == 0)
        std::cout << "MultiSize total: " << total_sec << " s for "
                  << x_sizes.size() << " sizes" << std::endl;
}

void TracerPOTotal::TraceEulerQuadratureMultiSize(
    int nBeta, int nGamma, double betaSym, double gammaSym,
    const std::vector<double> &x_sizes, double x_ref)
{
    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) {
        std::cerr << "Error! Handler is not HandlerPO in TraceEulerQuadratureMultiSize" << std::endl;
        throw std::exception();
    }

    double D_ref = m_particle->MaximalDimention();
    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    std::string baseName = m_resultDirName;
    int nOrient = nBeta * nGamma;

    if (m_mpiRank == 0)
        std::cout << "MultiSize Euler quadrature: " << x_sizes.size()
                  << " sizes, " << nBeta << " x " << nGamma
                  << " = " << nOrient << " orientations each" << std::endl;

    auto t_total_start = std::chrono::high_resolution_clock::now();

    for (size_t s = 0; s < x_sizes.size(); ++s)
    {
        double D_target = x_sizes[s] / x_ref * D_ref;
        m_particle->Resize(D_target);

        if (m_mpiRank == 0)
            std::cout << "  Size " << (s + 1) << "/" << x_sizes.size()
                      << ": x=" << std::setprecision(12) << x_sizes[s]
                      << " (D=" << m_particle->MaximalDimention() << ")" << std::endl;

        handlerPO->M = Arr2D(nAz + 1, nZen + 1, 4, 4);
        handlerPO->M_noshadow = Arr2D(nAz + 1, nZen + 1, 4, 4);
        m_incomingEnergy = 0;

        std::string savedName = m_resultDirName;
        std::string sizeName = baseName + "_x" + SizeFileLabel(x_sizes[s]);
        m_resultDirName = sizeName;
        TraceFromEulerQuadrature(nBeta, nGamma, betaSym, gammaSym);
        m_resultDirName = savedName;
    }

    m_particle->Resize(D_ref);

    auto t_total_end = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(t_total_end - t_total_start).count();

    if (m_mpiRank == 0)
        std::cout << "MultiSize Euler quadrature total: " << total_sec
                  << " s for " << x_sizes.size() << " sizes" << std::endl;
}

void TracerPOTotal::TraceLatticeMultiSize(
    int nOrient, double betaSym, double gammaSym,
    const std::vector<double> &x_sizes, double x_ref, int generator)
{
    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO) {
        std::cerr << "Error! Handler is not HandlerPO in TraceLatticeMultiSize" << std::endl;
        throw std::exception();
    }

    double D_ref = m_particle->MaximalDimention();
    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;
    std::string baseName = m_resultDirName;

    if (m_mpiRank == 0)
        std::cout << "MultiSize lattice: " << x_sizes.size()
                  << " sizes, " << nOrient << " orientations each" << std::endl;

    auto t_total_start = std::chrono::high_resolution_clock::now();

    for (size_t s = 0; s < x_sizes.size(); ++s)
    {
        double D_target = x_sizes[s] / x_ref * D_ref;
        m_particle->Resize(D_target);

        if (m_mpiRank == 0)
            std::cout << "  Size " << (s + 1) << "/" << x_sizes.size()
                      << ": x=" << std::setprecision(12) << x_sizes[s]
                      << " (D=" << m_particle->MaximalDimention() << ")" << std::endl;

        handlerPO->M = Arr2D(nAz + 1, nZen + 1, 4, 4);
        handlerPO->M_noshadow = Arr2D(nAz + 1, nZen + 1, 4, 4);
        m_incomingEnergy = 0;

        std::string savedName = m_resultDirName;
        std::string sizeName = baseName + "_x" + SizeFileLabel(x_sizes[s]);
        m_resultDirName = sizeName;
        if (generator > 0)
            TraceFromLatticeGenerator(nOrient, generator, betaSym, gammaSym);
        else
            TraceFromLattice(nOrient, betaSym, gammaSym);
        m_resultDirName = savedName;
    }

    m_particle->Resize(D_ref);

    auto t_total_end = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(t_total_end - t_total_start).count();

    if (m_mpiRank == 0)
        std::cout << "MultiSize lattice total: " << total_sec
                  << " s for " << x_sizes.size() << " sizes" << std::endl;
}

static void BuildGaussLegendreInterval(int n, double a, double b,
                                       std::vector<double> &x,
                                       std::vector<double> &w)
{
    if (n <= 0)
        throw std::invalid_argument("Gauss-Legendre order must be positive");

    x.assign(n, 0.0);
    w.assign(n, 0.0);

    const double center = 0.5 * (a + b);
    const double half = 0.5 * (b - a);
    const int m = (n + 1) / 2;
    const double eps = 1e-14;

    for (int i = 0; i < m; ++i)
    {
        double z = std::cos(M_PI * (i + 0.75) / (n + 0.5));
        double zPrev = 0.0;
        double p1 = 0.0;
        double p2 = 0.0;
        double pp = 0.0;

        do
        {
            p1 = 1.0;
            p2 = 0.0;
            for (int j = 1; j <= n; ++j)
            {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            zPrev = z;
            z = zPrev - p1 / pp;
        } while (std::fabs(z - zPrev) > eps);

        const double nodeLo = center - half * z;
        const double nodeHi = center + half * z;
        const double weight = 2.0 * half / ((1.0 - z * z) * pp * pp);

        x[i] = nodeLo;
        x[n - 1 - i] = nodeHi;
        w[i] = weight;
        w[n - 1 - i] = weight;
    }
}

static double RadicalInverseBase2(uint32_t value)
{
    value = (value << 16) | (value >> 16);
    value = ((value & 0x00ff00ffu) << 8) | ((value & 0xff00ff00u) >> 8);
    value = ((value & 0x0f0f0f0fu) << 4) | ((value & 0xf0f0f0f0u) >> 4);
    value = ((value & 0x33333333u) << 2) | ((value & 0xccccccccu) >> 2);
    value = ((value & 0x55555555u) << 1) | ((value & 0xaaaaaaaau) >> 1);
    return (double)value / 4294967296.0;
}

static double RadicalInverseBase(uint32_t value, uint32_t base)
{
    double invBase = 1.0 / (double)base;
    double inv = invBase;
    double result = 0.0;
    while (value > 0)
    {
        result += (double)(value % base) * inv;
        value /= base;
        inv *= invBase;
    }
    return result;
}

static int GcdInt(int a, int b)
{
    a = std::abs(a);
    b = std::abs(b);
    while (b != 0)
    {
        int t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static int PickLatticeGenerator(int n)
{
    if (n <= 1)
        return 1;
    const double invPhi = 0.6180339887498948482;
    int z = (int)std::floor(n * invPhi + 0.5);
    z = std::max(1, std::min(n - 1, z));
    if ((z & 1) == 0)
        ++z;
    if (z >= n)
        z = n - 1;
    while (z > 1 && GcdInt(z, n) != 1)
        z -= 2;
    if (GcdInt(z, n) != 1)
        z = 1;
    return z;
}

static int NormalizeLatticeGenerator(int n, int generator)
{
    if (n <= 1)
        return 1;
    generator %= n;
    if (generator < 0)
        generator += n;
    if (generator == 0)
        generator = 1;
    if (GcdInt(generator, n) != 1)
        throw std::invalid_argument("lattice generator must be coprime with N");
    return generator;
}

std::vector<TracerPOTotal::WeightedOrientation>
TracerPOTotal::BuildSobolOrientations(int nOrient, double betaSym,
                                      double gammaSym,
                                      unsigned int seed) const
{
    if (nOrient <= 0)
        throw std::invalid_argument("adaptive probe orientation count must be positive");
    Sobol2D sobol(seed);
    std::vector<double> u;
    std::vector<double> v;
    sobol.generate(nOrient, u, v);
    const double cosBetaSym = std::cos(betaSym);
    const double weight = 1.0 / nOrient;
    std::vector<WeightedOrientation> result;
    result.reserve(nOrient);
    for (int i = 0; i < nOrient; ++i)
    {
        const double beta = std::acos(1.0
            - (1.0 - cosBetaSym) * u[i]);
        const double gamma = gammaSym * v[i];
        result.push_back(WeightedOrientation(beta, gamma, weight));
    }
    return result;
}

std::vector<TracerPOTotal::WeightedOrientation>
TracerPOTotal::BuildEulerOrientations(int nBeta, int nGamma,
                                      double betaSym,
                                      double gammaSym) const
{
    if (nBeta <= 0 || nGamma <= 0)
        throw std::invalid_argument("adaptive Euler counts must be positive");

    double muMin = std::cos(betaSym);
    double muMax = 1.0;
    if (muMin > muMax)
        std::swap(muMin, muMax);
    std::vector<double> mu;
    std::vector<double> muWeights;
    BuildGaussLegendreInterval(nBeta, muMin, muMax, mu, muWeights);

    const double muSpan = muMax - muMin;
    std::vector<WeightedOrientation> result;
    result.reserve((size_t)nBeta * (size_t)nGamma);
    for (int ib = 0; ib < nBeta; ++ib)
    {
        const double beta = std::acos(
            std::max(-1.0, std::min(1.0, mu[ib])));
        const double betaWeight = muSpan > 0.0
            ? muWeights[ib] / muSpan : 1.0 / nBeta;
        for (int ig = 0; ig < nGamma; ++ig)
        {
            const double gamma = gammaSym * (ig + 0.5) / nGamma;
            result.push_back(WeightedOrientation(
                beta, gamma, betaWeight / nGamma));
        }
    }

    double weightSum = 0.0;
    for (const WeightedOrientation &orientation : result)
        weightSum += orientation.weight;
    if (weightSum > 0.0)
        for (WeightedOrientation &orientation : result)
            orientation.weight /= weightSum;
    return result;
}

void TracerPOTotal::PrepareOrientationProbe(
    const std::vector<WeightedOrientation> &orientations,
    std::vector<PreparedOrientation> &prepared,
    double &incomingEnergy, double &outputEnergy)
{
    HandlerPO *handler = dynamic_cast<HandlerPO*>(m_handler);
    if (!handler)
        throw std::runtime_error("adaptive PO probe requires HandlerPO");

    prepared.assign(orientations.size(), PreparedOrientation());
    std::vector<double> incoming(orientations.size(), 0.0);
    std::vector<double> outgoing(orientations.size(), 0.0);
    m_scattering->PrepareForParallelTrace();
    ParallelExceptionState parallelError;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        Particle localParticle = *m_particle;
        Scattering *localScatter = m_scattering->CloneFor(
            &localParticle, &m_incidentLight);
        HandlerPO localHandler(&localParticle, &m_incidentLight,
                               handler->nTheta, m_scattering->m_wave);
        localHandler.ConfigureForThreadLocalPrepare(*handler, localScatter);
        std::vector<Beam> beams;

#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 4)
#endif
        for (int i = 0; i < (int)orientations.size(); ++i)
        {
            if (parallelError.Failed())
                continue;
            try
            {
            const WeightedOrientation &orientation = orientations[i];
            if (orientation.useQuaternion)
                localParticle.RotateQuaternion(
                    orientation.qx, orientation.qy,
                    orientation.qz, orientation.qw);
            else
                localParticle.Rotate(orientation.beta, orientation.gamma,
                                     orientation.alpha);
            if (!shadowOff)
                localScatter->FormShadowBeam(beams);
            const bool ok = localScatter->ScatterLight(0, 0, beams);
            if (ok)
                localHandler.PrepareBeams(beams, orientation.weight,
                                          prepared[i]);
            else
                prepared[i].sinZenith = orientation.weight;
            incoming[i] = localScatter->GetIncedentEnergy()
                * orientation.weight;
            outgoing[i] = PreparedOutputEnergy(
                prepared[i], 1.0, localHandler.AbsorptionCoefficient());
            beams.clear();
            }
            catch (...)
            {
                parallelError.Capture();
            }
        }
        delete localScatter;
    }
    parallelError.Rethrow();

    incomingEnergy = 0.0;
    outputEnergy = 0.0;
    for (size_t i = 0; i < orientations.size(); ++i)
    {
        incomingEnergy += incoming[i];
        outputEnergy += outgoing[i];
    }
}

void TracerPOTotal::TraceFromSobol(int nOrient, double betaSym, double gammaSym)
{
    TraceFromSobolSeed(nOrient, 42u, betaSym, gammaSym);
}

void TracerPOTotal::TraceFromSO3Quaternion(int nOrient)
{
    if (nOrient <= 0)
        throw std::invalid_argument("SO(3) quaternion orientation count must be positive");

    const double weight = 1.0 / nOrient;
    std::vector<WeightedOrientation> orientations;
    orientations.reserve(nOrient);

    for (int i = 0; i < nOrient; ++i)
    {
        const double u1 = (i + 0.5) / nOrient;
        const double u2 = RadicalInverseBase2((uint32_t)i);
        const double u3 = RadicalInverseBase((uint32_t)i, 3);
        const double r1 = std::sqrt(std::max(0.0, 1.0 - u1));
        const double r2 = std::sqrt(u1);
        const double a = 2.0 * M_PI * u2;
        const double b = 2.0 * M_PI * u3;

        // Shoemake uniform unit quaternion on SO(3).  Stored directly and
        // converted to the particle rotation matrix immediately before tracing.
        const double qx = r1 * std::sin(a);
        const double qy = r1 * std::cos(a);
        const double qz = r2 * std::sin(b);
        const double qw = r2 * std::cos(b);
        orientations.push_back(WeightedOrientation(qx, qy, qz, qw, weight));
    }

    if (m_mpiRank == 0)
        std::cout << "SO(3) quaternion Hammersley: " << nOrient
                  << " orientations, full rotation group, no beta/gamma symmetry"
                  << std::endl;

    TraceWeightedOrientations(orientations, "SO(3) quaternion", M_PI, 2.0*M_PI);
}

void TracerPOTotal::TraceFromSobolSeed(int nOrient, unsigned int seed,
                                       double betaSym, double gammaSym)
{
    // Generate Sobol orientations mapped to [0, betaSym] x [0, gammaSym]
    // using the correct solid-angle measure:
    //   beta = arccos(1 - (1-cos(betaSym)) * u)  for uniform dOmega
    //   gamma = gammaSym * v
    if (nOrient <= 0)
        throw std::invalid_argument("Sobol orientation count must be positive");

    Sobol2D sobol(seed);
    std::vector<double> su, sv;
    sobol.generate(nOrient, su, sv);

    double cosBetaSym = cos(betaSym);
    std::vector<WeightedOrientation> orientations(nOrient);
    double weight = 1.0 / nOrient;
    for (int i = 0; i < nOrient; ++i)
    {
        double beta  = acos(1.0 - (1.0 - cosBetaSym) * su[i]);
        double gamma = gammaSym * sv[i];
        orientations[i] = {beta, gamma, weight};
    }

    std::ostringstream label;
    label << "Sobol Owen(seed=" << seed << ")";
    TraceWeightedOrientations(orientations, label.str(), betaSym, gammaSym);
}

void TracerPOTotal::TraceFromSobolRing(int nBeta, int nGamma,
                                       double betaSym, double gammaSym)
{
    if (nBeta <= 0 || nGamma <= 0)
        throw std::invalid_argument("--sobol_ring requires positive Nbeta and Ngamma");

    Sobol2D sobol(42);
    std::vector<double> su, sv;
    sobol.generate(nBeta, su, sv);

    double cosBetaSym = cos(betaSym);
    int nOrient = nBeta * nGamma;
    double weight = 1.0 / nOrient;
    std::vector<WeightedOrientation> orientations;
    orientations.reserve(nOrient);

    for (int ib = 0; ib < nBeta; ++ib)
    {
        double beta = acos(1.0 - (1.0 - cosBetaSym) * su[ib]);
        double gammaShift = sv[ib];
        for (int ig = 0; ig < nGamma; ++ig)
        {
            double gammaUnit = (ig + gammaShift) / nGamma;
            gammaUnit -= std::floor(gammaUnit);
            double gamma = gammaSym * gammaUnit;
            orientations.push_back({beta, gamma, weight});
        }
    }

    if (m_mpiRank == 0)
    {
        std::cout << "Sobol ring: Sobol beta x shifted periodic gamma, Nbeta="
                  << nBeta << ", Ngamma=" << nGamma
                  << ", total=" << orientations.size() << " orientations"
                  << std::endl;
    }

    TraceWeightedOrientations(orientations, "Sobol ring", betaSym, gammaSym);
}

void TracerPOTotal::TraceFromHammersley(int nOrient, double betaSym,
                                        double gammaSym)
{
    if (nOrient <= 0)
        throw std::invalid_argument("Hammersley orientation count must be positive");

    double cosBetaSym = cos(betaSym);
    double weight = 1.0 / nOrient;
    std::vector<WeightedOrientation> orientations;
    orientations.reserve(nOrient);

    for (int i = 0; i < nOrient; ++i)
    {
        double u = (i + 0.5) / nOrient;
        double v = RadicalInverseBase2((uint32_t)i);
        double beta = acos(1.0 - (1.0 - cosBetaSym) * u);
        double gamma = gammaSym * v;
        orientations.push_back({beta, gamma, weight});
    }

    if (m_mpiRank == 0)
        std::cout << "Hammersley: " << nOrient << " orientations" << std::endl;

    TraceWeightedOrientations(orientations, "Hammersley", betaSym, gammaSym);
}

void TracerPOTotal::TraceFromLattice(int nOrient, double betaSym,
                                     double gammaSym)
{
    TraceFromLatticeGenerator(nOrient, PickLatticeGenerator(nOrient),
                              betaSym, gammaSym);
}

void TracerPOTotal::TraceFromLatticeGenerator(int nOrient, int generator,
                                              double betaSym,
                                              double gammaSym)
{
    if (nOrient <= 0)
        throw std::invalid_argument("Lattice orientation count must be positive");

    double cosBetaSym = cos(betaSym);
    double weight = 1.0 / nOrient;
    generator = NormalizeLatticeGenerator(nOrient, generator);
    std::vector<WeightedOrientation> orientations;
    orientations.reserve(nOrient);

    for (int i = 0; i < nOrient; ++i)
    {
        double u = (i + 0.5) / nOrient;
        double v = (((long long)i * generator) % nOrient + 0.5) / nOrient;
        double beta = acos(1.0 - (1.0 - cosBetaSym) * u);
        double gamma = gammaSym * v;
        orientations.push_back({beta, gamma, weight});
    }

    if (m_mpiRank == 0)
        std::cout << "Rank-1 lattice: " << nOrient
                  << " orientations, generator=" << generator << std::endl;

    TraceWeightedOrientations(orientations, "Rank-1 lattice", betaSym, gammaSym);
}

void TracerPOTotal::TraceFromEulerQuadrature(int nBeta, int nGamma,
                                             double betaSym, double gammaSym)
{
    if (nBeta <= 0 || nGamma <= 0)
        throw std::invalid_argument("--euler_quad requires positive Nbeta and Ngamma");

    std::vector<WeightedOrientation> orientations = BuildEulerOrientations(
        nBeta, nGamma, betaSym, gammaSym);

    if (m_mpiRank == 0)
    {
        std::cout << "Euler quadrature: Gauss-Legendre in cos(beta) x "
                  << "periodic midpoint gamma, Nbeta=" << nBeta
                  << ", Ngamma=" << nGamma
                  << ", algebraic beta order=" << (2 * nBeta - 1)
                  << ", total=" << orientations.size() << " orientations"
                  << std::endl;
    }

    TraceWeightedOrientations(orientations, "Euler quadrature", betaSym, gammaSym);
}

void TracerPOTotal::TraceFromEulerAdaptiveGamma(int nBeta, int nGammaMax,
                                                double betaSym, double gammaSym)
{
    if (nBeta <= 0 || nGammaMax <= 0)
        throw std::invalid_argument("--euler_adapt requires positive Nbeta and NgammaMax");

    double muMin = std::cos(betaSym);
    double muMax = 1.0;
    if (muMin > muMax)
        std::swap(muMin, muMax);

    std::vector<double> mu;
    std::vector<double> muWeights;
    BuildGaussLegendreInterval(nBeta, muMin, muMax, mu, muWeights);

    const double muSpan = muMax - muMin;
    std::vector<double> betas(nBeta, 0.0);
    double maxRing = 0.0;
    for (int ib = 0; ib < nBeta; ++ib)
    {
        double clampedMu = std::max(-1.0, std::min(1.0, mu[ib]));
        betas[ib] = std::acos(clampedMu);
        maxRing = std::max(maxRing, std::sin(betas[ib]));
    }
    if (maxRing <= DBL_EPSILON)
        maxRing = 1.0;

    int minGamma = std::max(6, nGammaMax / 4);
    minGamma = ((minGamma + 5) / 6) * 6;
    nGammaMax = ((nGammaMax + 5) / 6) * 6;
    if (nGammaMax < minGamma)
        nGammaMax = minGamma;

    std::vector<WeightedOrientation> orientations;
    orientations.reserve((size_t)nBeta * (size_t)nGammaMax);
    std::vector<int> gammaCounts(nBeta, 0);

    for (int ib = 0; ib < nBeta; ++ib)
    {
        double betaWeight = (muSpan > 0.0)
            ? muWeights[ib] / muSpan
            : 1.0 / nBeta;

        int nGamma = (int)std::ceil(nGammaMax * std::sin(betas[ib]) / maxRing);
        nGamma = std::max(minGamma, nGamma);
        nGamma = ((nGamma + 5) / 6) * 6;
        nGamma = std::min(nGammaMax, nGamma);
        gammaCounts[ib] = nGamma;

        for (int ig = 0; ig < nGamma; ++ig)
        {
            double gamma = gammaSym * (ig + 0.5) / nGamma;
            orientations.push_back(WeightedOrientation(
                betas[ib], gamma, betaWeight / nGamma));
        }
    }

    double weightSum = 0.0;
    for (const WeightedOrientation &orientation : orientations)
        weightSum += orientation.weight;
    if (weightSum > 0.0)
        for (WeightedOrientation &orientation : orientations)
            orientation.weight /= weightSum;

    if (m_mpiRank == 0)
    {
        int minUsed = gammaCounts.empty() ? 0 : gammaCounts[0];
        int maxUsed = gammaCounts.empty() ? 0 : gammaCounts[0];
        long long sumUsed = 0;
        for (int n : gammaCounts)
        {
            minUsed = std::min(minUsed, n);
            maxUsed = std::max(maxUsed, n);
            sumUsed += n;
        }
        std::cout << "Euler adaptive gamma: Gauss-Legendre in cos(beta), "
                  << "Nbeta=" << nBeta << ", NgammaMax=" << nGammaMax
                  << ", Ngamma range=" << minUsed << ".." << maxUsed
                  << ", total=" << orientations.size()
                  << " orientations, saved="
                  << (100.0 * (1.0 - (double)orientations.size()
                              / std::max(1.0, (double)nBeta * nGammaMax)))
                  << "% vs fixed" << std::endl;
    }

    TraceWeightedOrientations(orientations, "Euler adaptive gamma",
                              betaSym, gammaSym);
}

void TracerPOTotal::TraceWeightedOrientations(
    const std::vector<WeightedOrientation> &orientations,
    const std::string &label, double betaSym, double gammaSym)
{
    int nOrient = (int)orientations.size();
    if (nOrient <= 0)
        throw std::invalid_argument("orientation list must be non-empty");

    // MPI: each rank processes a subset of orientations
    int myStart = m_mpiRank * nOrient / m_mpiSize;
    int myEnd = (m_mpiRank + 1) * nOrient / m_mpiSize;
    int myCount = myEnd - myStart;

    if (m_mpiRank == 0)
        std::cout << label << ": " << nOrient << " orientations"
                  << (m_mpiSize > 1 ? " (" + std::to_string(myCount) + " per rank, " + std::to_string(m_mpiSize) + " ranks)" : "")
                  << ", beta_sym=" << RadToDeg(betaSym) << " deg, gamma_sym="
                  << RadToDeg(gammaSym) << " deg" << std::endl;

    CalcTimer timer;
    timer.Start();
    if (m_mpiRank == 0) OutputStartTime(timer);

    m_handler->SetNormIndex(1);

    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO)
    {
        std::cerr << "Error! Handler is not HandlerPO in TraceWeightedOrientations" << std::endl;
        throw std::exception();
    }
    const bool computeNoShadow = handlerPO->ComputeNoShadow();

    m_incomingEnergy = 0;
    handlerPO->m_outputEnergy = 0;
    handlerPO->m_extinctionCrossSectionOt = 0;
    handlerPO->m_hasExtinctionOt = false;
    handlerPO->M.ClearArr();
    handlerPO->M_noshadow.ClearArr();
    handlerPO->CleanJ();

    // Phase 1 traces and preprocesses orientations in parallel using
    // thread-local Particle/Scattering/HandlerPO state. Phase 2 then performs
    // the diffraction loops in parallel with local Mueller accumulators.
    int nAz = handlerPO->m_sphere.nAzimuth;
    int nZen = handlerPO->m_sphere.nZenith;

    // Chunked streaming: auto-size chunks based on available RAM. The value is
    // capped by MBS_HOST_MEM_BUDGET_MB when launched by the parallel scheduler.
    long long availMB = EffectiveMemAvailableMb();
    long long beamBudget = std::max(100LL, availMB / 2);
    int chunkSize = std::max(32, std::min(4096, std::min(myCount, (int)(beamBudget * 1024 / 350))));
    if (m_sobolChunkSize > 0)
        chunkSize = std::max(1, std::min(chunkSize, m_sobolChunkSize));
    int nChunks = (myCount + chunkSize - 1) / chunkSize;
    int nThreads = 1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        nThreads = omp_get_num_threads();
    }
#endif

    if (m_mpiRank == 0)
    {
        std::ostringstream log;
        log << "Memory: " << availMB << " MB available";
        if (HostMemoryBudgetOverrideKb() > 0)
            log << " (shared budget)";
        log << ", chunk="
            << chunkSize << " orientations (" << nChunks << " chunks), "
            << nThreads << " threads";
        std::cerr << log.str() << std::endl;
        AppendTextLog(log.str() + "\n");
    }

    double phase1_total = 0, phase2_total = 0;
    std::vector<Beam> outBeams;
    long long count = 0;
    m_scattering->PrepareForParallelTrace();

    for (int chunk = 0; chunk < nChunks; ++chunk)
    {
        int iStart = myStart + chunk * chunkSize;
        int iEnd = std::min(iStart + chunkSize, myEnd);
        int thisChunk = iEnd - iStart;

        // Phase 1: trace and preprocess this chunk in parallel. Each thread
        // has its own Particle, Scattering and HandlerPO scratch state.
        auto tp1 = std::chrono::high_resolution_clock::now();
        std::vector<PreparedOrientation> chunkPrepared(thisChunk);
        std::vector<double> chunkEnergies(thisChunk, 0);
        std::vector<double> chunkOutputEnergies(thisChunk, 0);
        ParallelExceptionState parallelError;

        #pragma omp parallel
        {
            // Thread-local copies of all mutable tracing/preprocessing state.
            Particle localParticle = *m_particle;
            Scattering *localScatter = m_scattering->CloneFor(&localParticle, &m_incidentLight);

            HandlerPO localHandler(&localParticle, &m_incidentLight,
                                   handlerPO->nTheta, m_scattering->m_wave);
            localHandler.ConfigureForThreadLocalPrepare(*handlerPO, localScatter);

            std::vector<Beam> localBeams;

            #pragma omp for schedule(dynamic, 4)
            for (int i = 0; i < thisChunk; ++i)
            {
                if (parallelError.Failed())
                    continue;
                try
                {
                int idx = iStart + i;
                const WeightedOrientation &orientation = orientations[idx];
                double weight = orientation.weight;
                if (orientation.useQuaternion)
                    localParticle.RotateQuaternion(orientation.qx, orientation.qy,
                                                   orientation.qz, orientation.qw);
                else
                    localParticle.Rotate(orientation.beta, orientation.gamma,
                                         orientation.alpha);
                if (!shadowOff) localScatter->FormShadowBeam(localBeams);
                bool ok = localScatter->ScatterLight(0, 0, localBeams);
                if (ok)
                {
                    double beforeOutput = localHandler.m_outputEnergy;
                    localHandler.PrepareBeams(localBeams, weight, chunkPrepared[i]);
                    chunkOutputEnergies[i] = localHandler.m_outputEnergy - beforeOutput;
                }
                else
                {
                    chunkPrepared[i].sinZenith = weight;
                }
                chunkEnergies[i] = localScatter->GetIncedentEnergy() * weight;
                localBeams.clear();
                }
                catch (...)
                {
                    parallelError.Capture();
                }
            }

            delete localScatter;
        }
        parallelError.Rethrow();

        // Accumulate scalar counters sequentially after the parallel region.
        for (int i = 0; i < thisChunk; ++i) {
            m_incomingEnergy += chunkEnergies[i];
            handlerPO->m_outputEnergy += chunkOutputEnergies[i];
            handlerPO->m_extinctionCrossSectionOt +=
                chunkPrepared[i].extinctionOt;
            handlerPO->m_hasExtinctionOt = true;
            count++;
        }
        phase1_total += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tp1).count();

        // Phase 2: parallel diffraction for this chunk
        auto tp2 = std::chrono::high_resolution_clock::now();

        if (handlerPO->IsGpuEnabled())
        {
            Arr2D localM(nAz + 1, nZen + 1, 4, 4);
            localM.ClearArr();
            Arr2D localM_ns(nAz + 1, nZen + 1, 4, 4);
            localM_ns.ClearArr();

            for (int gpuStart = 0; gpuStart < thisChunk; )
            {
                int gpuBatchSize = handlerPO->SelectGpuOrientationBatchSize(
                    chunkPrepared, gpuStart, thisChunk - gpuStart);
                int gpuEnd = std::min(gpuStart + gpuBatchSize, thisChunk);
                bool ok = handlerPO->IsFftEnabled()
                    ? handlerPO->HandleOrientationsToLocalGpuFftPhi(
                        chunkPrepared, gpuStart, gpuEnd - gpuStart,
                        localM, localM_ns)
                    : handlerPO->HandleOrientationsToLocalGpu(
                        chunkPrepared, gpuStart, gpuEnd - gpuStart,
                        localM, localM_ns);
                if (!ok)
                {
                    std::cerr << "ERROR: --gpu requested but GPU diffraction backend "
                              << "could not process this chunk." << std::endl;
                    throw std::runtime_error("GPU diffraction backend failed");
                }
                gpuStart = gpuEnd;
            }

            for (int p = 0; p < nAz; ++p)
                for (int t = 0; t <= nZen; ++t)
                {
                    handlerPO->M.insert(p, t, localM(p, t));
                    if (computeNoShadow)
                        handlerPO->M_noshadow.insert(
                            p, t, localM_ns(p, t));
                }
        }
        else
        {
            #pragma omp parallel
            {
                Arr2D localM(nAz + 1, nZen + 1, 4, 4);
                localM.ClearArr();
                Arr2D localM_ns(nAz + 1, nZen + 1, 4, 4);
                localM_ns.ClearArr();
                std::vector<Arr2DC> localJ, localJ_ns;
                if (handlerPO->isCoh) {
                    Arr2DC tmp(nAz + 1, nZen + 1, 2, 2);
                    tmp.ClearArr();
                    localJ.push_back(tmp);
                    Arr2DC tmp2(nAz + 1, nZen + 1, 2, 2);
                    tmp2.ClearArr();
                    localJ_ns.push_back(tmp2);
                }

                #pragma omp for schedule(dynamic, 1)
                for (int i = 0; i < thisChunk; ++i)
                {
                    if (!chunkPrepared[i].beams.empty())
                        handlerPO->HandleBeamsToLocal(chunkPrepared[i], localM, localJ,
                                                       handlerPO->isCoh ? &localJ_ns : nullptr);
                    if (handlerPO->isCoh && !localJ.empty()) {
                        double w = chunkPrepared[i].sinZenith;
                        HandlerPO::AddToMuellerLocal(localJ, w, localM, nAz, nZen);
                        HandlerPO::AddToMuellerLocal(localJ_ns, w, localM_ns, nAz, nZen);
                        localJ[0].ClearArr();
                        localJ_ns[0].ClearArr();
                    }
                }

                #pragma omp critical
                {
                    for (int p = 0; p < nAz; ++p)
                        for (int t = 0; t <= nZen; ++t) {
                            handlerPO->M.insert(p, t, localM(p, t));
                            if (computeNoShadow)
                                handlerPO->M_noshadow.insert(
                                    p, t, localM_ns(p, t));
                        }
                }
            }
        }
        phase2_total += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tp2).count();

        if (m_mpiRank == 0)
            OutputProgress(nOrient, count, iEnd - 1, chunk + 1, timer, -1);

        chunkPrepared.clear();
        chunkPrepared.shrink_to_fit();
    }

    // MPI: reduce Mueller matrices from all ranks to rank 0
    MPI_ReduceMueller(handlerPO, nAz, nZen, m_incomingEnergy, m_mpiRank);
    if (m_mpiRank == 0 && m_mirrorGamma)
    {
        ApplyMirrorGammaMueller(handlerPO->M, nAz, nZen);
        if (handlerPO->ComputeNoShadow())
            ApplyMirrorGammaMueller(handlerPO->M_noshadow, nAz, nZen);
    }

    EraseConsoleLine(60);
    if (m_mpiRank == 0) {
        std::cout << "Phase 1 (tracing): " << std::fixed
                  << std::setprecision(2) << phase1_total << " s" << std::endl;
        std::cout << "Phase 2 (diffraction, "
                  << (handlerPO->IsGpuEnabled() ? "CUDA" : "OpenMP")
                  << "): " << phase2_total << " s" << std::endl;
        std::cout << "Total: " << phase1_total + phase2_total << " s" << std::endl;

        m_handler->WriteTotalMatricesToFile(m_resultDirName);
        m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);

        if (handlerPO->ComputeNoShadow())
        {
            std::swap(handlerPO->M, handlerPO->M_noshadow);
            std::string nsName = m_resultDirName + "_noshadow";
            handlerPO->WriteMatricesToFile(nsName, m_incomingEnergy);
            std::swap(handlerPO->M, handlerPO->M_noshadow);
        }

        OutputStatisticsPO(timer, nOrient, m_resultDirName);
    }
}

double TracerPOTotal::TraceFromSobolVariablePhi(int nOrient, double betaSym,
                                                double gammaSym,
                                                const std::vector<int> &rowNphi,
                                                const std::vector<int> &rowNorient,
                                                int outputNphi,
                                                double fftEps,
                                                const std::vector<unsigned int> &owenSeeds,
                                                const std::vector<double> &rowBeamCutoff)
{
    HandlerPO *handlerPO = dynamic_cast<HandlerPO*>(m_handler);
    if (!handlerPO)
    {
        std::cerr << "Error! Handler is not HandlerPO in TraceFromSobolVariablePhi" << std::endl;
        throw std::exception();
    }

    ScatteringRange outputSphere = handlerPO->m_sphere;
    outputSphere.nAzimuth = outputNphi;
    outputSphere.azinuthStep = M_2PI / outputNphi;
    outputSphere.ComputeSphereDirections(m_incidentLight);
    handlerPO->SetScatteringSphere(outputSphere);

    const int nZen = outputSphere.nZenith;
    std::map<AutoFullThetaWorkKey, std::vector<int>> rowsByWork;
    const bool canUseFftAverage =
        AutoFullLowPhiAverageEnabled(handlerPO, fftEps);
    const bool directFullPhiEnabled =
        AutoFullDirectFullPhiEnabled(handlerPO);
    const int minDirectPhi = canUseFftAverage
        ? AutoFullLowPhiDirectFloor(outputNphi, fftEps)
        : 2;
    const int minRowPhi = canUseFftAverage
        ? std::min(outputNphi, RoundUpToMultiple(minDirectPhi, 6))
        : 2;
    const bool thetaZonesEnabled = AutoFullThetaZonesEnabled();
    const bool rowOrientEnabled = AutoFullRowOrientEnabled();
    const double baseImportanceCutoff =
        handlerPO->m_beamCutoffImportanceRel > 0.0
            ? handlerPO->m_beamCutoffImportanceRel
            : 0.0;
    std::vector<int> rowPhiUsed(nZen + 1, outputNphi);
    std::vector<int> rowZone(nZen + 1, 1);
    std::vector<int> rowOrientLimit(nZen + 1, nOrient);
    std::vector<int> rowOrientSource(nZen + 1, 0);
    std::vector<double> rowExtraImportanceCutoff(nZen + 1, 0.0);
    std::vector<int> rowDirectFullPhi(nZen + 1, 0);
    auto lookupRowOrient = [&](int row, int fallbackZone, int *source) {
        if (!rowOrientEnabled)
        {
            if (source) *source = 0;
            return AutoFullThetaZoneOrientLimit(nOrient, fallbackZone);
        }

        if (row < (int)rowNorient.size() && rowNorient[row] > 0)
        {
            if (source) *source = 1;
            return rowNorient[row];
        }

        int left = 0;
        for (int i = std::min(row - 1, (int)rowNorient.size() - 1);
             i >= 0; --i)
        {
            if (rowNorient[i] > 0)
            {
                left = rowNorient[i];
                break;
            }
        }
        int right = 0;
        for (int i = row + 1; i < (int)rowNorient.size(); ++i)
        {
            if (rowNorient[i] > 0)
            {
                right = rowNorient[i];
                break;
            }
        }
        int value = std::max(left, right);
        if (value > 0)
        {
            if (thetaZonesEnabled)
                value = std::max(value,
                    AutoFullThetaZoneOrientLimit(nOrient, fallbackZone));
            if (source) *source = 2;
            return value;
        }

        if (source) *source = 0;
        return AutoFullThetaZoneOrientLimit(nOrient, fallbackZone);
    };
    for (int t = 0; t <= nZen; ++t)
    {
        int nPhi = (t < (int)rowNphi.size() && rowNphi[t] > 0)
            ? rowNphi[t] : outputNphi;
        nPhi = std::max(2, std::min(nPhi, outputNphi));
        nPhi = RoundUpToMultiple(nPhi, 6);
        nPhi = std::max(2, std::min(nPhi, outputNphi));
        nPhi = std::max(nPhi, minRowPhi);

        const double thetaDeg = RadToDeg(outputSphere.GetZenith(t));
        const int physicalZone = AutoFullThetaZone(thetaDeg);
        const int fallbackZone = thetaZonesEnabled ? physicalZone : 2;
        int orientSource = 0;
        int orientLimit = lookupRowOrient(t, fallbackZone, &orientSource);
        orientLimit = std::max(64, std::min(orientLimit, nOrient));
        orientLimit = RoundDownToMultiple(orientLimit, 16);
        orientLimit = std::max(64, std::min(orientLimit, nOrient));

        double extraImportanceCutoff = 0.0;
        if (t < (int)rowBeamCutoff.size()
            && rowBeamCutoff[t] > baseImportanceCutoff * (1.0 + 1e-12))
        {
            extraImportanceCutoff = rowBeamCutoff[t];
        }
        else if (thetaZonesEnabled && baseImportanceCutoff > 0.0)
        {
            const double zoneCutoff =
                baseImportanceCutoff
                * AutoFullThetaZoneCutoffMultiplier(physicalZone);
            if (zoneCutoff > baseImportanceCutoff * (1.0 + 1e-12))
                extraImportanceCutoff = zoneCutoff;
        }

        rowPhiUsed[t] = nPhi;
        rowZone[t] = physicalZone;
        rowOrientLimit[t] = orientLimit;
        rowOrientSource[t] = orientSource;
        rowExtraImportanceCutoff[t] = extraImportanceCutoff;
        rowDirectFullPhi[t] =
            (directFullPhiEnabled && nPhi >= outputNphi) ? 1 : 0;
        const int workZone = extraImportanceCutoff > 0.0 ? physicalZone : -1;
        AutoFullThetaWorkKey key = {
            nPhi, orientLimit, workZone, extraImportanceCutoff,
            rowDirectFullPhi[t] != 0
        };
        rowsByWork[key].push_back(t);
    }
    rowsByWork.clear();
    if (nZen >= 1)
    {
        auto safeSharedCutoff = [](double a, double b)
        {
            if (a <= 0.0 || b <= 0.0)
                return 0.0;
            return std::min(a, b);
        };
        auto tiePoleToNeighbor = [&](int pole, int neighbor)
        {
            const int sharedPhi =
                std::max(rowPhiUsed[pole], rowPhiUsed[neighbor]);
            const int sharedOrient =
                std::max(rowOrientLimit[pole], rowOrientLimit[neighbor]);
            const double sharedCutoff = safeSharedCutoff(
                rowExtraImportanceCutoff[pole],
                rowExtraImportanceCutoff[neighbor]);
            rowPhiUsed[pole] = rowPhiUsed[neighbor] = sharedPhi;
            rowOrientLimit[pole] = rowOrientLimit[neighbor] = sharedOrient;
            rowExtraImportanceCutoff[pole] =
                rowExtraImportanceCutoff[neighbor] = sharedCutoff;
            const int sharedDirectFull =
                (directFullPhiEnabled && sharedPhi >= outputNphi) ? 1 : 0;
            rowDirectFullPhi[pole] = rowDirectFullPhi[neighbor] =
                sharedDirectFull;
        };
        tiePoleToNeighbor(0, 1);
        tiePoleToNeighbor(nZen, nZen - 1);
    }
    for (int t = 0; t <= nZen; ++t)
    {
        const int workZone =
            rowExtraImportanceCutoff[t] > 0.0 ? rowZone[t] : -1;
        AutoFullThetaWorkKey key = {
            rowPhiUsed[t], rowOrientLimit[t], workZone,
            rowExtraImportanceCutoff[t],
            rowDirectFullPhi[t] != 0
        };
        rowsByWork[key].push_back(t);
    }

    double work = 0.0;
    double directWork = 0.0;
    double orientWork = 0.0;
    int directFullPhiRows = 0;
    for (const auto &entry : rowsByWork)
    {
        const int nPhi = entry.first.nPhi;
        const int orientLimit = entry.first.nOrient;
        if (entry.first.directFullPhi)
            directFullPhiRows += (int)entry.second.size();
        work += (double)nPhi * (double)entry.second.size();
        int directPhi = AutoFullVariableCalcPhi(
            nPhi, fftEps, minDirectPhi, canUseFftAverage,
            entry.first.directFullPhi);
        directWork += (double)directPhi * (double)entry.second.size();
        orientWork += (double)directPhi * (double)entry.second.size()
                    * ((double)orientLimit / std::max(1, nOrient));
    }
    double fullWork = (double)outputNphi * (double)(nZen + 1);
    std::cout << "Autofull variable-phi final: " << rowsByWork.size()
              << " theta work groups, estimated phi-row work "
              << std::fixed << std::setprecision(1)
              << (100.0 * work / std::max(fullWork, 1.0)) << "% of rectangular"
              << std::endl;
    if (canUseFftAverage)
    {
        std::cout << "  FFT-average direct phi work "
                  << std::fixed << std::setprecision(1)
                  << (100.0 * directWork / std::max(fullWork, 1.0))
                  << "% of rectangular"
                  << ", direct N_phi floor=" << minDirectPhi
                  << " (MBS_AUTOFULL_FINAL_LOWPHI_AVG=0 disables)"
                  << std::endl;
        if (directFullPhiRows > 0)
        {
            std::cout << "  Direct full-phi fallback rows: "
                      << directFullPhiRows << "/" << (nZen + 1)
                      << " (MBS_AUTOFULL_DIRECT_FULL_NPHI=0 disables)"
                      << std::endl;
        }
    }
    if (thetaZonesEnabled || rowOrientEnabled)
    {
        std::cout << "  Theta-zone orientation/direct work "
                  << std::fixed << std::setprecision(1)
                  << (100.0 * orientWork / std::max(fullWork, 1.0))
                  << "% of full rectangular-orientation work"
                  << " (MBS_AUTOFULL_THETA_ZONES=0 disables zones,"
                  << " MBS_AUTOFULL_ROW_ORIENT=0 disables adaptive row N)"
                  << std::endl;
    }
    for (const auto &entry : rowsByWork)
    {
        const int nPhi = entry.first.nPhi;
        int directPhi = AutoFullVariableCalcPhi(
            nPhi, fftEps, minDirectPhi, canUseFftAverage,
            entry.first.directFullPhi);
        double extraCutoff = 0.0;
        if (!entry.second.empty())
            extraCutoff = rowExtraImportanceCutoff[entry.second.front()];
        std::cout << "  zone=" << AutoFullThetaZoneName(entry.first.zone)
                  << " N_orient=" << entry.first.nOrient
                  << " N_phi=" << nPhi << ": "
                  << entry.second.size() << " theta rows";
        if (entry.first.directFullPhi)
            std::cout << ", direct full N_phi";
        if (directPhi != nPhi)
            std::cout << ", averaged from direct N_phi=" << directPhi;
        if (extraCutoff > 0.0)
            std::cout << ", beam_importance_cutoff=" << std::scientific
                      << std::setprecision(1) << extraCutoff << std::fixed
                      << std::setprecision(1);
        std::cout << std::endl;
    }
    if (m_mpiRank == 0)
    {
        const std::string thetaWorkName =
            m_resultDirName + "_autofull_theta_work.dat";
        std::ofstream thetaWork(thetaWorkName.c_str(), std::ios::out);
        thetaWork << std::setprecision(10);
        thetaWork << "# theta_deg zone needed_Nphi needed_Norient"
                  << " orient_source extra_beam_importance_cutoff_rel"
                  << " direct_full_phi\n";
        for (int t = 0; t <= nZen; ++t)
        {
            const char *source = rowOrientSource[t] == 1 ? "adaptive"
                : (rowOrientSource[t] == 2 ? "neighbor" : "fallback");
            thetaWork << RadToDeg(outputSphere.GetZenith(t)) << ' '
                      << AutoFullThetaZoneName(rowZone[t]) << ' '
                      << rowPhiUsed[t] << ' '
                      << rowOrientLimit[t] << ' '
                      << source << ' '
                      << rowExtraImportanceCutoff[t] << ' '
                      << rowDirectFullPhi[t] << '\n';
        }
        thetaWork.close();
        std::cout << "  Theta work map written: " << thetaWorkName
                  << std::endl;
    }

    std::vector<unsigned int> seeds = owenSeeds;
    if (seeds.empty())
        seeds.push_back(42u);
    const int seedCount = (int)seeds.size();
    const long long totalOrientationWork = (long long)nOrient * seedCount;
    if (m_mpiRank == 0)
    {
        std::cout << "Autofull final Owen averaging: " << seedCount
                  << " seed" << (seedCount == 1 ? "" : "s")
                  << " x " << nOrient << " orientations";
        if (seedCount > 1)
        {
            std::cout << " (";
            for (size_t i = 0; i < seeds.size(); ++i)
            {
                if (i) std::cout << ',';
                std::cout << seeds[i];
            }
            std::cout << ")";
        }
        std::cout << std::endl;
    }

    int myStart = m_mpiRank * nOrient / m_mpiSize;
    int myEnd = (m_mpiRank + 1) * nOrient / m_mpiSize;
    int myCount = myEnd - myStart;

    CalcTimer timer;
    timer.Start();
    if (m_mpiRank == 0)
        OutputStartTime(timer);

    std::vector<matrix> rows;
    std::vector<matrix> rowsNoShadow;
    rows.reserve(nZen + 1);
    rowsNoShadow.reserve(nZen + 1);
    for (int t = 0; t <= nZen; ++t)
    {
        matrix m(4, 4);
        m.Fill(0.0);
        rows.push_back(m);
        rowsNoShadow.push_back(m);
    }
    std::vector<std::vector<matrix>> seedRows;
    if (seedCount > 1)
    {
        seedRows.resize(seedCount);
        for (int s = 0; s < seedCount; ++s)
        {
            seedRows[s].reserve(nZen + 1);
            for (int t = 0; t <= nZen; ++t)
            {
                matrix m(4, 4);
                m.Fill(0.0);
                seedRows[s].push_back(m);
            }
        }
    }

    m_incomingEnergy = 0.0;
    handlerPO->m_outputEnergy = 0.0;
    handlerPO->m_extinctionCrossSectionOt = 0.0;
    handlerPO->m_hasExtinctionOt = false;
    const bool computeNoShadow = handlerPO->ComputeNoShadow();
    const double weight = 1.0 / ((double)nOrient * (double)seedCount);

    long long availMB = EffectiveMemAvailableMb();
    long long beamBudget = std::max(100LL, availMB / 2);
    int chunkSize = std::max(32, std::min(4096, std::min(myCount, (int)(beamBudget * 1024 / 350))));
    if (m_sobolChunkSize > 0)
        chunkSize = std::max(1, std::min(chunkSize, m_sobolChunkSize));
    int nChunks = (myCount + chunkSize - 1) / chunkSize;
    int nThreads = 1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        nThreads = omp_get_num_threads();
    }
#endif
    if (m_mpiRank == 0)
    {
        std::cerr << "Variable-phi memory: " << availMB << " MB available";
        if (HostMemoryBudgetOverrideKb() > 0)
            std::cerr << " (shared budget)";
        std::cerr << ", chunk="
                  << chunkSize << " orientations (" << nChunks << " chunks), "
                  << nThreads << " threads" << std::endl;
    }

    double phase1 = 0.0;
    double phase2 = 0.0;
    long long count = 0;
    std::vector<Beam> outBeams;
    m_scattering->PrepareForParallelTrace();

    double cosBetaSym = cos(betaSym);
    for (size_t seedIndex = 0; seedIndex < seeds.size(); ++seedIndex)
    {
        Sobol2D sobol(seeds[seedIndex]);
        std::vector<double> su, sv;
        sobol.generate(nOrient, su, sv);
        std::vector<std::pair<double,double>> orientations(nOrient);
        for (int i = 0; i < nOrient; ++i)
        {
            double beta = acos(1.0 - (1.0 - cosBetaSym) * su[i]);
            double gamma = gammaSym * sv[i];
            orientations[i] = {beta, gamma};
        }

        if (m_mpiRank == 0 && seedCount > 1)
            std::cout << "  Owen seed " << seeds[seedIndex]
                      << " (" << (seedIndex + 1) << "/" << seedCount << ")"
                      << std::endl;

        for (int chunk = 0; chunk < nChunks; ++chunk)
        {
            int iStart = myStart + chunk * chunkSize;
            int iEnd = std::min(iStart + chunkSize, myEnd);
            int thisChunk = iEnd - iStart;

            auto tp1 = std::chrono::high_resolution_clock::now();
            std::vector<PreparedOrientation> chunkPrepared(thisChunk);
            std::vector<double> chunkEnergies(thisChunk, 0.0);
            std::vector<double> chunkOutputEnergies(thisChunk, 0.0);
            ParallelExceptionState parallelError;

            #pragma omp parallel
            {
                Particle localParticle = *m_particle;
                Scattering *localScatter = m_scattering->CloneFor(&localParticle, &m_incidentLight);
                HandlerPO localHandler(&localParticle, &m_incidentLight,
                                       handlerPO->nTheta, m_scattering->m_wave);
                localHandler.ConfigureForThreadLocalPrepare(*handlerPO, localScatter);
                std::vector<Beam> localBeams;

                #pragma omp for schedule(dynamic, 4)
                for (int i = 0; i < thisChunk; ++i)
                {
                    if (parallelError.Failed())
                        continue;
                    try
                    {
                    int idx = iStart + i;
                    localParticle.Rotate(orientations[idx].first, orientations[idx].second, 0);
                    if (!shadowOff)
                        localScatter->FormShadowBeam(localBeams);
                    bool ok = localScatter->ScatterLight(0, 0, localBeams);
                    if (ok)
                    {
                        double beforeOutput = localHandler.m_outputEnergy;
                        localHandler.PrepareBeams(localBeams, weight, chunkPrepared[i]);
                        chunkOutputEnergies[i] = localHandler.m_outputEnergy - beforeOutput;
                    }
                    else
                    {
                        chunkPrepared[i].sinZenith = weight;
                    }
                    chunkEnergies[i] = localScatter->GetIncedentEnergy() * weight;
                    localBeams.clear();
                    }
                    catch (...)
                    {
                        parallelError.Capture();
                    }
                }

                delete localScatter;
            }
            parallelError.Rethrow();

            for (int i = 0; i < thisChunk; ++i)
            {
                m_incomingEnergy += chunkEnergies[i];
                handlerPO->m_outputEnergy += chunkOutputEnergies[i];
                handlerPO->m_extinctionCrossSectionOt +=
                    chunkPrepared[i].extinctionOt;
                handlerPO->m_hasExtinctionOt = true;
                ++count;
            }
            phase1 += std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - tp1).count();

            auto tp2 = std::chrono::high_resolution_clock::now();
            for (const auto &entry : rowsByWork)
            {
            const int groupOrientLimit = entry.first.nOrient;
            if (iStart >= groupOrientLimit)
                continue;
            const int activeCount =
                std::min(thisChunk, groupOrientLimit - iStart);
            if (activeCount <= 0)
                continue;

            const int nPhi = entry.first.nPhi;
            const int calcNphi = AutoFullVariableCalcPhi(
                nPhi, fftEps, minDirectPhi, canUseFftAverage,
                entry.first.directFullPhi);
            const std::vector<int> &rowIds = entry.second;
            const int rowCount = (int)rowIds.size();
            ScatteringRange rowSphere = outputSphere;
            rowSphere.nAzimuth = calcNphi;
            rowSphere.azinuthStep = M_2PI / calcNphi;
            rowSphere.isNonUniform = true;
            rowSphere.thetaValues.clear();
            rowSphere.thetaValues.reserve(rowCount);
            for (int row : rowIds)
                rowSphere.thetaValues.push_back(outputSphere.GetZenith(row));
            rowSphere.nZenith = rowCount - 1;
            rowSphere.zenithStart = rowSphere.thetaValues.front();
            rowSphere.zenithEnd = rowSphere.thetaValues.back();
            rowSphere.zenithStep = rowSphere.nZenith > 0
                ? (rowSphere.zenithEnd - rowSphere.zenithStart) / rowSphere.nZenith
                : 0.0;
            rowSphere.ComputeSphereDirections(m_incidentLight);

            Arr2D localM(calcNphi + 1, rowCount, 4, 4);
            localM.ClearArr();
            Arr2D localMns(calcNphi + 1, rowCount, 4, 4);
            localMns.ClearArr();

            ScatteringRange savedSphere = handlerPO->m_sphere;
            handlerPO->m_sphere = rowSphere;
            const double weightScale =
                (double)nOrient / (double)std::max(1, groupOrientLimit);
            double extraCutoff = 0.0;
            if (!rowIds.empty())
                extraCutoff = rowExtraImportanceCutoff[rowIds.front()];
            std::vector<PreparedOrientation> filteredPrepared;
            const std::vector<PreparedOrientation> *preparedForGroup =
                &chunkPrepared;
            int preparedStart = 0;
            int preparedCount = activeCount;
            if (extraCutoff > 0.0)
            {
                filteredPrepared = BuildPreparedSliceForThetaZone(
                    chunkPrepared, 0, activeCount, 1.0, extraCutoff);
                preparedForGroup = &filteredPrepared;
                preparedStart = 0;
                preparedCount = (int)filteredPrepared.size();
            }
            if (handlerPO->IsGpuEnabled())
            {
                for (int gpuStart = preparedStart;
                     gpuStart < preparedStart + preparedCount; )
                {
                    int gpuBatchSize = handlerPO->SelectGpuOrientationBatchSize(
                        *preparedForGroup, gpuStart,
                        preparedStart + preparedCount - gpuStart);
                    int gpuEnd = std::min(gpuStart + gpuBatchSize,
                                          preparedStart + preparedCount);
                    const bool useFftForGroup =
                        handlerPO->IsFftEnabled()
                        && !canUseFftAverage
                        && !entry.first.directFullPhi;
                    bool ok = useFftForGroup
                        ? handlerPO->HandleOrientationsToLocalGpuFftPhi(
                            *preparedForGroup, gpuStart, gpuEnd - gpuStart,
                            localM, localMns)
                        : handlerPO->HandleOrientationsToLocalGpu(
                            *preparedForGroup, gpuStart, gpuEnd - gpuStart,
                            localM, localMns);
                    if (!ok)
                    {
                        handlerPO->m_sphere = savedSphere;
                        throw std::runtime_error("variable-phi GPU diffraction failed");
                    }
                    gpuStart = gpuEnd;
                }
            }
            else
            {
                #pragma omp parallel
                {
                    Arr2D threadM(calcNphi + 1, rowCount, 4, 4);
                    threadM.ClearArr();
                    Arr2D threadMns(calcNphi + 1, rowCount, 4, 4);
                    threadMns.ClearArr();
                    std::vector<Arr2DC> localJ, localJns;
                    if (handlerPO->isCoh)
                    {
                        Arr2DC tmp(calcNphi + 1, rowCount, 2, 2);
                        tmp.ClearArr();
                        localJ.push_back(tmp);
                        Arr2DC tmpNs(calcNphi + 1, rowCount, 2, 2);
                        tmpNs.ClearArr();
                        localJns.push_back(tmpNs);
                    }

                    #pragma omp for schedule(dynamic, 1)
                    for (int i = 0; i < preparedCount; ++i)
                    {
                        const PreparedOrientation &prepared =
                            (*preparedForGroup)[preparedStart + i];
                        if (!prepared.beams.empty())
                            handlerPO->HandleBeamsToLocal(
                                prepared, threadM, localJ,
                                handlerPO->isCoh ? &localJns : nullptr);
                        if (handlerPO->isCoh && !localJ.empty())
                        {
                            double w = prepared.sinZenith;
                            HandlerPO::AddToMuellerLocal(localJ, w, threadM,
                                                         calcNphi, rowCount - 1);
                            HandlerPO::AddToMuellerLocal(localJns, w, threadMns,
                                                         calcNphi, rowCount - 1);
                            localJ[0].ClearArr();
                            localJns[0].ClearArr();
                        }
                    }

                    #pragma omp critical
                    {
                        for (int p = 0; p < calcNphi; ++p)
                            for (int t = 0; t < rowCount; ++t)
                            {
                                localM.insert(p, t, threadM(p, t));
                                if (computeNoShadow)
                                    localMns.insert(p, t, threadMns(p, t));
                            }
                    }
                }
            }
            handlerPO->m_sphere = savedSphere;

            if (m_mirrorGamma)
            {
                ApplyMirrorGammaMueller(localM, calcNphi, rowCount - 1);
                if (computeNoShadow)
                    ApplyMirrorGammaMueller(localMns, calcNphi, rowCount - 1);
            }

            for (int localRow = 0; localRow < rowCount; ++localRow)
            {
                int globalRow = rowIds[localRow];
                matrix averaged = AzimuthAverageOutputMatrix(localM, rowSphere,
                                                             calcNphi, localRow);
                if (std::fabs(weightScale - 1.0) > 1e-15)
                    averaged *= weightScale;
                rows[globalRow] += averaged;
                if (!seedRows.empty())
                    seedRows[seedIndex][globalRow] += averaged;
                if (computeNoShadow)
                {
                    matrix averagedNoShadow = AzimuthAverageOutputMatrix(
                        localMns, rowSphere, calcNphi, localRow);
                    if (std::fabs(weightScale - 1.0) > 1e-15)
                        averagedNoShadow *= weightScale;
                    rowsNoShadow[globalRow] += averagedNoShadow;
                }
            }
        }
        phase2 += std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - tp2).count();

            if (m_mpiRank == 0)
            {
                int progressIndex =
                    (int)((long long)seedIndex * nOrient + iEnd - 1);
                OutputProgress((int)totalOrientationWork, count, progressIndex,
                               chunk + 1, timer, -1);
            }
        }
    }

#ifdef USE_MPI
    {
        std::vector<double> send((nZen + 1) * 16, 0.0);
        std::vector<double> recv((nZen + 1) * 16, 0.0);
        for (int t = 0; t <= nZen; ++t)
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    send[(t * 16) + r * 4 + c] = rows[t][r][c];
        MPI_Reduce(send.data(), recv.data(), (int)send.size(),
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (m_mpiRank == 0)
            for (int t = 0; t <= nZen; ++t)
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c)
                        rows[t][r][c] = recv[(t * 16) + r * 4 + c];

        if (computeNoShadow)
        {
            std::fill(send.begin(), send.end(), 0.0);
            std::fill(recv.begin(), recv.end(), 0.0);
            for (int t = 0; t <= nZen; ++t)
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c)
                        send[(t * 16) + r * 4 + c] = rowsNoShadow[t][r][c];
            MPI_Reduce(send.data(), recv.data(), (int)send.size(),
                       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (m_mpiRank == 0)
                for (int t = 0; t <= nZen; ++t)
                    for (int r = 0; r < 4; ++r)
                        for (int c = 0; c < 4; ++c)
                            rowsNoShadow[t][r][c] = recv[(t * 16) + r * 4 + c];
        }

        if (!seedRows.empty())
        {
            const size_t rowBlock = (size_t)(nZen + 1) * 16;
            std::vector<double> sendSeeds(seedRows.size() * rowBlock, 0.0);
            std::vector<double> recvSeeds(seedRows.size() * rowBlock, 0.0);
            for (size_t s = 0; s < seedRows.size(); ++s)
                for (int t = 0; t <= nZen; ++t)
                    for (int r = 0; r < 4; ++r)
                        for (int c = 0; c < 4; ++c)
                        {
                            size_t idx = s * rowBlock + (size_t)t * 16
                                       + (size_t)r * 4 + c;
                            sendSeeds[idx] = seedRows[s][t][r][c];
                        }
            MPI_Reduce(sendSeeds.data(), recvSeeds.data(),
                       (int)sendSeeds.size(), MPI_DOUBLE, MPI_SUM, 0,
                       MPI_COMM_WORLD);
            if (m_mpiRank == 0)
                for (size_t s = 0; s < seedRows.size(); ++s)
                    for (int t = 0; t <= nZen; ++t)
                        for (int r = 0; r < 4; ++r)
                            for (int c = 0; c < 4; ++c)
                            {
                                size_t idx = s * rowBlock + (size_t)t * 16
                                           + (size_t)r * 4 + c;
                                seedRows[s][t][r][c] = recvSeeds[idx];
                            }
        }

        double totalEnergy = 0.0;
        MPI_Reduce(&m_incomingEnergy, &totalEnergy, 1, MPI_DOUBLE,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        if (m_mpiRank == 0)
            m_incomingEnergy = totalEnergy;
        double totalOutput = 0.0;
        MPI_Reduce(&handlerPO->m_outputEnergy, &totalOutput, 1, MPI_DOUBLE,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        if (m_mpiRank == 0)
            handlerPO->m_outputEnergy = totalOutput;
        double totalExtOt = 0.0;
        MPI_Reduce(&handlerPO->m_extinctionCrossSectionOt, &totalExtOt, 1,
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (m_mpiRank == 0)
        {
            handlerPO->m_extinctionCrossSectionOt = totalExtOt;
            handlerPO->m_hasExtinctionOt = true;
        }
    }
#endif

    EraseConsoleLine(60);
    handlerPO->m_sphere = outputSphere;
    if (m_mpiRank == 0)
    {
        std::cout << "Variable-phi Phase 1 (tracing): " << std::fixed
                  << std::setprecision(2) << phase1 << " s" << std::endl;
        const char *diffractionMode = "OpenMP";
        if (handlerPO->IsGpuEnabled())
        {
            if (canUseFftAverage)
            {
                diffractionMode =
                    (directFullPhiRows == nZen + 1)
                        ? "CUDA direct full-phi"
                        : (directFullPhiRows > 0
                            ? "CUDA low-phi average + direct full-phi"
                            : "CUDA low-phi average");
            }
            else
            {
                diffractionMode = handlerPO->IsFftEnabled()
                    ? "CUDA FFT"
                    : "CUDA";
            }
        }
        std::cout << "Variable-phi Phase 2 (diffraction, "
                  << diffractionMode << "): " << phase2 << " s"
                  << std::endl;
        std::cout << "Variable-phi total: " << phase1 + phase2 << " s" << std::endl;

        const bool projectSymmetryInError =
            AutoFullSymmetryErrorProjectionEnabled(betaSym, gammaSym);
        if (AutoFullRobustOwenOutputEnabled() && seedRows.size() >= 3)
        {
            ReplaceRowsByRobustOwenMedian(rows, seedRows);
            std::cout << "Autofull robust Owen output: per-element median"
                      << " over " << seedRows.size()
                      << " final seeds"
                      << " (MBS_AUTOFULL_ROBUST_OWEN_OUTPUT=0 disables)"
                      << std::endl;
        }
        if (projectSymmetryInError && !seedRows.empty())
            ProjectSeedRowsRandomOrientationMuellerSymmetry(seedRows);

        if (AutoFullSymmetryOutputProjectionEnabled(betaSym, gammaSym))
        {
            ProjectRowsRandomOrientationMuellerSymmetry(rows);
            if (computeNoShadow)
                ProjectRowsRandomOrientationMuellerSymmetry(rowsNoShadow);
            std::cout << "Autofull output symmetry projection: zeroed"
                      << " M13/M14/M23/M24/M31/M32/M41/M42"
                      << " (MBS_AUTOFULL_SYMMETRY_OUTPUT_PROJECT=0 disables)"
                      << std::endl;
        }
        else if (projectSymmetryInError)
        {
            std::cout << "Legacy symmetry projection explicitly applied to convergence"
                      << " and Owen diagnostics only; output Mueller is not"
                      << " cosmetically zeroed"
                      << std::endl;
        }

        WriteAveragedRowsFile(m_resultDirName, outputSphere, rows,
                              m_incomingEnergy,
                              handlerPO->m_extinctionCrossSectionOt,
                              handlerPO->m_hasExtinctionOt,
                              handlerPO->HasAbsorptionAccounting(),
                              handlerPO->m_integralSummary);
        double finalMeanError = 0.0;
        if (!seedRows.empty())
        {
            std::string spreadName =
                m_resultDirName + "_autofull_owen_spread_by_theta.dat";
            finalMeanError =
                WriteOwenSpreadByThetaFile(spreadName, outputSphere, seedRows);
        }
        if (computeNoShadow)
        {
            std::string nsName = m_resultDirName + "_noshadow";
            std::string dummySummary;
            WriteAveragedRowsFile(nsName, outputSphere, rowsNoShadow,
                                  m_incomingEnergy,
                                  handlerPO->m_extinctionCrossSectionOt,
                                  handlerPO->m_hasExtinctionOt,
                                  false, dummySummary);
        }
        OutputStatisticsPO(timer, totalOrientationWork, m_resultDirName);
        return finalMeanError;
    }
    return 0.0;
}

int TracerPOTotal::TraceAdaptivePhi(double eps, int nOrient,
                                    double betaSym, double gammaSym)
{
    HandlerPO *handler = dynamic_cast<HandlerPO*>(m_handler);
    if (!handler)
        throw std::runtime_error("adaptive phi requires the PO handler");
    if (!(eps > 0.0 && eps < 1.0))
        throw std::invalid_argument("adaptive phi tolerance must be in (0, 1)");
    m_convergenceTarget = eps;

    const int minPhi = std::max(6,
        ((m_adaptiveLimits.minPhiPoints + 5) / 6) * 6);
    const int maxPhi = std::max(minPhi,
        (m_adaptiveLimits.maxPhiPoints / 6) * 6);
    std::vector<int> candidates;
    for (int value = minPhi; value < maxPhi; )
    {
        candidates.push_back(value);
        long long grown = (long long)std::ceil(value * 1.5);
        long long rounded = ((grown + 5) / 6) * 6;
        int next = rounded > maxPhi ? maxPhi : (int)rounded;
        if (next <= value)
            next = value + 6;
        value = next;
    }
    if (candidates.empty() || candidates.back() != maxPhi)
        candidates.push_back(maxPhi);

    std::vector<WeightedOrientation> orientations = BuildSobolOrientations(
        std::max(16, nOrient), betaSym, gammaSym, 42u);
    std::vector<PreparedOrientation> prepared;
    double incomingEnergy = 0.0;
    double outputEnergy = 0.0;
    PrepareOrientationProbe(orientations, prepared,
                            incomingEnergy, outputEnergy);

    const ScatteringRange baseSphere = handler->m_sphere;
    PreparedProbeResult previous;
    int previousPhi = 0;
    int selectedPhi = 0;
    int passStreak = 0;
    const int stablePasses = std::max(1, m_adaptiveLimits.stablePasses);

    std::cout << "Adaptive phi: target=" << eps * 100.0
              << "%, pilot orientations=" << orientations.size()
              << ", range=" << minPhi << ".." << maxPhi << std::endl;
    std::cout << "  N_phi represents the full laboratory-alpha average; "
              << "particle axial symmetry is already applied to gamma="
              << RadToDeg(gammaSym) << " deg and does not generally shorten alpha"
              << std::endl;
    for (int phi : candidates)
    {
        ScatteringRange sphere = baseSphere;
        sphere.nAzimuth = phi;
        sphere.azinuthStep = M_2PI / phi;
        PreparedProbeResult current = EvaluatePreparedProbe(
            handler, prepared, sphere, m_incidentLight, m_mirrorGamma);

        double error = 1.0;
        double worstTheta = 0.0;
        int worstElement = 0;
        if (!previous.rows.empty())
            error = ProbeRowsRelativeError(
                current.rows, previous.rows, sphere,
                current.scatteringIntegral, previous.scatteringIntegral,
                &worstTheta, &worstElement);

        const bool pass = !previous.rows.empty() && error <= eps;
        passStreak = pass ? passStreak + 1 : 0;
        const std::string status = previous.rows.empty()
            ? "baseline" : (pass ? "pass" : "refine");
        RecordConvergenceStep("phi", 0, phi, (int)orientations.size(),
                              error, worstTheta, worstElement, status);
        std::cout << "  N_phi=" << phi << " delta=" << std::fixed
                  << std::setprecision(4) << error * 100.0 << "%"
                  << " worst_theta=" << std::setprecision(3)
                  << worstTheta << " deg, element=" << worstElement
                  << (current.usedGpu ? " GPU" : " CPU")
                  << (pass ? " pass" : "") << std::endl;

        if (passStreak >= stablePasses)
        {
            selectedPhi = previousPhi;
            break;
        }
        previous = current;
        previousPhi = phi;
    }

    if (selectedPhi <= 0)
    {
        RecordConvergenceStep("phi", 0, maxPhi,
                              (int)orientations.size(), 1.0, 0.0, -1,
                              "limit");
        std::ostringstream message;
        message << "N_phi did not converge to " << eps
                << " before --max-phi-points " << maxPhi
                << "; increase the limit or loosen the tolerance";
        throw std::runtime_error(message.str());
    }

    ScatteringRange selectedSphere = baseSphere;
    selectedSphere.nAzimuth = selectedPhi;
    selectedSphere.azinuthStep = M_2PI / selectedPhi;
    selectedSphere.ComputeSphereDirections(m_incidentLight);
    handler->SetScatteringSphere(selectedSphere);
    if (handler->IsFftEnabled()
        && handler->FftPhiFactor() == 0
        && !HandlerPO::HasNumericFftPhiFactorOverride())
        handler->AutoSelectFftPhiFactor(eps);
    std::cout << "Adaptive phi selected N_phi=" << selectedPhi << std::endl;
    RecordConvergenceStep("phi", 0, selectedPhi,
                          (int)orientations.size(), 0.0, 0.0, -1,
                          "selected");
    return selectedPhi;
}

int TracerPOTotal::TraceAdaptiveReflections(double eps, int nOrient,
                                            double betaSym,
                                            double gammaSym,
                                            int startReflection)
{
    HandlerPO *handler = dynamic_cast<HandlerPO*>(m_handler);
    if (!handler)
        throw std::runtime_error("adaptive reflections require the PO handler");
    if (!(eps > 0.0 && eps < 1.0))
        throw std::invalid_argument(
            "adaptive reflection tolerance must be in (0, 1)");
    m_convergenceTarget = eps;

    const int maxReflection = std::max(1,
        m_adaptiveLimits.maxReflections);
    const int requestedStart = startReflection > 0
        ? startReflection : m_adaptiveLimits.minReflections;
    const int first = std::max(
        m_adaptiveLimits.minReflections,
        std::max(1, std::min(requestedStart, maxReflection)));
    const int stablePasses = std::max(1, m_adaptiveLimits.stablePasses);
    const std::vector<WeightedOrientation> orientations =
        BuildSobolOrientations(std::max(16, nOrient), betaSym, gammaSym, 42u);

    PreparedProbeResult previous;
    double previousOutput = 0.0;
    int previousN = 0;
    int selectedN = 0;
    int passStreak = 0;
    std::cout << "Adaptive reflections: target=" << eps * 100.0
              << "%, pilot orientations=" << orientations.size()
              << ", range=" << first << ".." << maxReflection << std::endl;

    for (int n = first; n <= maxReflection; ++n)
    {
        m_scattering->SetMaxReflections(n);
        std::vector<PreparedOrientation> prepared;
        double incomingEnergy = 0.0;
        double outputEnergy = 0.0;
        PrepareOrientationProbe(orientations, prepared,
                                incomingEnergy, outputEnergy);
        PreparedProbeResult current = EvaluatePreparedProbe(
            handler, prepared, handler->m_sphere,
            m_incidentLight, m_mirrorGamma);

        double error = 1.0;
        double worstTheta = 0.0;
        int worstElement = 0;
        if (!previous.rows.empty())
        {
            error = ProbeRowsRelativeError(
                current.rows, previous.rows, handler->m_sphere,
                current.scatteringIntegral, previous.scatteringIntegral,
                &worstTheta, &worstElement);
            const double energyDen = std::max(
                std::max(std::fabs(outputEnergy),
                         std::fabs(previousOutput)), 1e-30);
            const double energyError = std::fabs(
                outputEnergy - previousOutput) / energyDen;
            if (std::isfinite(energyError) && energyError > error)
            {
                error = energyError;
                worstElement = 17; // Outgoing-energy criterion.
            }
        }

        const bool pass = !previous.rows.empty() && error <= eps;
        passStreak = pass ? passStreak + 1 : 0;
        const std::string status = previous.rows.empty()
            ? "baseline" : (pass ? "pass" : "refine");
        RecordConvergenceStep("reflections", 0, n,
                              (int)orientations.size(), error,
                              worstTheta, worstElement, status);
        std::cout << "  n=" << n << " delta=" << std::fixed
                  << std::setprecision(4) << error * 100.0 << "%"
                  << " worst_theta=" << std::setprecision(3)
                  << worstTheta << " deg, element=" << worstElement
                  << (current.usedGpu ? " GPU" : " CPU")
                  << (pass ? " pass" : "") << std::endl;

        if (passStreak >= stablePasses)
        {
            selectedN = previousN;
            break;
        }
        previous = current;
        previousOutput = outputEnergy;
        previousN = n;
    }

    if (selectedN <= 0)
    {
        m_scattering->SetMaxReflections(maxReflection);
        RecordConvergenceStep("reflections", 0, maxReflection,
                              (int)orientations.size(), 1.0, 0.0, -1,
                              "limit");
        std::ostringstream message;
        message << "reflection depth did not converge to " << eps
                << " before --max-reflections " << maxReflection
                << "; increase the limit or loosen the tolerance";
        throw std::runtime_error(message.str());
    }

    m_scattering->SetMaxReflections(selectedN);
    RecordConvergenceStep("reflections", 0, selectedN,
                          (int)orientations.size(), 0.0, 0.0, -1,
                          "selected");
    std::cout << "Adaptive reflections selected n=" << selectedN
              << std::endl;
    return selectedN;
}

void TracerPOTotal::TraceAdaptiveTheta(int nOrient, double betaSym,
                                       double gammaSym, double eps,
                                       int maxDepth, bool gridOnly)
{
    HandlerPO *handler = dynamic_cast<HandlerPO*>(m_handler);
    if (!handler)
        throw std::runtime_error("adaptive theta requires the PO handler");
    if (!(eps > 0.0 && eps < 1.0))
        throw std::invalid_argument(
            "adaptive theta tolerance must be in (0, 1)");
    m_convergenceTarget = eps;

    const int minPoints = std::max(17, m_adaptiveLimits.minThetaPoints);
    const int maxPoints = std::max(minPoints, m_adaptiveLimits.maxThetaPoints);
    const double diameter = std::max(
        m_particle->MaximalDimention(), 1e-12);
    const double physicsStep = std::max(
        m_scattering->m_wave / (12.0 * diameter),
        M_PI / (2.0 * maxPoints));
    const int baseIntervals = std::min(minPoints - 1, maxPoints - 1);

    std::vector<WeightedOrientation> orientations = BuildSobolOrientations(
        std::max(16, nOrient), betaSym, gammaSym, 42u);
    std::vector<PreparedOrientation> prepared;
    double incomingEnergy = 0.0;
    double outputEnergy = 0.0;
    PrepareOrientationProbe(orientations, prepared,
                            incomingEnergy, outputEnergy);

    const ScatteringRange savedSphere = handler->m_sphere;
    std::map<double, MuellerRowAverage> values;
    std::vector<double> seed;
    seed.reserve(baseIntervals + 32);
    for (int i = 0; i <= baseIntervals; ++i)
        seed.push_back(M_PI * i / baseIntervals);

    // Generic pole refinement follows the diffraction scale only. No
    // particle-specific halo or rainbow angles are embedded in this seed.
    for (double angle = physicsStep;
         angle < M_PI / baseIntervals; angle *= 2.0)
    {
        seed.push_back(angle);
        seed.push_back(M_PI - angle);
    }
    std::sort(seed.begin(), seed.end());
    seed.erase(std::unique(seed.begin(), seed.end(),
        [](double a, double b) { return std::fabs(a - b) < 1e-13; }),
        seed.end());
    if ((int)seed.size() > maxPoints)
        throw std::runtime_error(
            "--max-theta-points is too small for the initial adaptive grid");

    auto evaluate = [&](const std::vector<double> &theta)
    {
        ScatteringRange sphere = savedSphere;
        sphere.isNonUniform = true;
        sphere.thetaValues = theta;
        sphere.nZenith = (int)theta.size() - 1;
        sphere.zenithStart = theta.front();
        sphere.zenithEnd = theta.back();
        sphere.zenithStep = sphere.nZenith > 0
            ? (sphere.zenithEnd - sphere.zenithStart) / sphere.nZenith
            : 0.0;
        return EvaluatePreparedProbe(handler, prepared, sphere,
                                     m_incidentLight, m_mirrorGamma);
    };

    PreparedProbeResult initial = evaluate(seed);
    for (size_t i = 0; i < seed.size(); ++i)
        values[seed[i]] = initial.rows[i];

    const double baseStep = M_PI / baseIntervals;
    int requiredDepth = 1;
    if (physicsStep < baseStep)
        requiredDepth += (int)std::ceil(
            std::log(baseStep / physicsStep) / std::log(2.0));
    const int depthLimit = std::max(maxDepth,
        std::min(20, requiredDepth + 2));
    bool converged = false;
    double finalMaxError = 1.0;

    std::cout << "Adaptive theta: target=" << eps * 100.0
              << "%, pilot orientations=" << orientations.size()
              << ", initial points=" << values.size()
              << ", configured range=" << minPoints << ".." << maxPoints
              << ", minimum step=" << RadToDeg(physicsStep)
              << " deg" << std::endl;

    for (int depth = 0; depth < depthLimit; ++depth)
    {
        struct IntervalProbe
        {
            double left;
            double right;
            double first;
            double second;
            MuellerRowAverage leftValue;
            MuellerRowAverage rightValue;
        };
        std::vector<IntervalProbe> intervals;
        std::vector<double> theta;
        for (auto right = std::next(values.begin());
             right != values.end(); ++right)
        {
            auto left = std::prev(right);
            const double width = right->first - left->first;
            if (width <= physicsStep)
                continue;
            IntervalProbe probe;
            probe.left = left->first;
            probe.right = right->first;
            probe.first = probe.left + width / 3.0;
            probe.second = probe.left + 2.0 * width / 3.0;
            probe.leftValue = left->second;
            probe.rightValue = right->second;
            intervals.push_back(probe);
            theta.push_back(probe.first);
            theta.push_back(probe.second);
        }

        if (intervals.empty())
        {
            converged = true;
            finalMaxError = 0.0;
            break;
        }

        PreparedProbeResult actual = evaluate(theta);
        std::vector<std::pair<double, MuellerRowAverage>> additions;
        finalMaxError = 0.0;
        double worstTheta = 0.0;
        int worstElement = 0;
        for (size_t i = 0; i < intervals.size(); ++i)
        {
            double intervalError = 0.0;
            const double fractions[2] = { 1.0 / 3.0, 2.0 / 3.0 };
            const double angles[2] = {
                intervals[i].first, intervals[i].second
            };
            for (int sample = 0; sample < 2; ++sample)
            {
                MuellerRowAverage interpolated;
                for (int k = 0; k < 16; ++k)
                    interpolated.v[k] =
                        intervals[i].leftValue.v[k] * (1.0 - fractions[sample])
                        + intervals[i].rightValue.v[k] * fractions[sample];
                int element = 0;
                const double error = MuellerRowRelativeErrorValues(
                    actual.rows[2 * i + sample].v,
                    interpolated.v, &element);
                if (std::isfinite(error) && error > intervalError)
                    intervalError = error;
                if (std::isfinite(error) && error > finalMaxError)
                {
                    finalMaxError = error;
                    worstTheta = RadToDeg(angles[sample]);
                    worstElement = element;
                }
            }
            if (intervalError > eps)
            {
                additions.push_back(std::make_pair(
                    intervals[i].first, actual.rows[2 * i]));
                additions.push_back(std::make_pair(
                    intervals[i].second, actual.rows[2 * i + 1]));
            }
        }

        std::map<double, MuellerRowAverage> fullyRefined = values;
        for (size_t i = 0; i < intervals.size(); ++i)
        {
            fullyRefined[intervals[i].first] = actual.rows[2 * i];
            fullyRefined[intervals[i].second] = actual.rows[2 * i + 1];
        }
        const double coarseIntegral = AdaptiveThetaIntegral(values);
        const double refinedIntegral = AdaptiveThetaIntegral(fullyRefined);
        const double integralDen = std::max(
            std::max(std::fabs(coarseIntegral),
                     std::fabs(refinedIntegral)), 1e-30);
        const double integralError = std::fabs(
            refinedIntegral - coarseIntegral) / integralDen;
        if (std::isfinite(integralError) && integralError > finalMaxError)
        {
            finalMaxError = integralError;
            worstTheta = 0.0;
            worstElement = 16;
        }
        if (std::isfinite(integralError) && integralError > eps)
        {
            additions.clear();
            for (size_t i = 0; i < intervals.size(); ++i)
            {
                additions.push_back(std::make_pair(
                    intervals[i].first, actual.rows[2 * i]));
                additions.push_back(std::make_pair(
                    intervals[i].second, actual.rows[2 * i + 1]));
            }
        }

        const int available = maxPoints - (int)values.size();
        if ((int)additions.size() > available)
        {
            RecordConvergenceStep("theta", depth, (int)values.size(),
                                  (int)additions.size(), finalMaxError,
                                  worstTheta, worstElement, "limit");
            std::ostringstream message;
            message << "theta grid did not converge to " << eps
                    << " before --max-theta-points " << maxPoints
                    << "; increase the limit or loosen the tolerance";
            throw std::runtime_error(message.str());
        }
        for (const auto &addition : additions)
            values[addition.first] = addition.second;

        const std::string status = additions.empty() ? "pass" : "refine";
        RecordConvergenceStep("theta", depth, (int)values.size(),
                              (int)intervals.size(), finalMaxError,
                              worstTheta, worstElement, status);
        std::cout << "  theta pass " << depth + 1
                  << ": points=" << values.size()
                  << ", added=" << additions.size()
                  << ", max interpolation error=" << std::fixed
                  << std::setprecision(4) << finalMaxError * 100.0
                  << "% at " << std::setprecision(3) << worstTheta
                  << " deg, element=" << worstElement << std::endl;

        if (additions.empty())
        {
            converged = true;
            break;
        }
    }

    if (!converged)
    {
        RecordConvergenceStep("theta", depthLimit, (int)values.size(), 0,
                              finalMaxError, 0.0, -1, "depth_limit");
        std::ostringstream message;
        message << "theta grid still changed after " << depthLimit
                << " refinements; increase --max-theta-points or loosen"
                << " --auto-theta-grid";
        throw std::runtime_error(message.str());
    }

    ScatteringRange finalSphere = savedSphere;
    finalSphere.isNonUniform = true;
    finalSphere.thetaValues.clear();
    finalSphere.thetaValues.reserve(values.size());
    for (const auto &entry : values)
        finalSphere.thetaValues.push_back(entry.first);
    finalSphere.nZenith = (int)finalSphere.thetaValues.size() - 1;
    finalSphere.zenithStart = finalSphere.thetaValues.front();
    finalSphere.zenithEnd = finalSphere.thetaValues.back();
    finalSphere.zenithStep = finalSphere.nZenith > 0
        ? (finalSphere.zenithEnd - finalSphere.zenithStart)
            / finalSphere.nZenith
        : 0.0;
    finalSphere.ComputeSphereDirections(m_incidentLight);
    handler->SetScatteringSphere(finalSphere);
    RecordConvergenceStep("theta", depthLimit,
                          (int)finalSphere.thetaValues.size(),
                          (int)orientations.size(), finalMaxError,
                          0.0, -1, "selected");

    std::cout << "Adaptive theta selected "
              << finalSphere.thetaValues.size() << " points" << std::endl;
    if (!gridOnly)
        TraceFromSobol(nOrient, betaSym, gammaSym);
}

void TracerPOTotal::TraceAdaptiveEuler(double eps, double betaSym,
                                       double gammaSym,
                                       int maxOrientOverride,
                                       int &nBeta, int &nGamma,
                                       bool runFinalDiffraction)
{
    HandlerPO *handler = dynamic_cast<HandlerPO*>(m_handler);
    if (!handler)
        throw std::runtime_error("adaptive Euler requires the PO handler");
    if (!(eps > 0.0 && eps < 1.0))
        throw std::invalid_argument(
            "adaptive Euler tolerance must be in (0, 1)");
    m_convergenceTarget = eps;

    AutoFullOldautoOrientEstimate startEstimate = AutoFullOldautoEstimate(
        betaSym, gammaSym, m_scattering->m_wave,
        m_particle->MaximalDimention(), m_ringPoints, 16);
    AutoFullOldautoOrientEstimate capEstimate = AutoFullOldautoEstimate(
        betaSym, gammaSym, m_scattering->m_wave,
        m_particle->MaximalDimention(), m_ringPoints, 1);

    const int betaMin = std::max(2, m_adaptiveLimits.minBetaPoints);
    const int gammaMin = std::max(6,
        ((m_adaptiveLimits.minGammaPoints + 5) / 6) * 6);
    const int betaCap = std::max(betaMin, m_adaptiveLimits.maxBetaPoints);
    const int gammaCap = std::max(gammaMin,
        (m_adaptiveLimits.maxGammaPoints / 6) * 6);
    const long long configuredProduct = (long long)betaCap * gammaCap;
    const long long physicsTotal = std::max<long long>(
        64, capEstimate.total);
    const long long defaultTotalCap = std::min(
        configuredProduct, std::max(4096LL, physicsTotal));
    const long long totalCap = maxOrientOverride > 0
        ? maxOrientOverride : defaultTotalCap;

    nBeta = std::max(betaMin, std::min(betaCap,
        std::max(betaMin, startEstimate.nBeta)));
    nGamma = std::max(gammaMin, std::min(gammaCap,
        ((std::max(gammaMin, startEstimate.nGamma) + 5) / 6) * 6));
    while ((long long)nBeta * nGamma > totalCap && nGamma > gammaMin)
        nGamma = std::max(gammaMin, (nGamma / 2 / 6) * 6);
    while ((long long)nBeta * nGamma > totalCap && nBeta > betaMin)
        nBeta = std::max(betaMin, nBeta / 2);
    if ((long long)nBeta * nGamma > totalCap)
        throw std::runtime_error(
            "adaptive Euler minimum beta/gamma grid exceeds the orientation cap");

    typedef std::pair<int, int> EulerKey;
    std::map<EulerKey, PreparedProbeResult> cache;
    auto evaluate = [&](int betaCount, int gammaCount)
        -> const PreparedProbeResult &
    {
        const EulerKey key(betaCount, gammaCount);
        std::map<EulerKey, PreparedProbeResult>::iterator found =
            cache.find(key);
        if (found != cache.end())
            return found->second;
        if ((long long)betaCount * gammaCount > totalCap)
        {
            std::ostringstream message;
            message << "adaptive Euler needs more than " << totalCap
                    << " orientations; increase --max-orientations";
            throw std::runtime_error(message.str());
        }
        std::vector<WeightedOrientation> orientations =
            BuildEulerOrientations(betaCount, gammaCount,
                                   betaSym, gammaSym);
        std::vector<PreparedOrientation> prepared;
        double incomingEnergy = 0.0;
        double outputEnergy = 0.0;
        PrepareOrientationProbe(orientations, prepared,
                                incomingEnergy, outputEnergy);
        PreparedProbeResult result = EvaluatePreparedProbe(
            handler, prepared, handler->m_sphere,
            m_incidentLight, m_mirrorGamma);
        return cache.insert(std::make_pair(key, result)).first->second;
    };

    auto nextBeta = [&](int value)
    {
        return std::min(betaCap, value > betaCap / 2
            ? betaCap : value * 2);
    };
    auto nextGamma = [&](int value)
    {
        int next = value > gammaCap / 2 ? gammaCap : value * 2;
        next = std::min(gammaCap, ((next + 5) / 6) * 6);
        return next;
    };

    const int stablePasses = std::max(1, m_adaptiveLimits.stablePasses);
    const double betaEps = ResolveAdaptiveTolerance(
        m_adaptiveLimits.betaTolerance, eps);
    const double gammaEps = ResolveAdaptiveTolerance(
        m_adaptiveLimits.gammaTolerance, eps);
    const double jointEps = std::min(betaEps, gammaEps);
    std::cout << "Adaptive Euler: beta target=" << betaEps * 100.0
              << "%, gamma target=" << gammaEps * 100.0
              << "%, joint target=" << jointEps * 100.0
              << "%, start=" << nBeta << " x " << nGamma
              << ", ranges=" << betaMin << ".." << betaCap
              << " x " << gammaMin << ".." << gammaCap
              << ", total cap=" << totalCap << std::endl;

    auto convergeAxis = [&](bool betaAxis, int sweep)
    {
        int passStreak = 0;
        const PreparedProbeResult *previous = &evaluate(nBeta, nGamma);
        int previousPrimary = betaAxis ? nBeta : nGamma;
        while (passStreak < stablePasses)
        {
            const int candidateBeta = betaAxis ? nextBeta(nBeta) : nBeta;
            const int candidateGamma = betaAxis ? nGamma : nextGamma(nGamma);
            if (candidateBeta == nBeta && candidateGamma == nGamma)
            {
                std::ostringstream message;
                message << (betaAxis ? "N_beta" : "N_gamma")
                        << " did not converge before "
                        << (betaAxis ? "--max-beta-points "
                                     : "--max-gamma-points ")
                        << (betaAxis ? betaCap : gammaCap);
                throw std::runtime_error(message.str());
            }
            const PreparedProbeResult &current = evaluate(
                candidateBeta, candidateGamma);
            double worstTheta = 0.0;
            int worstElement = 0;
            const double error = ProbeRowsRelativeError(
                current.rows, previous->rows, handler->m_sphere,
                current.scatteringIntegral, previous->scatteringIntegral,
                &worstTheta, &worstElement);
            const double axisEps = betaAxis ? betaEps : gammaEps;
            const bool pass = error <= axisEps;
            passStreak = pass ? passStreak + 1 : 0;
            m_convergenceTarget = axisEps;
            RecordConvergenceStep(betaAxis ? "beta" : "gamma", sweep,
                                  betaAxis ? candidateBeta : candidateGamma,
                                  betaAxis ? candidateGamma : candidateBeta,
                                  error, worstTheta, worstElement,
                                  pass ? "pass" : "refine");
            std::cout << "  " << (betaAxis ? "N_beta" : "N_gamma")
                      << " probe " << candidateBeta << " x "
                      << candidateGamma << ": delta=" << std::fixed
                      << std::setprecision(4) << error * 100.0 << "%"
                      << (pass ? " pass" : "") << std::endl;

            if (passStreak >= stablePasses)
            {
                if (betaAxis)
                    nBeta = previousPrimary;
                else
                    nGamma = previousPrimary;
                break;
            }
            nBeta = candidateBeta;
            nGamma = candidateGamma;
            previousPrimary = betaAxis ? nBeta : nGamma;
            previous = &current;
        }
    };

    bool jointConverged = false;
    const int maxEulerSweeps = std::max(1, m_adaptiveLimits.maxEulerSweeps);
    for (int sweep = 1; sweep <= maxEulerSweeps; ++sweep)
    {
        convergeAxis(true, sweep);
        convergeAxis(false, sweep);

        const PreparedProbeResult *previous = &evaluate(nBeta, nGamma);
        int jointBeta = nBeta;
        int jointGamma = nGamma;
        int jointPasses = 0;
        while (jointPasses < stablePasses)
        {
            const int candidateBeta = nextBeta(jointBeta);
            const int candidateGamma = nextGamma(jointGamma);
            if (candidateBeta == jointBeta && candidateGamma == jointGamma)
                break;
            const PreparedProbeResult &current = evaluate(
                candidateBeta, candidateGamma);
            double worstTheta = 0.0;
            int worstElement = 0;
            const double error = ProbeRowsRelativeError(
                current.rows, previous->rows, handler->m_sphere,
                current.scatteringIntegral, previous->scatteringIntegral,
                &worstTheta, &worstElement);
            const bool pass = error <= jointEps;
            jointPasses = pass ? jointPasses + 1 : 0;
            m_convergenceTarget = jointEps;
            RecordConvergenceStep("beta_gamma_joint", sweep,
                                  candidateBeta, candidateGamma, error,
                                  worstTheta, worstElement,
                                  pass ? "pass" : "refine");
            std::cout << "  joint beta/gamma probe " << candidateBeta
                      << " x " << candidateGamma << ": delta="
                      << std::fixed << std::setprecision(4)
                      << error * 100.0 << "%"
                      << (pass ? " pass" : "") << std::endl;
            if (pass && jointPasses >= stablePasses)
            {
                nBeta = jointBeta;
                nGamma = jointGamma;
                jointConverged = true;
                break;
            }
            if (!pass)
            {
                nBeta = candidateBeta;
                nGamma = candidateGamma;
                break;
            }
            jointBeta = candidateBeta;
            jointGamma = candidateGamma;
            previous = &current;
        }
        if (jointConverged)
            break;
    }

    if (!jointConverged)
    {
        std::ostringstream message;
        message << "N_beta/N_gamma did not pass joint convergence after "
                << maxEulerSweeps
                << " sweeps; increase controller.max_euler_sweeps or the limits";
        throw std::runtime_error(message.str());
    }

    m_lastOldAutoFullGridValid = true;
    m_lastOldAutoFullN = m_scattering->GetMaxReflections();
    m_lastOldAutoFullNphi = handler->m_sphere.nAzimuth;
    m_lastOldAutoFullNBeta = nBeta;
    m_lastOldAutoFullNGamma = nGamma;
    m_lastOldAutoFullDiv = 0;
    m_lastOldAutoFullBetaSym = betaSym;
    m_lastOldAutoFullGammaSym = gammaSym;
    RecordConvergenceStep("beta_gamma_joint", maxEulerSweeps + 1,
                          nBeta, nGamma,
                          0.0, 0.0, -1, "selected");
    std::cout << "Adaptive Euler selected " << nBeta << " x "
              << nGamma << " = " << (long long)nBeta * nGamma
              << " orientations" << std::endl;

    if (runFinalDiffraction)
        TraceFromEulerQuadrature(nBeta, nGamma, betaSym, gammaSym);
}

void TracerPOTotal::EvaluateOwenControlBatch(
    int orientStart, int orientCount, int sobolCount, unsigned int seed,
    double sampleWeight, double betaSym, double gammaSym,
    const ScatteringRange &ctrlSphere, const std::vector<int> &ctrlIdx,
    int nAz, std::vector<double> &ctrlAvg, double &energyAvg)
{
    HandlerPO *hp = dynamic_cast<HandlerPO*>(m_handler);
    if (!hp)
        throw std::runtime_error("EvaluateOwenControlBatch requires HandlerPO");

    const int ctrlRows = (int)ctrlIdx.size();
    const int ctrlValues = ctrlRows * 16;
    ctrlAvg.assign(ctrlValues, 0.0);
    energyAvg = 0.0;
    if (orientCount <= 0 || ctrlRows <= 0)
        return;
    sobolCount = std::max(sobolCount, orientStart + orientCount);

    Sobol2D sobol(seed);
    std::vector<double> su, sv;
    sobol.generate(sobolCount, su, sv);

    int myStart = m_mpiRank * orientCount / m_mpiSize;
    int myEnd = (m_mpiRank + 1) * orientCount / m_mpiSize;
    int myCount = myEnd - myStart;
    int subChunkMax = std::min(4096, std::max(1, myCount));
    double cosBetaSym = cos(betaSym);
    double weight = sampleWeight;
    m_scattering->PrepareForParallelTrace();

    for (int sc = 0; sc < myCount; sc += subChunkMax)
    {
        int scSize = std::min(subChunkMax, myCount - sc);
        std::vector<PreparedOrientation> chunkPrep(scSize);
        std::vector<double> chunkEnergy(scSize, 0.0);
        ParallelExceptionState parallelError;

#pragma omp parallel
        {
            Particle localParticle = *m_particle;
            Scattering *localScatter =
                m_scattering->CloneFor(&localParticle, &m_incidentLight);
            HandlerPO localHandler(&localParticle, &m_incidentLight,
                                   hp->nTheta, m_scattering->m_wave);
            localHandler.ConfigureForThreadLocalPrepare(*hp, localScatter);
            std::vector<Beam> localBeams;

#pragma omp for schedule(dynamic, 4)
            for (int i = 0; i < scSize; ++i)
            {
                if (parallelError.Failed())
                    continue;
                try
                {
                int idx = orientStart + myStart + sc + i;
                double beta = acos(1.0 - (1.0 - cosBetaSym) * su[idx]);
                double gamma = gammaSym * sv[idx];
                localParticle.Rotate(beta, gamma, 0);
                if (!shadowOff)
                    localScatter->FormShadowBeam(localBeams);
                bool ok = localScatter->ScatterLight(0, 0, localBeams);
                if (ok)
                    localHandler.PrepareBeams(localBeams, weight, chunkPrep[i]);
                else
                    chunkPrep[i].sinZenith = weight;
                chunkEnergy[i] =
                    localScatter->GetIncedentEnergy() * weight;
                localBeams.clear();
                }
                catch (...)
                {
                    parallelError.Capture();
                }
            }

            delete localScatter;
        }
        parallelError.Rethrow();
        for (double e : chunkEnergy)
            energyAvg += e;

        Arr2D ctrlM(nAz + 1, ctrlRows, 4, 4);
        ctrlM.ClearArr();
        Arr2D ctrlMns(nAz + 1, ctrlRows, 4, 4);
        ctrlMns.ClearArr();

        ScatteringRange savedSphere = hp->m_sphere;
        hp->m_sphere = ctrlSphere;
        bool usedGpuCtrl = false;
        if (hp->IsGpuEnabled())
        {
            usedGpuCtrl = true;
            for (int gpuStart = 0; gpuStart < scSize; )
            {
                int gpuBatchSize = hp->SelectGpuOrientationBatchSize(
                    chunkPrep, gpuStart, scSize - gpuStart);
                int gpuEnd = std::min(gpuStart + gpuBatchSize, scSize);
                if (!hp->HandleOrientationsToLocalGpu(
                        chunkPrep, gpuStart, gpuEnd - gpuStart,
                        ctrlM, ctrlMns))
                {
                    usedGpuCtrl = false;
                    ctrlM.ClearArr();
                    ctrlMns.ClearArr();
                    break;
                }
                gpuStart = gpuEnd;
            }
        }

        if (!usedGpuCtrl)
        {
            std::vector<Arr2DC> localJ, localJns;
            if (hp->isCoh)
            {
                Arr2DC tmp(nAz + 1, ctrlRows, 2, 2);
                tmp.ClearArr();
                localJ.push_back(tmp);
                Arr2DC tmpNs(nAz + 1, ctrlRows, 2, 2);
                tmpNs.ClearArr();
                localJns.push_back(tmpNs);
            }

            for (int i = 0; i < scSize; ++i)
            {
                if (!chunkPrep[i].beams.empty())
                    hp->HandleBeamsToLocal(chunkPrep[i], ctrlM, localJ,
                                           hp->isCoh ? &localJns : nullptr);
                if (hp->isCoh && !localJ.empty())
                {
                    double w = chunkPrep[i].sinZenith;
                    HandlerPO::AddToMuellerLocal(localJ, w, ctrlM,
                                                 nAz, ctrlRows - 1);
                    localJ[0].ClearArr();
                    localJns[0].ClearArr();
                }
            }
        }
        hp->m_sphere = savedSphere;

        if (m_mirrorGamma)
            ApplyMirrorGammaMueller(ctrlM, nAz, ctrlRows - 1);

        for (int row = 0; row < ctrlRows; ++row)
        {
            matrix averaged = AzimuthAverageOutputMatrix(
                ctrlM, ctrlSphere, nAz, row);
            double *dst = &ctrlAvg[row * 16];
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    dst[r * 4 + c] += averaged[r][c];
        }
    }

#ifdef USE_MPI
    if (m_mpiSize > 1)
    {
        std::vector<double> send(ctrlAvg.size() + 1, 0.0);
        std::vector<double> recv(ctrlAvg.size() + 1, 0.0);
        for (size_t k = 0; k < ctrlAvg.size(); ++k)
            send[k] = ctrlAvg[k];
        send.back() = energyAvg;
        MPI_Allreduce(send.data(), recv.data(), (int)send.size(),
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (size_t k = 0; k < ctrlAvg.size(); ++k)
            ctrlAvg[k] = recv[k];
        energyAvg = recv.back();
    }
#endif
}

void TracerPOTotal::EvaluateOwenControlSample(
    int nOrient, unsigned int seed, double betaSym, double gammaSym,
    const ScatteringRange &ctrlSphere, const std::vector<int> &ctrlIdx,
    int nAz, std::vector<double> &ctrlAvg, double &energyAvg)
{
    const double weight = nOrient > 0 ? 1.0 / nOrient : 0.0;
    EvaluateOwenControlBatch(0, nOrient, nOrient, seed, weight,
                             betaSym, gammaSym, ctrlSphere, ctrlIdx, nAz,
                             ctrlAvg, energyAvg);
}

int TracerPOTotal::TraceAdaptive(double eps, double betaSym, double gammaSym,
                                  int maxOrientOverride,
                                  bool runFinalDiffraction,
                                  bool strictConvergence)
{
    m_convergenceTarget = eps;
    const int requiredStablePasses = std::max(
        1, m_adaptiveLimits.stablePasses);
    std::cout << "Adaptive mode: target relative change = " << eps << std::endl;
    std::cout << "  Convergence criteria: incoming energy and all Mueller elements"
              << " on control theta rows" << std::endl;
    std::cout << "  beta_sym=" << RadToDeg(betaSym) << " deg, gamma_sym="
              << RadToDeg(gammaSym) << " deg" << std::endl;
    std::vector<std::string> adaptiveLogLines;
    if (m_mpiRank == 0)
    {
        std::ostringstream log;
        log << "===== ADAPTIVE SOBOL =====\n";
        log << "target relative change = " << eps << "\n";
        log << "controls: incoming energy and all Mueller elements on control theta rows\n";
        log << "beta_sym=" << RadToDeg(betaSym) << " deg, gamma_sym="
            << RadToDeg(gammaSym) << " deg\n";
        adaptiveLogLines.push_back(log.str());
    }

    HandlerPO *hp = dynamic_cast<HandlerPO*>(m_handler);
    if (!hp)
    {
        std::cerr << "Error: handler is not HandlerPO" << std::endl;
        TraceFromSobol(1024, betaSym, gammaSym);
        return 1024;
    }

    double cosBetaSym = cos(betaSym);
    int nZen = hp->m_sphere.nZenith;
    int nAz = hp->m_sphere.nAzimuth;
    const int fullControlAz = nAz;
    const int controlAz = AutoFullLowPhiAverageEnabled(hp, eps)
        ? std::max(2, AutoFullFinalAveragePhi(fullControlAz, eps))
        : fullControlAz;
    std::vector<int> ctrlIdx = BuildThetaControlIndices(hp->m_sphere, eps);
    if (ctrlIdx.empty())
        ctrlIdx.push_back(nZen);
    std::vector<int> backConeCtrlRows =
        BuildBackscatterConeControlRows(hp->m_sphere, ctrlIdx);
    std::vector<int> nonBackConeCtrlRows;
    std::vector<char> isBackConeCtrlRow(ctrlIdx.size(), 0);
    if (!backConeCtrlRows.empty())
    {
        for (int row : backConeCtrlRows)
            if (row >= 0 && row < (int)isBackConeCtrlRow.size())
                isBackConeCtrlRow[row] = 1;
    }
    for (int row = 0; row < (int)ctrlIdx.size(); ++row)
        if (row >= (int)isBackConeCtrlRow.size() || !isBackConeCtrlRow[row])
            nonBackConeCtrlRows.push_back(row);
    ScatteringRange fullSphereForControls = hp->m_sphere;
    ScatteringRange ctrlSphere = fullSphereForControls;
    ctrlSphere.nAzimuth = controlAz;
    ctrlSphere.azinuthStep = 2.0 * M_PI / controlAz;
    ctrlSphere.isNonUniform = true;
    ctrlSphere.thetaValues.clear();
    ctrlSphere.thetaValues.reserve(ctrlIdx.size());
    for (int idx : ctrlIdx)
        ctrlSphere.thetaValues.push_back(fullSphereForControls.GetZenith(idx));
    ctrlSphere.nZenith = (int)ctrlSphere.thetaValues.size() - 1;
    ctrlSphere.zenithStart = ctrlSphere.thetaValues.front();
    ctrlSphere.zenithEnd = ctrlSphere.thetaValues.back();
    ctrlSphere.zenithStep = ctrlSphere.nZenith > 0
        ? (ctrlSphere.zenithEnd - ctrlSphere.zenithStart) / ctrlSphere.nZenith
        : 0.0;
    ctrlSphere.ComputeSphereDirections(m_incidentLight);
    nAz = controlAz;
    const int ctrlRows = (int)ctrlIdx.size();
    const int ctrlValues = ctrlRows * 16;
    const bool projectSymmetryInError =
        AutoFullSymmetryErrorProjectionEnabled(betaSym, gammaSym);
    m_lastAdaptiveRowOrient.assign(nZen + 1, 0);
    m_lastAdaptiveRowOrientError.assign(nZen + 1, 0.0);
    m_lastAdaptiveRowOrientElement.assign(nZen + 1, 0);
    m_lastAdaptiveAcceptedOldautoCap = false;
    if (m_mpiRank == 0)
    {
        std::cout << "  Control theta rows (" << ctrlIdx.size() << "):";
        for (int idx : ctrlIdx)
            std::cout << ' ' << std::fixed << std::setprecision(2)
                      << RadToDeg(hp->m_sphere.GetZenith(idx));
        std::cout << " deg" << std::endl;
        if (controlAz != fullControlAz)
            std::cout << "  Control diffraction uses direct N_phi="
                      << controlAz << " matching final low-phi averaging"
                      << std::endl;
        if (projectSymmetryInError)
            std::cout << "  Legacy symmetry-projected error explicitly enabled: ignoring"
                      << " M13/M14/M23/M24/M31/M32/M41/M42"
                      << " (unset MBS_AUTOFULL_SYMMETRY_ERROR_PROJECT for strict checks)"
                      << std::endl;
        if (!backConeCtrlRows.empty())
        {
            std::cout << "  Backscatter convergence uses cone estimator theta>="
                      << AutoFullBackscatterConeMinDeg() << " deg ("
                      << backConeCtrlRows.size()
                      << " control rows; legacy opt-in through"
                      << " MBS_AUTOFULL_BACK_CONE_CONTROL)" << std::endl;
        }
    }

    // =========================================================================
    // Incremental adaptive: each iteration adds NEW orientations only.
    // Sobol property: first 2^m points ⊂ first 2^(m+1) points.
    // So doubling N means adding N new points (indices N..2N-1).
    //
    // Mueller is accumulated with weight=1 per batch, then after each
    // doubling: M_total = (M_old + M_new) / 2 = average of equal-size batches.
    // Since Mueller is additive, we store M_accumulated (sum of all batches)
    // and divide by number of batches when extracting results.
    //
    // Cost: N + N + N + ... = 2N (geometric series) instead of N+2N+4N = 7N.
    // Speedup vs restart: ~3.5× for convergence at the same N.
    // =========================================================================
    // Physics-based orientation estimates use the same ceil/full-grid formula
    // as --oldauto.  This matters for large particles: the oldauto div2 grid is
    // the practical upper reference used by the production column baselines.
    double Dmax = m_particle->MaximalDimention();
    AutoFullOldautoOrientEstimate oldauto16 = AutoFullOldautoEstimate(
        betaSym, gammaSym, m_scattering->m_wave, Dmax, m_ringPoints, 16);
    const int minOrientations = std::max(
        16, m_adaptiveLimits.minOrientations);
    int nStartRaw = oldauto16.total;
    int nOrient = 1;
    while (nOrient * 2 <= nStartRaw) nOrient *= 2;
    if (nOrient < minOrientations) nOrient = minOrientations;
    std::cout << "  Start estimate (oldauto div16): "
              << oldauto16.nBeta << " x " << oldauto16.nGamma
              << " = " << nStartRaw << " -> " << nOrient << " (power of 2)" << std::endl;
    if (m_mpiRank == 0)
    {
        std::ostringstream line;
        line << "Start estimate (oldauto div16): "
             << oldauto16.nBeta << " x " << oldauto16.nGamma
             << " = " << nStartRaw << " -> " << nOrient << " (power of 2)\n";
        adaptiveLogLines.push_back(line.str());
    }

    long long availMB_ad = 2048;
#ifdef __linux__
    {
        std::ifstream meminfo("/proc/meminfo");
        std::string line;
        while (std::getline(meminfo, line)) {
            if (line.find("MemAvailable:") == 0) {
                long long kb = 0;
                sscanf(line.c_str(), "MemAvailable: %lld", &kb);
                if (kb > 0) availMB_ad = kb / 1024;
                break;
            }
        }
    }
#endif
    AutoFullOldautoOrientEstimate oldauto1 = AutoFullOldautoEstimate(
        betaSym, gammaSym, m_scattering->m_wave, Dmax, m_ringPoints, 1);
    int maxFromPhysics = oldauto1.total;
    int maxP2 = 1;
    while (maxP2 <= std::numeric_limits<int>::max() / 2
           && maxP2 * 2 <= maxFromPhysics)
    {
        maxP2 *= 2;
    }

    const int oldautoMaxDiv = strictConvergence
        ? ReadEnvInt("MBS_AUTOFULL_MAX_OLDAUTO_DIV", 2, 1, 64)
        : 0;
    AutoFullOldautoOrientEstimate oldautoCapEstimate =
        oldautoMaxDiv > 0
            ? AutoFullOldautoEstimate(betaSym, gammaSym,
                                      m_scattering->m_wave, Dmax,
                                      m_ringPoints, oldautoMaxDiv)
            : AutoFullOldautoOrientEstimate();
    int oldautoCap = oldautoMaxDiv > 0
        ? std::max(1024, oldautoCapEstimate.total)
        : 0;

    if (maxOrientOverride <= 0 && m_adaptiveLimits.maxOrientations > 0)
        maxOrientOverride = m_adaptiveLimits.maxOrientations;
    int maxOrient;
    bool userMaxOrient = maxOrientOverride > 0;
    bool maxOrientFromOldautoCap = false;
    if (maxOrientOverride > 0) {
        maxOrient = 1;
        while (maxOrient <= std::numeric_limits<int>::max() / 2
               && maxOrient < maxOrientOverride)
        {
            maxOrient *= 2;
        }
        if (maxOrient > maxOrientOverride) maxOrient /= 2;
        if (maxOrient < minOrientations) maxOrient = minOrientations;
        // An explicit user ceiling is authoritative. The oldauto estimate is
        // a default safety cap, not permission to silently lower --maxorient.
    } else {
        if (oldautoCap > 0)
        {
            maxOrient = oldautoCap;
            maxOrientFromOldautoCap = true;
        }
        else
        {
            int maxDiv2 = std::max(1024, maxP2 / 2);
            int p2 = 1;
            while (p2 <= std::numeric_limits<int>::max() / 2
                   && p2 * 2 <= maxDiv2)
            {
                p2 *= 2;
            }
            maxOrient = p2;
        }
    }
    std::cout << "  Oldauto div1 estimate: "
              << oldauto1.nBeta << " x " << oldauto1.nGamma
              << " = " << maxFromPhysics << std::endl;
    if (oldautoCap > 0)
    {
        std::cout << "  Max estimate (oldauto div" << oldautoMaxDiv
                  << "): " << oldautoCapEstimate.nBeta << " x "
                  << oldautoCapEstimate.nGamma << " = "
                  << oldautoCapEstimate.total
                  << (userMaxOrient ? ", with user cap -> " : " -> ")
                  << maxOrient << std::endl;
    }
    else
    {
        std::cout << "  Max estimate (power cap): "
                  << (userMaxOrient ? "user cap -> " : "div2 cap -> ")
                  << maxOrient << std::endl;
    }
    if (nOrient > maxOrient)
    {
        std::cout << "  Start estimate clipped to max orientations: "
                  << nOrient << " -> " << maxOrient << std::endl;
        nOrient = maxOrient;
    }
    std::cerr << "Adaptive: max orientations = " << maxOrient
              << " (" << availMB_ad << " MB available)" << std::endl;
    if (m_mpiRank == 0)
    {
        std::ostringstream line;
        line << "Oldauto div1 estimate: "
             << oldauto1.nBeta << " x " << oldauto1.nGamma
             << " = " << maxFromPhysics << "\n";
        if (oldautoCap > 0)
        {
            line << "Max estimate (oldauto div" << oldautoMaxDiv
                 << "): " << oldautoCapEstimate.nBeta << " x "
                 << oldautoCapEstimate.nGamma << " = "
                 << oldautoCapEstimate.total
                 << (userMaxOrient ? ", with user cap -> " : " -> ")
                 << maxOrient << "\n";
        }
        else
        {
            line << "Max estimate (power cap): "
                 << (userMaxOrient ? "user cap -> " : "div2 cap -> ")
                 << maxOrient << "\n";
        }
        if (nOrient == maxOrient && nStartRaw > maxOrient)
            line << "Start estimate clipped to max orientations\n";
        line << "Adaptive: max orientations = " << maxOrient
             << " (" << availMB_ad << " MB available)\n";
        adaptiveLogLines.push_back(line.str());
    }

    auto t_start = std::chrono::high_resolution_clock::now();

    if (m_owenAverageSeeds.size() >= 2)
    {
        if (m_mpiRank == 0)
        {
            std::ostringstream line;
            line << "Adaptive Owen "
                 << (m_owenAverageSeeds.size() >= 3
                     ? "mean-error"
                     : "inter-seed")
                 << " mode: " << m_owenAverageSeeds.size() << " seeds";
            std::cout << line.str() << std::endl;
            adaptiveLogLines.push_back(line.str() + "\n");
        }

        int totalOrient = 0;
        int convergedCount = 0;
        bool converged = false;
        double lastDMax = std::numeric_limits<double>::infinity();
        double lastDCtrl = std::numeric_limits<double>::infinity();
        int lastCtrlIdx = ctrlIdx.empty() ? 0 : ctrlIdx[0];
        int lastCtrlElement = 0;
        std::vector<int> rowOrientRequirement(nZen + 1, 0);
        std::vector<double> lastRowCtrlErrors(ctrlRows, 0.0);
        std::vector<int> lastRowCtrlElements(ctrlRows, 0);
        std::vector<double> previousSeedMean;
        int previousSeedMeanOrient = 0;
        const bool incrementalOwen =
            EnvFlagEnabled("MBS_AUTOFULL_OWEN_INCREMENTAL", true);
        std::vector<std::vector<double>> seedControlSums;
        std::vector<double> seedEnergySums;
        int accumulatedOrient = 0;
        if (incrementalOwen)
        {
            seedControlSums.assign(m_owenAverageSeeds.size(),
                                   std::vector<double>(ctrlValues, 0.0));
            seedEnergySums.assign(m_owenAverageSeeds.size(), 0.0);
            if (m_mpiRank == 0)
            {
                std::cout << "  Incremental Owen control batches enabled"
                          << " (MBS_AUTOFULL_OWEN_INCREMENTAL=0 disables)"
                          << std::endl;
                adaptiveLogLines.push_back(
                    "Incremental Owen control batches enabled\n");
            }
        }
        for (int iter = 0; iter < 15; ++iter)
        {
            totalOrient = nOrient;
            const int batchStart = incrementalOwen ? accumulatedOrient : 0;
            const int batchCount = incrementalOwen
                ? totalOrient - accumulatedOrient
                : totalOrient;
            std::vector<std::vector<double>> seedCtrl;
            std::vector<double> seedEnergy;
            seedCtrl.reserve(m_owenAverageSeeds.size());
            seedEnergy.reserve(m_owenAverageSeeds.size());

            for (size_t seedIndex = 0; seedIndex < m_owenAverageSeeds.size();
                 ++seedIndex)
            {
                const unsigned int seed = m_owenAverageSeeds[seedIndex];
                if (m_mpiRank == 0)
                {
                    std::cout << "  control seed " << seed
                              << ", N=" << totalOrient;
                    if (incrementalOwen)
                        std::cout << " (new " << std::max(0, batchCount)
                                  << ")";
                    std::cout << std::endl;
                }
                auto seedTimeStart =
                    std::chrono::high_resolution_clock::now();
                std::vector<double> ctrl;
                double energy = 0.0;
                if (incrementalOwen)
                {
                    EvaluateOwenControlBatch(
                        batchStart, batchCount, totalOrient, seed, 1.0,
                        betaSym, gammaSym, ctrlSphere, ctrlIdx, nAz,
                        ctrl, energy);
                    for (int i = 0; i < ctrlValues
                         && i < (int)ctrl.size(); ++i)
                    {
                        seedControlSums[seedIndex][i] += ctrl[i];
                    }
                    seedEnergySums[seedIndex] += energy;
                    ctrl.assign(ctrlValues, 0.0);
                    const double invOrient =
                        totalOrient > 0 ? 1.0 / totalOrient : 0.0;
                    for (int i = 0; i < ctrlValues; ++i)
                        ctrl[i] = seedControlSums[seedIndex][i] * invOrient;
                    energy = seedEnergySums[seedIndex] * invOrient;
                }
                else
                {
                    EvaluateOwenControlSample(
                        totalOrient, seed, betaSym, gammaSym,
                        ctrlSphere, ctrlIdx, nAz, ctrl, energy);
                }
                if (projectSymmetryInError)
                    ProjectMuellerControlValues(ctrl);
                if (m_mpiRank == 0)
                {
                    double seedSec = std::chrono::duration<double>(
                        std::chrono::high_resolution_clock::now()
                        - seedTimeStart).count();
                    std::cout << "    control seed time: "
                              << std::fixed << std::setprecision(2)
                              << seedSec << " s" << std::endl;
                }
                seedCtrl.push_back(ctrl);
                seedEnergy.push_back(energy);
            }
            if (incrementalOwen)
                accumulatedOrient = totalOrient;

            double dEnergyMax = 0.0;
            double dCtrlMax = 0.0;
            int dCtrlMaxIdx = ctrlIdx.empty() ? 0 : ctrlIdx[0];
            int dCtrlMaxElement = 0;
            std::vector<double> rowCtrlErrors(ctrlRows, 0.0);
            std::vector<int> rowCtrlElements(ctrlRows, 0);
            std::vector<double> seedMean(ctrlValues, 0.0);
            if (!seedCtrl.empty())
            {
                for (const std::vector<double> &ctrl : seedCtrl)
                    for (int i = 0; i < ctrlValues && i < (int)ctrl.size(); ++i)
                        seedMean[i] += ctrl[i];
                const double invSeeds = 1.0 / (double)seedCtrl.size();
                for (double &value : seedMean)
                    value *= invSeeds;
            }
            std::vector<double> rowStepErrors(
                ctrlRows, std::numeric_limits<double>::infinity());
            if ((int)previousSeedMean.size() == ctrlValues)
            {
                for (int row = 0; row < ctrlRows; ++row)
                {
                    int element = 0;
                    rowStepErrors[row] = MuellerRowRelativeErrorValues(
                        &seedMean[row * 16],
                        &previousSeedMean[row * 16],
                        &element);
                }
                if (!backConeCtrlRows.empty())
                {
                    std::vector<double> coneMean =
                        AggregateControlRows(seedMean, backConeCtrlRows);
                    std::vector<double> conePrev =
                        AggregateControlRows(previousSeedMean,
                                             backConeCtrlRows);
                    int coneElement = 0;
                    const double coneStep = MuellerRowRelativeErrorValues(
                        coneMean.data(), conePrev.data(), &coneElement);
                    for (int row : backConeCtrlRows)
                    {
                        if (row >= 0 && row < (int)rowStepErrors.size())
                            rowStepErrors[row] = coneStep;
                    }
                }
            }

            const bool meanErrorCriterion = seedCtrl.size() >= 3;
            if (meanErrorCriterion)
            {
                dEnergyMax = OwenMeanScalarRelativeError(seedEnergy);
                int dCtrlMaxRow = 0;
                if (backConeCtrlRows.empty())
                {
                    dCtrlMax = OwenMeanMuellerControlError(
                        seedCtrl, ctrlRows, &rowCtrlErrors, &rowCtrlElements,
                        &dCtrlMaxRow, &dCtrlMaxElement);
                }
                else
                {
                    dCtrlMax = OwenMeanMuellerSelectedRowsError(
                        seedCtrl, ctrlRows, nonBackConeCtrlRows,
                        &rowCtrlErrors, &rowCtrlElements,
                        &dCtrlMaxRow, &dCtrlMaxElement);
                    std::vector<std::vector<double>> coneSeedCtrl =
                        AggregateSeedControlRows(seedCtrl, backConeCtrlRows);
                    int coneElement = 0;
                    const double coneErr = OwenMeanMuellerControlError(
                        coneSeedCtrl, 1, nullptr, nullptr, nullptr,
                        &coneElement);
                    for (int row : backConeCtrlRows)
                    {
                        if (row >= 0 && row < ctrlRows)
                        {
                            rowCtrlErrors[row] = coneErr;
                            rowCtrlElements[row] = coneElement;
                        }
                    }
                    if (std::isfinite(coneErr) && coneErr > dCtrlMax)
                    {
                        dCtrlMax = coneErr;
                        dCtrlMaxRow = backConeCtrlRows.back();
                        dCtrlMaxElement = coneElement;
                    }
                }
                dCtrlMaxIdx = ctrlIdx[std::max(0,
                    std::min(dCtrlMaxRow, (int)ctrlIdx.size() - 1))];
            }
            else
            {
                std::vector<std::vector<double>> coneSeedCtrl;
                if (!backConeCtrlRows.empty())
                    coneSeedCtrl =
                        AggregateSeedControlRows(seedCtrl, backConeCtrlRows);
                for (size_t a = 0; a < seedCtrl.size(); ++a)
                {
                    for (size_t b = a + 1; b < seedCtrl.size(); ++b)
                    {
                        double denomE = std::max(
                            0.5 * (std::fabs(seedEnergy[a])
                                   + std::fabs(seedEnergy[b])),
                            1e-30);
                        dEnergyMax = std::max(
                            dEnergyMax,
                            std::fabs(seedEnergy[a] - seedEnergy[b]) / denomE);

                        for (int row = 0; row < ctrlRows; ++row)
                        {
                            if (row < (int)isBackConeCtrlRow.size()
                                && isBackConeCtrlRow[row])
                                continue;
                            int worstElement = 0;
                            double d = MuellerRowRelativeErrorValues(
                                &seedCtrl[a][row * 16],
                                &seedCtrl[b][row * 16],
                                &worstElement);
                            if (d > dCtrlMax)
                            {
                                dCtrlMax = d;
                                dCtrlMaxIdx = ctrlIdx[row];
                                dCtrlMaxElement = worstElement;
                            }
                            if (d > rowCtrlErrors[row])
                            {
                                rowCtrlErrors[row] = d;
                                rowCtrlElements[row] = worstElement;
                            }
                        }
                        if (!coneSeedCtrl.empty())
                        {
                            int worstElement = 0;
                            double d = MuellerRowRelativeErrorValues(
                                coneSeedCtrl[a].data(),
                                coneSeedCtrl[b].data(),
                                &worstElement);
                            for (int row : backConeCtrlRows)
                            {
                                if (row >= 0 && row < ctrlRows
                                    && d > rowCtrlErrors[row])
                                {
                                    rowCtrlErrors[row] = d;
                                    rowCtrlElements[row] = worstElement;
                                }
                            }
                            if (d > dCtrlMax)
                            {
                                dCtrlMax = d;
                                dCtrlMaxIdx = ctrlIdx[backConeCtrlRows.back()];
                                dCtrlMaxElement = worstElement;
                            }
                        }
                    }
                }
            }

            for (int row = 0; row < ctrlRows; ++row)
            {
                const int globalRow = ctrlIdx[row];
                if (globalRow < 0 || globalRow >= (int)rowOrientRequirement.size())
                    continue;
                const double thetaDeg =
                    RadToDeg(hp->m_sphere.GetZenith(globalRow));
                const int thetaZone = AutoFullThetaZone(thetaDeg);
                const double rowOrientEps =
                    eps * AutoFullRowOrientSafetyFactor(thetaZone);
                if (rowOrientRequirement[globalRow] == 0
                    && rowOrientEps > 0.0
                    && previousSeedMeanOrient > 0
                    && rowCtrlErrors[row] <= rowOrientEps
                    && rowStepErrors[row] <= rowOrientEps)
                {
                    rowOrientRequirement[globalRow] = previousSeedMeanOrient;
                }
            }
            previousSeedMean.swap(seedMean);
            previousSeedMeanOrient = totalOrient;
            lastRowCtrlErrors = rowCtrlErrors;
            lastRowCtrlElements = rowCtrlElements;

            double dMax = std::max(dEnergyMax, dCtrlMax);
            lastDMax = dMax;
            lastDCtrl = dCtrlMax;
            lastCtrlIdx = dCtrlMaxIdx;
            lastCtrlElement = dCtrlMaxElement;
            if (m_mpiRank == 0)
            {
                std::ostringstream line;
                line << std::fixed << std::setprecision(2);
                line << "  N=" << totalOrient
                     << (meanErrorCriterion
                         ? "  OwenMeanError dE="
                         : "  OwenInterseed dE=")
                     << dEnergyMax * 100.0 << "%"
                     << (meanErrorCriterion
                         ? "  dMuellerMean="
                         : "  dMuellerCtrl=")
                     << dCtrlMax * 100.0 << "%"
                     << "@" << RadToDeg(hp->m_sphere.GetZenith(dCtrlMaxIdx))
                     << "deg/" << MuellerElementName(dCtrlMaxElement)
                     << "  max=" << dMax * 100.0 << "%";
                std::cout << line.str() << std::endl;
                adaptiveLogLines.push_back(line.str() + "\n");
            }

            if (dMax <= eps)
                ++convergedCount;
            else
                convergedCount = 0;
            RecordConvergenceStep(
                "orientations", iter, totalOrient,
                (int)m_owenAverageSeeds.size(), dMax,
                RadToDeg(hp->m_sphere.GetZenith(dCtrlMaxIdx)),
                dCtrlMaxElement,
                dMax <= eps ? "pass" : "refine");

            if (convergedCount >= requiredStablePasses)
            {
                converged = true;
                if (m_mpiRank == 0)
                {
                    std::ostringstream line;
                    line << "Converged at N=" << totalOrient
                         << " by Owen "
                         << (meanErrorCriterion ? "mean-error"
                                                : "inter-seed")
                         << " controls within "
                         << eps * 100.0 << "%";
                    std::cout << line.str() << std::endl;
                    adaptiveLogLines.push_back(line.str() + "\n");
                }
                break;
            }

            long long nextOrient = (long long)nOrient * 2;
            if (nextOrient > maxOrient && totalOrient < maxOrient)
            {
                nOrient = maxOrient;
                continue;
            }
            nOrient = ClampOrientCount(nextOrient);
            if (nOrient > maxOrient)
            {
                if (m_mpiRank == 0)
                {
                    std::ostringstream line;
                    line << "WARNING: Max orientations reached (N="
                         << totalOrient << ", limit=" << maxOrient
                         << "). Owen "
                         << (m_owenAverageSeeds.size() >= 3
                             ? "mean-error"
                             : "inter-seed")
                         << " target accuracy "
                         << eps * 100.0 << "% may not be achieved.";
                    std::cout << line.str() << std::endl;
                    adaptiveLogLines.push_back(line.str() + "\n");
                    if (maxOrientFromOldautoCap)
                    {
                        std::cout << "  Reached oldauto-based maximum; not"
                                  << " extrapolating beyond that cap."
                                  << std::endl;
                        adaptiveLogLines.push_back(
                            "Reached oldauto-based maximum; not extrapolating beyond that cap.\n");
                    }
                    else
                    {
                        std::cout << "  To improve: use --maxorient "
                                  << maxOrient * 2 << std::endl;
                        adaptiveLogLines.push_back(
                            "  To improve: use --maxorient "
                            + std::to_string(maxOrient * 2) + "\n");
                    }
                }
                break;
            }
        }

        if (!converged && strictConvergence)
        {
            const bool acceptOldautoCap =
                maxOrientFromOldautoCap
                && totalOrient >= maxOrient
                && EnvFlagEnabled("MBS_AUTOFULL_ACCEPT_OLDAUTO_CAP", false);
            if (acceptOldautoCap)
            {
                m_lastAdaptiveAcceptedOldautoCap = true;
                std::ostringstream warn;
                warn << "WARNING: AUTOFULL reached oldauto-based orientation"
                     << " cap N=" << totalOrient << " before Owen target "
                     << eps * 100.0 << "% (last max="
                     << lastDMax * 100.0 << "%, control="
                     << lastDCtrl * 100.0 << "% @ "
                     << RadToDeg(hp->m_sphere.GetZenith(lastCtrlIdx))
                     << "deg/" << MuellerElementName(lastCtrlElement)
                     << "). Accepting oldauto cap as the maximum estimate.";
                std::cout << warn.str() << std::endl;
                adaptiveLogLines.push_back(warn.str() + "\n");
            }
            else
            {
                std::ostringstream err;
                err << "ERROR: AUTOFULL did not converge to target "
                    << eps * 100.0 << "% by N=" << totalOrient
                    << " (last max=" << lastDMax * 100.0
                    << "%, control=" << lastDCtrl * 100.0 << "% @ "
                    << RadToDeg(hp->m_sphere.GetZenith(lastCtrlIdx)) << "deg/"
                    << MuellerElementName(lastCtrlElement) << "). "
                    << "Increase --maxorient to at least " << maxOrient * 2
                    << " or loosen --autofull eps.";
                std::cerr << err.str() << std::endl;
                std::exit(2);
            }
        }

        for (int row = 0; row < ctrlRows; ++row)
        {
            const int globalRow = ctrlIdx[row];
            if (globalRow < 0 || globalRow >= (int)rowOrientRequirement.size())
                continue;
            if (rowOrientRequirement[globalRow] == 0)
                rowOrientRequirement[globalRow] = totalOrient;
            if (row < (int)lastRowCtrlErrors.size())
                m_lastAdaptiveRowOrientError[globalRow] = lastRowCtrlErrors[row];
            if (row < (int)lastRowCtrlElements.size())
                m_lastAdaptiveRowOrientElement[globalRow] =
                    lastRowCtrlElements[row];
        }
        m_lastAdaptiveRowOrient.swap(rowOrientRequirement);
        if (m_mpiRank == 0)
        {
            const std::string orientMapName =
                m_resultDirName + "_autofull_orient_by_theta.dat";
            std::ofstream orientMap(orientMapName.c_str(), std::ios::out);
            orientMap << std::setprecision(10);
            orientMap << "# theta_deg required_Norient last_orient_error_pct"
                      << " worst_element\n";
            for (int t = 0; t <= nZen; ++t)
            {
                orientMap << RadToDeg(hp->m_sphere.GetZenith(t)) << ' '
                          << m_lastAdaptiveRowOrient[t] << ' '
                          << m_lastAdaptiveRowOrientError[t] * 100.0 << ' '
                          << MuellerElementName(
                                 m_lastAdaptiveRowOrientElement[t])
                          << '\n';
            }
            orientMap.close();
            std::cout << "Adaptive orientation map written: "
                      << orientMapName << std::endl;
        }
        if (m_mpiRank == 0)
        {
            int knownRows = 0;
            int minRowOrient = totalOrient;
            long long sumRowOrient = 0;
            for (int value : m_lastAdaptiveRowOrient)
            {
                if (value <= 0)
                    continue;
                ++knownRows;
                minRowOrient = std::min(minRowOrient, value);
                sumRowOrient += value;
            }
            if (knownRows > 0)
            {
                std::ostringstream line;
                line << "Adaptive per-theta N_orient map: "
                     << knownRows << " rows, min=" << minRowOrient
                     << ", avg=" << std::fixed << std::setprecision(1)
                     << ((double)sumRowOrient / knownRows)
                     << ", max=" << totalOrient;
                std::cout << line.str() << std::endl;
                adaptiveLogLines.push_back(line.str() + "\n");
            }
        }

        if (!runFinalDiffraction)
        {
            if (m_mpiRank == 0)
                std::cout << "\nFinal diffraction deferred; selected N="
                          << totalOrient << std::endl;
            return totalOrient;
        }

        if (m_mpiRank == 0)
            std::cout << "\nFinal: full diffraction with N=" << totalOrient
                      << " (first Owen seed; averaging is autofull-only)"
                      << std::endl;
        hp->M.ClearArr();
        hp->M_noshadow.ClearArr();
        hp->CleanJ();
        m_incomingEnergy = 0;
        hp->m_outputEnergy = 0;
        hp->m_extinctionCrossSectionOt = 0;
        hp->m_hasExtinctionOt = false;
        TraceFromSobolSeed(totalOrient, m_owenAverageSeeds.front(),
                           betaSym, gammaSym);
        return totalOrient;
    }

    // Generate ALL Sobol points up to maxOrient at once (deterministic)
    Sobol2D sobol_gen(42);
    std::vector<double> su_all, sv_all;
    sobol_gen.generate(maxOrient, su_all, sv_all);

    // Accumulator: each batch is averaged internally and then combined by
    // its orientation count. This keeps incremental Sobol identical to a
    // single run over the same prefix of Sobol points.
    hp->M.ClearArr();
    hp->M_noshadow.ClearArr();
    hp->CleanJ();
    m_incomingEnergy = 0;
    hp->m_outputEnergy = 0;
    hp->m_extinctionCrossSectionOt = 0;
    hp->m_hasExtinctionOt = false;

    double prevEnergy = 0;
    std::vector<double> prevCtrl(ctrlValues, 0.0);
    int totalOrient = 0;
    int convergedCount = 0;
    bool converged = false;
    double lastDMax = std::numeric_limits<double>::infinity();
    double lastDCtrl = std::numeric_limits<double>::infinity();
    int lastCtrlIdx = ctrlIdx.empty() ? 0 : ctrlIdx[0];
    int lastCtrlElement = 0;
    std::vector<double> ctrlWeightedSum(ctrlValues, 0.0); // sum(batch_average * batch_size)
    double energyWeightedSum = 0.0;                           // sum(batch_average * batch_size)

    for (int iter = 0; iter < 15; ++iter)
    {
        int batchStart = totalOrient;
        int batchEnd = nOrient;  // nOrient = target total after this iteration
        int batchSize = batchEnd - batchStart;

        if (batchSize <= 0) { nOrient *= 2; continue; }

        // Trace and diffract only the NEW orientations [batchStart..batchEnd).
        // Each batch is averaged with weight=1/batchSize, then accumulated
        // with weight=batchSize/totalOrient. This is required because after
        // the first doubling the added Sobol batches have different sizes.

        // MPI: each rank processes a subset of this batch
        int myBatchStart = m_mpiRank * batchSize / m_mpiSize;
        int myBatchEnd = (m_mpiRank + 1) * batchSize / m_mpiSize;
        int myBatchSize = myBatchEnd - myBatchStart;

        double batchEnergy = 0;
        double weight = 1.0 / batchSize;  // global weight (not per-rank)
        std::vector<Beam> outBeams;

        std::vector<double> batchCtrl(ctrlValues, 0.0);

        // Process in sub-chunks to limit memory (max 4096 per chunk)
        int subChunkMax = std::min(4096, myBatchSize);
        for (int sc = 0; sc < myBatchSize; sc += subChunkMax)
        {
            int scSize = std::min(subChunkMax, myBatchSize - sc);
            std::vector<PreparedOrientation> chunkPrep(scSize);

            for (int i = 0; i < scSize; ++i)
            {
                int idx = batchStart + myBatchStart + sc + i;
                double beta  = acos(1.0 - (1.0 - cosBetaSym) * su_all[idx]);
                double gamma = gammaSym * sv_all[idx];
                m_particle->Rotate(beta, gamma, 0);
                if (!shadowOff) m_scattering->FormShadowBeam(outBeams);
                bool ok = m_scattering->ScatterLight(0, 0, outBeams);
                if (ok)
                    hp->PrepareBeams(outBeams, weight, chunkPrep[i]);
                else
                    chunkPrep[i].sinZenith = weight;
                batchEnergy += m_scattering->GetIncedentEnergy() * weight;
                outBeams.clear();
            }

            // Diffract selected control theta rows only, but keep all 16
            // Mueller elements. This uses the same local Mueller path as the
            // final calculation on a tiny non-uniform theta grid.
            Arr2D ctrlM(nAz + 1, ctrlRows, 4, 4);
            ctrlM.ClearArr();
            Arr2D ctrlMns(nAz + 1, ctrlRows, 4, 4);
            ctrlMns.ClearArr();

            ScatteringRange savedSphere = hp->m_sphere;
            hp->m_sphere = ctrlSphere;
            bool usedGpuCtrl = false;
            if (hp->IsGpuEnabled())
            {
                usedGpuCtrl = true;
                for (int gpuStart = 0; gpuStart < scSize; )
                {
                    int gpuBatchSize = hp->SelectGpuOrientationBatchSize(
                        chunkPrep, gpuStart, scSize - gpuStart);
                    int gpuEnd = std::min(gpuStart + gpuBatchSize, scSize);
                    if (!hp->HandleOrientationsToLocalGpu(
                            chunkPrep, gpuStart, gpuEnd - gpuStart,
                            ctrlM, ctrlMns))
                    {
                        usedGpuCtrl = false;
                        ctrlM.ClearArr();
                        ctrlMns.ClearArr();
                        break;
                    }
                    gpuStart = gpuEnd;
                }
            }

            if (!usedGpuCtrl)
            {
                std::vector<Arr2DC> localJ, localJns;
                if (hp->isCoh)
                {
                    Arr2DC tmp(nAz + 1, ctrlRows, 2, 2);
                    tmp.ClearArr();
                    localJ.push_back(tmp);
                    Arr2DC tmpNs(nAz + 1, ctrlRows, 2, 2);
                    tmpNs.ClearArr();
                    localJns.push_back(tmpNs);
                }

                for (int i = 0; i < scSize; ++i)
                {
                    if (!chunkPrep[i].beams.empty())
                        hp->HandleBeamsToLocal(chunkPrep[i], ctrlM, localJ,
                                               hp->isCoh ? &localJns : nullptr);
                    if (hp->isCoh && !localJ.empty())
                    {
                        double w = chunkPrep[i].sinZenith;
                        HandlerPO::AddToMuellerLocal(localJ, w, ctrlM,
                                                     nAz, ctrlRows - 1);
                        localJ[0].ClearArr();
                        localJns[0].ClearArr();
                    }
                }
            }
            hp->m_sphere = savedSphere;

            if (m_mirrorGamma)
                ApplyMirrorGammaMueller(ctrlM, nAz, ctrlRows - 1);

            for (int row = 0; row < ctrlRows; ++row)
            {
                matrix averaged = AzimuthAverageOutputMatrix(
                    ctrlM, ctrlSphere, nAz, row);
                double *dst = &batchCtrl[row * 16];
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c)
                        dst[r * 4 + c] += averaged[r][c];
            }
        } // sub-chunks

        // MPI: reduce only THIS BATCH's control points + energy.
#ifdef USE_MPI
        if (m_mpiSize > 1)
        {
            std::vector<double> send(batchCtrl.size() + 1, 0.0);
            std::vector<double> recv(batchCtrl.size() + 1, 0.0);
            for (size_t k = 0; k < batchCtrl.size(); ++k)
                send[k] = batchCtrl[k];
            send.back() = batchEnergy;
            MPI_Allreduce(send.data(), recv.data(), (int)send.size(),
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            for (size_t k = 0; k < batchCtrl.size(); ++k)
                batchCtrl[k] = recv[k];
            batchEnergy = recv.back();
        }
#endif

        for (size_t k = 0; k < ctrlWeightedSum.size(); ++k)
            ctrlWeightedSum[k] += batchCtrl[k] * batchSize;
        energyWeightedSum += batchEnergy * batchSize;

        totalOrient = batchEnd;

        // Control points from a tiny full-Mueller grid.
        double energyAvg = (totalOrient > 0) ? energyWeightedSum / totalOrient : 0;
        std::vector<double> ctrlAvg(ctrlValues, 0.0);
        for (size_t k = 0; k < ctrlAvg.size(); ++k)
            ctrlAvg[k] = (totalOrient > 0) ? ctrlWeightedSum[k] / totalOrient : 0.0;

        auto relChange = [](double current, double previous) {
            double denom = std::max(std::fabs(previous), 1e-30);
            return std::fabs(current - previous) / denom;
        };

        double dEnergy = (iter > 0) ? relChange(energyAvg, prevEnergy) : 1.0;
        double dCtrlMax = iter > 0 ? 0.0 : 1.0;
        int dCtrlMaxIdx = ctrlIdx.empty() ? 0 : ctrlIdx[0];
        int dCtrlMaxElement = 0;
        for (int row = 0; row < ctrlRows; ++row)
        {
            if (row < (int)isBackConeCtrlRow.size()
                && isBackConeCtrlRow[row])
                continue;
            int worstElement = 0;
            double d = (iter > 0)
                ? MuellerRowRelativeErrorValues(&ctrlAvg[row * 16],
                                                &prevCtrl[row * 16],
                                                &worstElement)
                : 1.0;
            if (d > dCtrlMax)
            {
                dCtrlMax = d;
                dCtrlMaxIdx = ctrlIdx[row];
                dCtrlMaxElement = worstElement;
            }
        }
        if (iter > 0 && !backConeCtrlRows.empty())
        {
            std::vector<double> coneCurrent =
                AggregateControlRows(ctrlAvg, backConeCtrlRows);
            std::vector<double> conePrevious =
                AggregateControlRows(prevCtrl, backConeCtrlRows);
            int worstElement = 0;
            double d = MuellerRowRelativeErrorValues(
                coneCurrent.data(), conePrevious.data(), &worstElement);
            if (d > dCtrlMax)
            {
                dCtrlMax = d;
                dCtrlMaxIdx = ctrlIdx[backConeCtrlRows.back()];
                dCtrlMaxElement = worstElement;
            }
        }
        double dMax = std::max(dEnergy, dCtrlMax);
        lastDMax = dMax;
        lastDCtrl = dCtrlMax;
        lastCtrlIdx = dCtrlMaxIdx;
        lastCtrlElement = dCtrlMaxElement;

        std::cout << std::fixed << std::setprecision(2);
        if (m_mpiRank == 0)
        {
            std::ostringstream line;
            line << std::fixed << std::setprecision(2);
            line << "  N=" << totalOrient << " (+" << batchSize << ")"
                 << "  dE=" << dEnergy*100 << "%"
                 << "  dMuellerCtrl=" << dCtrlMax*100 << "%"
                 << "@" << RadToDeg(hp->m_sphere.GetZenith(dCtrlMaxIdx))
                 << "deg/" << MuellerElementName(dCtrlMaxElement)
                 << "  max=" << dMax*100 << "%";
            std::cout << line.str() << std::endl;
            adaptiveLogLines.push_back(line.str() + "\n");
        }

        bool all_ok = (dMax <= eps && iter > 0);

        if (all_ok)
            convergedCount++;
        else
            convergedCount = 0;
        RecordConvergenceStep(
            "orientations", iter, totalOrient, batchSize, dMax,
            RadToDeg(hp->m_sphere.GetZenith(dCtrlMaxIdx)),
            dCtrlMaxElement, all_ok ? "pass" : "refine");

        if (convergedCount >= requiredStablePasses)
        {
            converged = true;
            if (m_mpiRank == 0)
            {
                std::ostringstream line;
                line << "Converged at N=" << totalOrient
                     << " (energy + all Mueller elements on " << ctrlIdx.size()
                     << " theta controls within "
                     << eps*100
                     << "% for " << requiredStablePasses
                     << " consecutive step(s))";
                std::cout << line.str() << std::endl;
                adaptiveLogLines.push_back(line.str() + "\n");
            }
            break;
        }

        prevEnergy = energyAvg;
        prevCtrl = ctrlAvg;
        nOrient *= 2;
        if (nOrient > maxOrient)
        {
            if (m_mpiRank == 0)
            {
                std::ostringstream line;
                if (all_ok)
                {
                    line << "WARNING: Max orientations reached (N=" << totalOrient
                         << ", limit=" << maxOrient << ") after one passing step. "
                         << "Two-step confirmation of target accuracy "
                         << eps*100 << "% was not completed.";
                }
                else
                {
                    line << "WARNING: Max orientations reached (N=" << totalOrient
                         << ", limit=" << maxOrient << "). Target accuracy "
                         << eps*100 << "% may not be achieved.";
                }
                std::cout << line.str() << std::endl;
                std::cout << "  To improve: use --maxorient " << maxOrient*2 << std::endl;
                adaptiveLogLines.push_back(line.str() + "\n");
                adaptiveLogLines.push_back("  To improve: use --maxorient " + std::to_string(maxOrient*2) + "\n");
            }
            break;
        }
    }

    if (!converged && strictConvergence)
    {
        std::ostringstream err;
        err << "ERROR: AUTOFULL did not converge to target "
            << eps * 100.0 << "% by N=" << totalOrient
            << " (last max=" << lastDMax * 100.0
            << "%, control=" << lastDCtrl * 100.0 << "% @ "
            << RadToDeg(hp->m_sphere.GetZenith(lastCtrlIdx)) << "deg/"
            << MuellerElementName(lastCtrlElement) << "). "
            << "Increase --maxorient to at least " << maxOrient * 2
            << " or loosen --autofull eps.";
        std::cerr << err.str() << std::endl;
        std::exit(2);
    }

    // =================================================================
    // FINAL PHASE: full diffraction via TraceFromSobol (re-traces).
    // Sobol is deterministic → same orientations. Phase 1 = ~5% overhead.
    // Memory = O(chunkSize), not O(totalOrient). Chunked + OpenMP + MPI.
    // =================================================================
    if (!runFinalDiffraction)
    {
        if (m_mpiRank == 0)
            std::cout << "\nFinal diffraction deferred; selected N="
                      << totalOrient << std::endl;
        return totalOrient;
    }

    if (m_mpiRank == 0)
        std::cout << "\nFinal: full diffraction with N=" << totalOrient
                  << " (re-tracing, chunked)" << std::endl;

    hp->M.ClearArr();
    hp->M_noshadow.ClearArr();
    hp->CleanJ();
    m_incomingEnergy = 0;
    hp->m_outputEnergy = 0;
    hp->m_extinctionCrossSectionOt = 0;
    hp->m_hasExtinctionOt = false;

    // TraceFromSobol handles file output, MPI reduce, and cleanup.
    // It re-traces all orientations (Sobol deterministic → same results)
    // using chunked streaming (O(chunkSize) memory).
    TraceFromSobol(totalOrient, betaSym, gammaSym);

    if (m_mpiRank == 0) {
        auto t_end = std::chrono::high_resolution_clock::now();
        double total_sec = std::chrono::duration<double>(t_end - t_start).count();
        std::cout << "Adaptive total time: " << std::fixed
                  << std::setprecision(1) << total_sec << " s" << std::endl;
        std::ostringstream log;
        log << "\n";
        for (const std::string &line : adaptiveLogLines)
            log << line;
        log << "Adaptive total time: " << std::fixed << std::setprecision(1)
            << total_sec << " s\n";
        AppendTextLog(log.str());
    }
    return totalOrient;
}

// =============================================================================
// TraceAutoFull: 3D sequential optimization n → N_phi → N_orient
// Step 1: Increase n until Q_sca converges (cheap: tracing only)
// Step 2: Increase N_phi until Q_sca converges (medium: new grid)
// Step 3: Increase N_orient until all control points converge (expensive)
// =============================================================================
static bool SameAdaptiveGrid(const ScatteringRange &left,
                             const ScatteringRange &right)
{
    if (left.nAzimuth != right.nAzimuth
        || left.nZenith != right.nZenith
        || left.isNonUniform != right.isNonUniform)
        return false;
    for (int row = 0; row <= left.nZenith; ++row)
        if (std::fabs(left.GetZenith(row) - right.GetZenith(row)) > 1e-12)
            return false;
    return true;
}

void TracerPOTotal::TraceAutoConverged(double eps, double betaSym,
                                       double gammaSym,
                                       int maxOrientOverride,
                                       bool tuneReflections,
                                       bool regularEulerFinal,
                                       ScatteringRange &conus,
                                       bool tunePhi,
                                       bool tuneTheta)
{
    HandlerPO *handler = dynamic_cast<HandlerPO*>(m_handler);
    if (!handler)
        throw std::runtime_error("automatic convergence requires HandlerPO");
    if (!(eps > 0.0 && eps < 1.0))
        throw std::invalid_argument(
            "automatic convergence tolerance must be in (0, 1)");

    const std::string mode = regularEulerFinal
        ? "diffraction-autofull"
        : (tuneReflections ? "autofull" : "auto");
    ResetConvergenceReport(mode, eps);

    const double reflectionEps = ResolveAdaptiveTolerance(
        m_adaptiveLimits.reflectionTolerance, eps);
    const double phiEps = ResolveAdaptiveTolerance(
        m_adaptiveLimits.phiTolerance, eps);
    const double thetaEps = ResolveAdaptiveTolerance(
        m_adaptiveLimits.thetaTolerance, eps);
    const double orientationEps = ResolveAdaptiveTolerance(
        m_adaptiveLimits.orientationTolerance, eps);

    if (tuneTheta
        && handler->m_sphere.nZenith + 1 < m_adaptiveLimits.minThetaPoints)
    {
        ScatteringRange startup(
            0.0, M_PI,
            tunePhi ? m_adaptiveLimits.minPhiPoints
                    : handler->m_sphere.nAzimuth,
            m_adaptiveLimits.minThetaPoints - 1);
        handler->SetScatteringSphere(startup);
    }
    else if (tunePhi
             && handler->m_sphere.nAzimuth < m_adaptiveLimits.minPhiPoints)
    {
        ScatteringRange startup = handler->m_sphere;
        startup.nAzimuth = m_adaptiveLimits.minPhiPoints;
        startup.azinuthStep = M_2PI / startup.nAzimuth;
        startup.ComputeSphereDirections(m_incidentLight);
        handler->SetScatteringSphere(startup);
    }

    int selectedReflection = m_scattering->GetMaxReflections();
    int selectedPhi = handler->m_sphere.nAzimuth;
    int selectedOrient = 0;
    int selectedBeta = 0;
    int selectedGamma = 0;
    int basePilot = std::max(
        m_adaptiveLimits.minPilotOrientations,
        std::min(256, AutoFullPilotOrientations(orientationEps)));
    ScatteringRange previousGrid = handler->m_sphere;
    int previousReflection = -1;
    int previousPhi = -1;
    int previousOrient = -1;
    int previousBeta = -1;
    int previousGamma = -1;
    bool jointStable = false;
    int pilot = basePilot;
    int jointPilotCap = m_adaptiveLimits.maxPilotOrientations;
    if (maxOrientOverride > 0)
        jointPilotCap = std::min(jointPilotCap,
                                 std::max(16, maxOrientOverride));
    basePilot = std::min(basePilot, jointPilotCap);
    pilot = basePilot;

    std::cout << "===== " << mode << ": unified convergence ====="
              << std::endl;
    std::cout << "Order: " << (tuneReflections ? "n -> " : "")
              << (tunePhi ? "N_alpha(N_phi) -> " : "N_alpha fixed -> ")
              << (tuneTheta ? "theta -> " : "theta fixed -> ")
              << "orientations; repeat until all settings stabilize" << std::endl;
    std::cout << "Targets: n=" << reflectionEps * 100.0
              << "%, alpha=" << phiEps * 100.0
              << "%, theta=" << thetaEps * 100.0
              << "%, orientations=" << orientationEps * 100.0
              << "%" << std::endl;

    const int maxJointSweeps = std::max(
        1, m_adaptiveLimits.maxJointSweeps);
    for (int sweep = 1; sweep <= maxJointSweeps; ++sweep)
    {
        std::cout << "--- joint sweep " << sweep
                  << ", pilot N=" << pilot << " ---" << std::endl;

        if (tuneReflections)
        {
            selectedReflection = TraceAdaptiveReflections(
                reflectionEps, pilot, betaSym, gammaSym,
                m_adaptiveLimits.minReflections);
        }
        if (tunePhi)
            selectedPhi = TraceAdaptivePhi(
                phiEps, pilot, betaSym, gammaSym);
        if (tuneTheta)
            TraceAdaptiveTheta(pilot, betaSym, gammaSym,
                               thetaEps, 12, true);
        conus = handler->m_sphere;

        if (regularEulerFinal)
        {
            TraceAdaptiveEuler(orientationEps, betaSym, gammaSym,
                               maxOrientOverride, selectedBeta,
                               selectedGamma, false);
            const long long count = (long long)selectedBeta * selectedGamma;
            selectedOrient = count > std::numeric_limits<int>::max()
                ? std::numeric_limits<int>::max() : (int)count;
        }
        else
        {
            selectedOrient = TraceAdaptive(
                orientationEps, betaSym, gammaSym,
                maxOrientOverride, false, true);
            selectedBeta = 0;
            selectedGamma = 0;
        }

        const int nextPilot = std::max(
            16, std::min(selectedOrient, jointPilotCap));
        const bool sameSettings = previousReflection == selectedReflection
            && previousPhi == selectedPhi
            && previousOrient == selectedOrient
            && previousBeta == selectedBeta
            && previousGamma == selectedGamma
            && SameAdaptiveGrid(previousGrid, handler->m_sphere);
        const bool pilotFixedPoint = nextPilot == pilot;
        const bool same = sameSettings || pilotFixedPoint;
        RecordConvergenceStep("joint_sweep", sweep, selectedReflection,
                              selectedPhi, same ? 0.0 : 1.0,
                              0.0, -1, same ? "pass" : "refine");
        RecordConvergenceStep("joint_orientations", sweep,
                              regularEulerFinal ? selectedBeta : selectedOrient,
                              regularEulerFinal ? selectedGamma : 0,
                              same ? 0.0 : 1.0, 0.0, -1,
                              same ? "pass" : "refine");
        if (same)
        {
            std::cout << "Joint settings accepted: "
                      << (pilotFixedPoint
                          ? "orientation pilot reached its coupled fixed point"
                          : "all selected values repeated on the refined pilot")
                      << std::endl;
            jointStable = true;
            break;
        }
        previousReflection = selectedReflection;
        previousPhi = selectedPhi;
        previousOrient = selectedOrient;
        previousBeta = selectedBeta;
        previousGamma = selectedGamma;
        previousGrid = handler->m_sphere;
        pilot = nextPilot;
    }

    if (!jointStable)
    {
        std::ostringstream message;
        message << "automatic n/phi/theta/orientation settings did not stabilize after "
                << maxJointSweeps
                << " joint sweeps; increase controller.max_joint_sweeps, adaptive limits, or tolerances";
        throw std::runtime_error(message.str());
    }

    if (regularEulerFinal)
    {
        TraceFromEulerQuadrature(selectedBeta, selectedGamma,
                                 betaSym, gammaSym);
    }
    else
    {
        RecordConvergenceStep("orientations", 0, selectedOrient, 0,
                              0.0, 0.0, -1, "selected");
        TraceFromSobol(selectedOrient, betaSym, gammaSym);
    }

    conus = handler->m_sphere;
    std::cout << "Unified convergence selected: n="
              << m_scattering->GetMaxReflections()
              << ", N_phi=" << handler->m_sphere.nAzimuth
              << ", N_theta=" << handler->m_sphere.nZenith + 1
              << std::endl;
    std::cout << "Convergence report: " << m_resultDirName
              << "_convergence.tsv" << std::endl;
}
