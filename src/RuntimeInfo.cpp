#include "RuntimeInfo.h"

#include "GpuSupport.h"

#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

extern char **environ;

namespace
{

std::string Trim(const std::string &value)
{
    const size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return std::string();
    const size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

long long ReadNamedKb(const char *path, const std::string &key)
{
    std::ifstream input(path);
    std::string line;
    while (std::getline(input, line))
    {
        std::istringstream fields(line);
        std::string name;
        long long value = 0;
        if (fields >> name >> value && name == key)
            return value;
    }
    return 0;
}

std::string ReadCpuModel()
{
    std::ifstream input("/proc/cpuinfo");
    std::string line;
    while (std::getline(input, line))
    {
        if (line.compare(0, 10, "model name") == 0
            || line.compare(0, 8, "Hardware") == 0)
        {
            const size_t colon = line.find(':');
            if (colon != std::string::npos)
                return Trim(line.substr(colon + 1));
        }
    }
    return "unknown";
}

int ReadPhysicalCoreCount()
{
    std::ifstream input("/proc/cpuinfo");
    std::set<std::string> cores;
    std::string line;
    std::string physicalId = "0";
    while (std::getline(input, line))
    {
        const size_t colon = line.find(':');
        if (colon == std::string::npos)
            continue;
        const std::string key = Trim(line.substr(0, colon));
        const std::string value = Trim(line.substr(colon + 1));
        if (key == "physical id")
            physicalId = value;
        else if (key == "core id")
            cores.insert(physicalId + ":" + value);
    }
    return (int)cores.size();
}

double ToMib(long long kib)
{
    return (double)kib / 1024.0;
}

bool StartsWith(const std::string &value, const char *prefix)
{
    const std::string p(prefix);
    return value.compare(0, p.size(), p) == 0;
}

bool IsExperimentalEnvironmentName(const std::string &name)
{
    static const char *const prefixes[] = {
        "MBS_AUTOFULL_", "MBS_OLDAUTO_", "MBS_OLDAUTOFULL_",
        "MBS_FFT_", "MBS_DEBUG_ORIENT_", "MBS_AGGREGATE_PART_SHADOW",
        "MBS_GPU_TRACE_"
    };
    for (const char *prefix : prefixes)
        if (StartsWith(name, prefix))
            return true;

    static const char *const exact[] = {
        "MBS_DISABLE_TRACK_IDS", "MBS_FORCE_TRACK_IDS", "MBS_GPU_EPS1",
        "MBS_GPU_MULTI_K_FULL", "MBS_MONTE_SEED",
        "MBS_SHARED_GPU_EXPLICIT_SCALE", "MBS_SHARED_PIPELINE",
        "MBS_TRACE_CPU_PREFILTER_MARGIN", "MBS_TRACE_MIN_GAMMA"
    };
    return std::find(std::begin(exact), std::end(exact), name)
        != std::end(exact);
}

} // namespace

RuntimeResourceSnapshot QueryRuntimeResourceSnapshot()
{
    RuntimeResourceSnapshot snapshot;
    char hostname[256] = {};
    if (gethostname(hostname, sizeof(hostname) - 1) == 0)
        snapshot.hostname = hostname;
    else
        snapshot.hostname = "unknown";
    snapshot.cpuModel = ReadCpuModel();
    const long online = sysconf(_SC_NPROCESSORS_ONLN);
    snapshot.logicalProcessors = online > 0
        ? (int)online : (int)std::thread::hardware_concurrency();
    snapshot.physicalCores = ReadPhysicalCoreCount();
    if (snapshot.physicalCores <= 0)
        snapshot.physicalCores = snapshot.logicalProcessors;
#ifdef _OPENMP
    snapshot.openmpMaxThreads = omp_get_max_threads();
#else
    snapshot.openmpMaxThreads = 1;
#endif
    snapshot.ramTotalKb = ReadNamedKb("/proc/meminfo", "MemTotal:");
    snapshot.ramAvailableKb = ReadNamedKb("/proc/meminfo", "MemAvailable:");
    snapshot.processRssKb = ReadNamedKb("/proc/self/status", "VmRSS:");
    snapshot.processPeakRssKb = ReadNamedKb("/proc/self/status", "VmHWM:");
    return snapshot;
}

RuntimeEnvironmentSnapshot QueryRuntimeEnvironmentSnapshot()
{
    RuntimeEnvironmentSnapshot snapshot;
    for (char **entry = environ; entry && *entry; ++entry)
    {
        const std::string item(*entry);
        const size_t equals = item.find('=');
        const std::string name = item.substr(0, equals);
        if (!StartsWith(name, "MBS_"))
            continue;
        const std::string value = equals == std::string::npos
            ? std::string() : item.substr(equals + 1);
        snapshot.active.push_back(std::make_pair(name, value));
        if (IsExperimentalEnvironmentName(name))
            snapshot.experimental.push_back(name);
    }
    std::sort(snapshot.active.begin(), snapshot.active.end());
    std::sort(snapshot.experimental.begin(), snapshot.experimental.end());
    snapshot.experimental.erase(
        std::unique(snapshot.experimental.begin(), snapshot.experimental.end()),
        snapshot.experimental.end());
    return snapshot;
}

std::string FormatRuntimeResourceReport(const std::string &stage,
                                        bool includeActiveGpu)
{
    const RuntimeResourceSnapshot snapshot = QueryRuntimeResourceSnapshot();
    std::ostringstream out;
    out << "\n===== RESOURCE REPORT: " << stage << " =====\n";
    out << "Host: " << snapshot.hostname << "\n";
    out << "CPU: " << snapshot.cpuModel << "\n";
    out << "CPU processors: " << snapshot.logicalProcessors
        << " logical, " << snapshot.physicalCores << " physical cores\n";
    out << "OpenMP max threads: " << snapshot.openmpMaxThreads << "\n";
    out << std::fixed << std::setprecision(1);
    out << "RAM: " << ToMib(snapshot.ramAvailableKb) << " MiB available / "
        << ToMib(snapshot.ramTotalKb) << " MiB total\n";
    out << "Process memory: RSS " << ToMib(snapshot.processRssKb)
        << " MiB, peak RSS " << ToMib(snapshot.processPeakRssKb) << " MiB\n";
    if (includeActiveGpu)
    {
        long long freeBytes = 0;
        long long totalBytes = 0;
        std::string error;
        if (QueryActiveGpuMemory(freeBytes, totalBytes, error))
        {
            const double mib = 1024.0 * 1024.0;
            out << "Active GPU VRAM: " << (totalBytes - freeBytes) / mib
                << " MiB used, " << freeBytes / mib << " MiB free, "
                << totalBytes / mib << " MiB total\n";
        }
        else
        {
            out << "Active GPU VRAM: unavailable (" << error << ")\n";
        }
    }
    out << "=======================================\n";
    return out.str();
}

std::string FormatRuntimeEnvironmentReport(
    const RuntimeEnvironmentSnapshot &snapshot)
{
    std::ostringstream out;
    out << "\n===== MBS ENVIRONMENT =====\n";
    if (snapshot.active.empty())
    {
        out << "Active MBS_* variables: none\n";
    }
    else
    {
        out << "Active MBS_* variables:\n";
        for (const auto &entry : snapshot.active)
            out << "  " << entry.first << '=' << entry.second << '\n';
    }
    if (!snapshot.experimental.empty())
    {
        out << "Explicitly allowed experimental variables:";
        for (const std::string &name : snapshot.experimental)
            out << ' ' << name;
        out << '\n';
    }
    out << "===========================\n";
    return out.str();
}
