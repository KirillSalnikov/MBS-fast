#include "GpuProcessLock.h"

#ifdef USE_CUDA
#include <cuda_runtime.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <mutex>
#include <string>
#include <sys/file.h>
#include <unistd.h>

namespace
{
struct ProcessLockState
{
    std::mutex mutex;
    int fd = -1;
    unsigned activeCalls = 0;
    bool unavailable = false;
};

ProcessLockState g_processLocks[32];
std::once_flag g_processLockWarning;

bool ProcessLockEnabled()
{
    const char *value = std::getenv("MBS_GPU_TRACE_PROCESS_LOCK");
    return !(value && value[0] == '0' && value[1] == '\0');
}

void WarnProcessLock(const char *reason)
{
    std::call_once(g_processLockWarning, [reason]() {
        std::fprintf(stderr,
            "WARNING: CUDA process lock is unavailable (%s); concurrent "
            "calculations on one physical GPU can exhaust kernel-stack "
            "memory.\n",
            reason);
    });
}

bool LockFile(int fd)
{
    int result;
    do
    {
        result = flock(fd, LOCK_EX);
    }
    while (result != 0 && errno == EINTR);
    return result == 0;
}

void UnlockFile(int fd)
{
    int result;
    do
    {
        result = flock(fd, LOCK_UN);
    }
    while (result != 0 && errno == EINTR);
}
} // namespace
#endif

GpuProcessLock::GpuProcessLock(int device)
    : m_state(nullptr)
{
#ifdef USE_CUDA
    if (!ProcessLockEnabled())
        return;
    if (device < 0 && cudaGetDevice(&device) != cudaSuccess)
    {
        WarnProcessLock("cannot determine the active CUDA device");
        return;
    }
    if (device < 0
        || device >= static_cast<int>(
            sizeof(g_processLocks) / sizeof(g_processLocks[0])))
    {
        WarnProcessLock("CUDA device index is outside the lock table");
        return;
    }

    ProcessLockState &state = g_processLocks[device];
    std::unique_lock<std::mutex> lock(state.mutex);
    if (state.unavailable)
        return;
    if (state.fd < 0)
    {
        char pciBusId[32] = {};
        if (cudaDeviceGetPCIBusId(pciBusId, sizeof(pciBusId), device)
            != cudaSuccess)
        {
            state.unavailable = true;
            WarnProcessLock("cannot determine the CUDA PCI address");
            return;
        }
        std::string suffix(pciBusId);
        for (char &value : suffix)
        {
            if (value == ':' || value == '.')
                value = '_';
        }
        const std::string path = "/tmp/mbs_gpu_trace_" + suffix + ".lock";
        state.fd = open(path.c_str(), O_CREAT | O_RDONLY | O_CLOEXEC, 0666);
        if (state.fd < 0)
        {
            state.unavailable = true;
            WarnProcessLock(std::strerror(errno));
            return;
        }
    }
    if (state.activeCalls == 0 && !LockFile(state.fd))
    {
        state.unavailable = true;
        WarnProcessLock(std::strerror(errno));
        return;
    }
    ++state.activeCalls;
    m_state = &state;
#else
    (void)device;
#endif
}

GpuProcessLock::~GpuProcessLock()
{
#ifdef USE_CUDA
    ProcessLockState *state = static_cast<ProcessLockState *>(m_state);
    if (state == nullptr)
        return;
    std::lock_guard<std::mutex> lock(state->mutex);
    if (--state->activeCalls == 0)
        UnlockFile(state->fd);
#endif
}
