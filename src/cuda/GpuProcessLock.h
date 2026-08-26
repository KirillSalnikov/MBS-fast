#pragma once

class GpuProcessLock
{
public:
    explicit GpuProcessLock(int device = -1);
    ~GpuProcessLock();

    GpuProcessLock(const GpuProcessLock &) = delete;
    GpuProcessLock &operator=(const GpuProcessLock &) = delete;

private:
    void *m_state;
};
