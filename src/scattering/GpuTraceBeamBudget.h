#pragma once

#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <mutex>

class GpuTraceBeamBudget
{
public:
    explicit GpuTraceBeamBudget(std::size_t limit) : m_limit(limit) {}

    std::size_t ReduceAfterOom()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_limit = std::max<std::size_t>(32, m_limit / 2);
        const auto reduced = m_limit;
        lock.unlock();
        m_available.notify_all();
        return reduced;
    }

    class Reservation
    {
    public:
        Reservation(GpuTraceBeamBudget &budget, std::size_t requested)
            : m_budget(budget), m_count(0)
        {
            std::unique_lock<std::mutex> lock(budget.m_mutex);
            budget.m_available.wait(lock, [&budget, requested]() {
                // Another thread can shrink the limit while this batch waits.
                return requested > budget.m_limit ||
                       requested <= budget.m_limit - std::min(budget.m_active, budget.m_limit);
            });
            if (requested <= budget.m_limit)
            {
                budget.m_active += requested;
                m_count = requested;
            }
        }

        bool Reserved() const { return m_count != 0; }

        ~Reservation()
        {
            if (m_count == 0)
                return;
            {
                std::lock_guard<std::mutex> lock(m_budget.m_mutex);
                m_budget.m_active -= m_count;
            }
            m_budget.m_available.notify_all();
        }

        Reservation(const Reservation &) = delete;
        Reservation &operator=(const Reservation &) = delete;

    private:
        GpuTraceBeamBudget &m_budget;
        std::size_t m_count;
    };

private:
    std::mutex m_mutex;
    std::condition_variable m_available;
    std::size_t m_active = 0;
    std::size_t m_limit;
};
