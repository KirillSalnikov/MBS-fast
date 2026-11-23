#include "../src/scattering/GpuTraceBeamBudget.h"

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <future>
#include <iostream>

int main()
{
    using Reservation = GpuTraceBeamBudget::Reservation;
    using namespace std::chrono;
    GpuTraceBeamBudget budget(128);
    std::future<bool> small;
    {
        Reservation holder(budget, 128);
        assert(holder.Reserved());
        auto large = std::async(std::launch::async, [&] {
            return Reservation(budget, 96).Reserved();
        });
        small = std::async(std::launch::async, [&] {
            return Reservation(budget, 32).Reserved();
        });
        assert(large.wait_for(milliseconds(25)) == std::future_status::timeout);
        assert(budget.ReduceAfterOom() == 64);
        if (large.wait_for(seconds(2)) != std::future_status::ready)
        {
            std::cerr << "Oversized waiter did not wake after OOM\n";
            std::abort();
        }
        assert(!large.get());
        assert(budget.ReduceAfterOom() == 32);
        assert(!Reservation(budget, 33).Reserved());
        assert(small.wait_for(milliseconds(25)) == std::future_status::timeout);
    }
    assert(small.get());
    assert(Reservation(budget, 32).Reserved());
}
