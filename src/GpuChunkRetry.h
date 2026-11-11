#pragma once
// A wide resident GPU packet may be legal on the device but unsafe to prepare
// on the host after a device-path refusal. Abort that native request before
// host preparation; fullauto may replay its identical sampling once, smaller.
#include <algorithm>
#include <climits>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace GpuChunkRetry {
constexpr int SafeUpper = 128;
constexpr int StatusCode = 76;

class Required : public std::runtime_error {
public:
    Required(size_t count, const char* reason)
        : std::runtime_error("resident GPU packet of " + std::to_string(count)
            + " orientations cannot safely fall back to CPU preparation: "
            + (reason ? reason : "device path unavailable")
            + ".\n  Fix: retry the identical sampling with --orientation-chunk 128 or smaller.") {}
};

inline int ExitStatus(bool persistentWorker) { return persistentWorker ? StatusCode : 1; }

inline int PositiveUpper(const std::string& value) {
    if (value.empty()) throw std::runtime_error("missing orientation chunk upper bound");
    int result = 0;
    for (char digit : value) {
        if (digit < '0' || digit > '9' || result > (INT_MAX - (digit - '0')) / 10)
            throw std::runtime_error("invalid orientation chunk upper bound");
        result = result * 10 + (digit - '0');
    }
    if (!result) throw std::runtime_error("orientation chunk upper bound must be positive");
    return result;
}

// Commands are the controller's canonical native argv, not shell text. Make a
// complete plan before changing any word, including duplicate explicit limits.
// Only this scheduling upper bound changes: all physics/sampling argv survive.
inline int LimitCommand(std::vector<std::string>& command, int dispatchUpper,
                        int configuredUpper = 0) {
    int upper = SafeUpper;
    if (dispatchUpper > 0) upper = std::min(upper, dispatchUpper);
    if (configuredUpper > 0) upper = std::min(upper, configuredUpper);
    std::vector<size_t> values;
    for (size_t i = 0; i < command.size(); ++i) {
        if (command[i] != "--orientation-chunk") continue;
        if (i + 1 == command.size()) throw std::runtime_error("missing orientation chunk upper bound");
        values.push_back(i + 1);
        upper = std::min(upper, PositiveUpper(command[++i]));
    }
    auto next = command;
    if (values.empty()) {
        next.push_back("--orientation-chunk");
        next.push_back(std::to_string(upper));
    } else {
        for (size_t index : values) next[index] = std::to_string(upper);
    }
    command.swap(next);
    return upper;
}
}
