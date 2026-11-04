#ifndef MBS_FULLAUTO_ROUND_PATH_H
#define MBS_FULLAUTO_ROUND_PATH_H

#include <cstdlib>
#include <string>
#include <sys/stat.h>

namespace FullAutoRoundPath {

// A restored checkpoint may retain an old pathname that now aliases its new
// filesystem. Accept the alias only when the actual directory remains below
// this result's round_snapshots directory. Do not accept symlink escapes.
inline bool Resolve(const std::string& root, const std::string& candidate,
                    std::string& resolved) {
    if (candidate.empty() || candidate[0] != '/'
            || candidate.find('\\') != std::string::npos
            || candidate.find("..") != std::string::npos) return false;
    char* canonicalRoot = realpath(root.c_str(), nullptr);
    if (!canonicalRoot) return false;
    const std::string prefix = std::string(canonicalRoot) + "/round_snapshots/";
    std::free(canonicalRoot);
    char* canonicalCandidate = realpath(candidate.c_str(), nullptr);
    if (!canonicalCandidate) return false;
    const std::string actual(canonicalCandidate);
    std::free(canonicalCandidate);
    struct stat info;
    if (actual.compare(0, prefix.size(), prefix) != 0
            || stat(actual.c_str(), &info) != 0 || !S_ISDIR(info.st_mode)) return false;
    resolved = actual;
    return true;
}

} // namespace FullAutoRoundPath
#endif
