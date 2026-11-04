#include "FullAutoRoundPath.h"
#include <cassert>
#include <fstream>
#include <unistd.h>

int main() {
    char pattern[] = "/tmp/mbs-round-path-XXXXXX";
    const char* directory = mkdtemp(pattern);
    assert(directory);
    const std::string base(directory), root = base + "/new", old = base + "/old";
    const std::string snapshots = root + "/round_snapshots", round = snapshots + "/round1";
    assert(mkdir(root.c_str(), 0700) == 0);
    assert(mkdir(snapshots.c_str(), 0700) == 0);
    assert(mkdir(round.c_str(), 0700) == 0);
    assert(symlink(root.c_str(), old.c_str()) == 0);
    std::string resolved;
    assert(FullAutoRoundPath::Resolve(root, round, resolved) && resolved == round);
    assert(FullAutoRoundPath::Resolve(root, old + "/round_snapshots/round1", resolved));
    assert(resolved == round);
    assert(FullAutoRoundPath::Resolve(old, round, resolved) && resolved == round);
    // The production call passes the same string as input and output.
    resolved = old + "/round_snapshots/round1";
    assert(FullAutoRoundPath::Resolve(root, resolved, resolved) && resolved == round);
    const std::string outside = base + "/outside", escape = snapshots + "/escape";
    assert(mkdir(outside.c_str(), 0700) == 0);
    assert(symlink(outside.c_str(), escape.c_str()) == 0);
    assert(!FullAutoRoundPath::Resolve(root, escape, resolved));
    assert(!FullAutoRoundPath::Resolve(root, outside, resolved));
    assert(!FullAutoRoundPath::Resolve(root, snapshots, resolved));
    assert(!FullAutoRoundPath::Resolve(root, snapshots + "/missing", resolved));
    assert(!FullAutoRoundPath::Resolve(root, snapshots + "/../round_snapshots/round1", resolved));
    assert(!FullAutoRoundPath::Resolve(root, "round_snapshots/round1", resolved));
    const std::string file = snapshots + "/file";
    std::ofstream(file.c_str()) << "not a snapshot directory";
    assert(!FullAutoRoundPath::Resolve(root, file, resolved));
    assert(unlink(file.c_str()) == 0);
    assert(unlink(escape.c_str()) == 0);
    assert(unlink(old.c_str()) == 0);
    assert(rmdir(round.c_str()) == 0);
    assert(rmdir(snapshots.c_str()) == 0);
    assert(rmdir(root.c_str()) == 0);
    assert(rmdir(outside.c_str()) == 0);
    assert(rmdir(base.c_str()) == 0);
}
