#include "../src/FullAutoPoleBudget.h"
#include <cassert>

int main() {
    using FullAutoPoleBudget::Count;
    assert(Count(512, 65536, 16, 1024, 0) == 512);
    assert(Count(512, 65536, 16, 30, 2) == 512);
    assert(Count(512, 65536, 1, 1024, 2) == 512);
    assert(Count(512, 65536, 16, 1024, 2) == 4096);
    assert(Count(512, 3000, 16, 1024, 2) == 2048);
    assert(Count(512, 65536, 4, 1024, 1) == 2048);
}
