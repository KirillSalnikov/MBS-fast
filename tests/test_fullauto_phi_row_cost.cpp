#include "FullAutoPhiRowCost.h"
#include <cassert>
int main() {
    using FullAutoPhiRowCost::Project;
    assert(Project(10,2,9,8,8)==80);
    assert(Project(10,2,9,8,2)==26);
    assert(Project(10,9,9,8,2)==72);
    assert(Project(10,9,0,8,2)==80);
    assert(Project(10,0,9,8,2)==80);
}
