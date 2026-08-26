#include "bigint/BigInteger.hh"

#include <cassert>
#include <climits>
#include <initializer_list>

int main()
{
    const unsigned long values[] = {0, 1, 2, 17, ULONG_MAX};
    for (unsigned long value : values)
    {
        for (unsigned long multiplier : values)
        {
            for (unsigned long addend : values)
            {
                BigUnsigned expected(value);
                expected = expected * BigUnsigned(multiplier)
                         + BigUnsigned(addend);
                BigUnsigned actual(value);
                actual.multiplyAdd(multiplier, addend);
                assert(actual == expected);
            }
        }
    }

    for (unsigned long base : {2UL, 17UL, 257UL})
    {
        BigInteger expected(0);
        BigInteger actual(0);
        for (int step = 0; step < 512; ++step)
        {
            const unsigned long digit =
                static_cast<unsigned long>(step % (base - 1) + 1);
            expected = (expected + BigInteger(digit)) * BigInteger(base);
            actual.multiplyAddNonnegative(base, digit*base);
            assert(actual == expected);
        }
    }
    return 0;
}
