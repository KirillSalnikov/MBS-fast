#include "matrix.hpp"
#include <cassert>
#include <iostream>

// Independent elementwise expectations catch reversed compound subtraction.
// Unequal rectangular matrices also catch square-only indexing assumptions.
int main()
{
    matrix a(2, 3), b(2, 3);
    matrixC ca(2, 3), cb(2, 3);
    for (unsigned i = 0; i < 2; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            a[i][j] = 7 + 3*i - double(j);
            b[i][j] = -2 + double(i) + 2*j;
            ca[i][j] = complex(a[i][j], 4 + i + j);
            cb[i][j] = complex(b[i][j], -3 + double(i) - j);
        }
    const matrix original(a), rhs(b);
    const matrixC coriginal(ca), crhs(cb);
    a -= b;
    ca -= cb;
    for (unsigned i = 0; i < 2; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            assert(a[i][j] == 9 + 2*i - 3*double(j));
            assert(real(ca[i][j]) == a[i][j]);
            assert(imag(ca[i][j]) == 7 + 2*j);
        }
    assert(b == rhs && cb == crhs);
    a += b;
    ca += cb;
    assert(a == original && ca == coriginal);
    a -= a;
    ca -= ca;
    assert(norm(a) == 0 && norm(ca) == 0);

    // A scattered-field accumulator must cancel identical incident fields
    // even when the aperture is split into several patches.
    a = original + rhs;
    ca = coriginal + crhs;
    a -= original;
    ca -= coriginal;
    a -= rhs;
    ca -= crhs;
    assert(norm(a) == 0 && norm(ca) == 0);
    std::cout << "Real/complex rectangular subtraction, aliasing and field cancellation passed\n";
}
