#include "Splitting.h"
#include <cassert>
#include <cmath>
#include <iostream>

int main()
{
    Splitting splitting(true);
    // The lossless effective squared index cannot depend on incidence, even
    // at grazing incidence. Include subnormal values in the angle sweep.
    for (double index : {0.0, 0.5, 1.0, 1.3116, 4.0, 100.0})
    {
        splitting.ComputeRiParams(complex(index, 0));
        for (double cosine : {1.0, 0.3, 1e-12, 1e-200, 0.0, -1e-12, -1.0})
        {
            splitting.cosA = cosine;
            assert(splitting.ComputeEffectiveReRi() == index*index);
        }
    }
    // At normal incidence the absorbing expression also reduces to n_real^2.
    // A nonzero imaginary part must retain its angle-dependent grazing limit.
    for (double loss : {1e-7, 0.00042407, 0.1, 2.0})
    {
        splitting.ComputeRiParams(complex(1.3116, loss));
        splitting.cosA = 1.0;
        const double normal = splitting.ComputeEffectiveReRi();
        assert(std::fabs(normal - 1.3116*1.3116) < 2e-15);
        splitting.cosA = 1e-12;
        const double grazing = splitting.ComputeEffectiveReRi();
        assert(std::isfinite(grazing) && grazing > normal);
    }
    // A purely imaginary index has no propagating real part in this model.
    splitting.ComputeRiParams(complex(0, 0.5));
    splitting.cosA = 0.0;
    assert(splitting.ComputeEffectiveReRi() == 0.0);
    std::cout << "lossless and absorbing splitting limits passed\n";
}
