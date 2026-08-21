#include "Mueller.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

int main()
{
    double worst = 0.0;
    for (int sample = 0; sample < 64; ++sample)
    {
        const double x = 0.125 * (sample + 1);
        const complex j00(std::sin(x), std::cos(0.7*x));
        const complex j01(std::cos(1.3*x), -std::sin(0.4*x));
        const complex j10(-std::sin(0.9*x), std::cos(1.1*x));
        const complex j11(std::cos(0.6*x), std::sin(1.7*x));
        const double weight = 0.03125 * (sample + 1);

        matrixC jones(2, 2);
        jones[0][0] = j00;
        jones[0][1] = j01;
        jones[1][0] = j10;
        jones[1][1] = j11;
        matrix reference = Mueller(jones);
        reference *= weight;

        double direct[16] = {};
        AddMuellerFromJones(direct, j00, j01, j10, j11, weight);
        for (int row = 0; row < 4; ++row)
        {
            for (int column = 0; column < 4; ++column)
            {
                const double error = std::fabs(
                    direct[4*row + column] - reference[row][column]);
                worst = std::max(worst, error);
            }
        }
    }

    if (worst > 2e-15)
    {
        std::cerr << "direct Jones-to-Mueller mismatch: " << worst << '\n';
        return 1;
    }
    std::cout << "direct Jones-to-Mueller regression passed: max_abs="
              << worst << '\n';
    return 0;
}
