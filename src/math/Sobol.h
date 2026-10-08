#pragma once

#include <cstdint>
#include <cmath>
#include <vector>

/**
 * Simple 2D Sobol quasi-random sequence generator.
 * Uses Van der Corput for dimension 1, Joe-Kuo direction numbers for dimension 2.
 * Supports nested Owen scrambling via seed-dependent per-prefix bit permutations.
 */
class Sobol2D
{
public:
    Sobol2D(uint32_t seed = 0);

    /// Reset state to beginning of sequence
    void reset();

    /// Generate next point in [0,1]^2
    void next(double &x, double &y);

    /// Generate n points in [0,1]^2
    void generate(int n, std::vector<double> &x, std::vector<double> &y);

private:
    uint32_t m_index;
    uint32_t m_seed;
    uint32_t m_state[2];
    uint32_t m_directions[2][32];

    void initDirections();
    uint32_t hash(uint32_t seed, uint32_t dim, uint32_t level,
                  uint32_t prefix) const;
    uint32_t scrambleOwen(uint32_t value, uint32_t dim) const;
};

// Three-dimensional Sobol sequence with the same first two dimensions and
// nested Owen permutations as Sobol2D. Dimension 3: x^2+x+1, m=(1,3).
class Sobol3D
{
public:
    explicit Sobol3D(uint32_t seed=0, bool scramble=true);
    void next(double &x, double &y, double &z);
private:
    uint32_t index, seed, state[3], directions[3][32];
    bool scramble;
    uint32_t Scramble(uint32_t value, uint32_t dim) const;
};
