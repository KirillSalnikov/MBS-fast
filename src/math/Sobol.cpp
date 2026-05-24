#include "Sobol.h"

Sobol2D::Sobol2D(uint32_t seed)
    : m_index(0), m_seed(seed)
{
    m_state[0] = 0;
    m_state[1] = 0;

    initDirections();
}

void Sobol2D::reset()
{
    m_index = 0;
    m_state[0] = 0;
    m_state[1] = 0;
}

uint32_t Sobol2D::hash(uint32_t seed, uint32_t dim, uint32_t level,
                       uint32_t prefix) const
{
    uint32_t h = seed ^ 0x9e3779b9u;
    h ^= dim + 0x85ebca6bu + (h << 6) + (h >> 2);
    h ^= level + 0xc2b2ae35u + (h << 6) + (h >> 2);
    h ^= prefix + 0x27d4eb2fu + (h << 6) + (h >> 2);
    h ^= h >> 16;
    h *= 0x7feb352du;
    h ^= h >> 15;
    h *= 0x846ca68bu;
    h ^= h >> 16;
    return h;
}

uint32_t Sobol2D::scrambleOwen(uint32_t value, uint32_t dim) const
{
    uint32_t out = 0;
    uint32_t prefix = 0;

    for (uint32_t level = 0; level < 32; ++level)
    {
        uint32_t shift = 31u - level;
        uint32_t bit = (value >> shift) & 1u;

        // Base-2 nested Owen scrambling: at every node of the binary
        // prefix tree, choose a seed-dependent permutation of the next bit.
        uint32_t flip = hash(m_seed, dim, level, prefix) & 1u;
        uint32_t scrambledBit = bit ^ flip;

        out |= scrambledBit << shift;
        prefix = (prefix << 1) | bit;
    }

    return out;
}

void Sobol2D::initDirections()
{
    // Dimension 1: Van der Corput (base 2)
    // direction[0][i] = 1 << (31 - i)
    for (int i = 0; i < 32; ++i)
        m_directions[0][i] = 1u << (31 - i);

    // Dimension 2: Joe-Kuo direction numbers for s=1 (primitive polynomial x+1, degree 1)
    // Initial direction number: m_1 = 1
    // Recurrence for degree-1 primitive polynomial p(x) = x + 1:
    //   v[i] = v[i-1] ^ (v[i-1] >> 1)
    // But stored as v[i] * 2^(32-i) (left-justified)
    m_directions[1][0] = 1u << 31; // v_1 = 1, shifted
    for (int i = 1; i < 32; ++i)
    {
        m_directions[1][i] = m_directions[1][i-1] ^ (m_directions[1][i-1] >> 1);
    }
}

void Sobol2D::next(double &x, double &y)
{
    // Gray code: find rightmost zero bit of m_index
    uint32_t c = 0;
    uint32_t val = m_index;
    while ((val & 1) != 0)
    {
        val >>= 1;
        ++c;
    }
    if (c >= 32) c = 0; // safety

    m_state[0] ^= m_directions[0][c];
    m_state[1] ^= m_directions[1][c];
    ++m_index;

    // Apply nested Owen scrambling and convert to [0,1).
    uint32_t s0 = scrambleOwen(m_state[0], 0);
    uint32_t s1 = scrambleOwen(m_state[1], 1);

    x = (double)s0 / 4294967296.0; // 2^32
    y = (double)s1 / 4294967296.0;
}

void Sobol2D::generate(int n, std::vector<double> &x, std::vector<double> &y)
{
    x.resize(n);
    y.resize(n);
    reset();
    for (int i = 0; i < n; ++i)
        next(x[i], y[i]);
}
