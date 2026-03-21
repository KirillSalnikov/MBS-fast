#include "Sobol.h"

Sobol2D::Sobol2D(uint32_t seed)
    : m_index(0)
{
    m_state[0] = 0;
    m_state[1] = 0;

    initDirections();

    // Owen scrambling: hash seed into per-dimension scramble values
    m_scramble[0] = hash(seed, 0);
    m_scramble[1] = hash(seed, 1);
}

void Sobol2D::reset()
{
    m_index = 0;
    m_state[0] = 0;
    m_state[1] = 0;
}

uint32_t Sobol2D::hash(uint32_t seed, uint32_t dim)
{
    // Simple hash for scrambling
    uint32_t h = seed * 2654435761u + dim * 2246822519u;
    h ^= h >> 16;
    h *= 0x45d9f3b;
    h ^= h >> 16;
    return h;
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

    // Apply scramble and convert to [0,1)
    uint32_t s0 = m_state[0] ^ m_scramble[0];
    uint32_t s1 = m_state[1] ^ m_scramble[1];

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
