#pragma once

#include "HandlerPO.h"
#include "PoleMueller.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace AdaptivePhi {

inline double Parameter(const char *name, double fallback)
{
    const char *value = std::getenv(name);
    if (!value) return fallback;
    char *end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || *end || !std::isfinite(parsed) || parsed <= 0)
        throw std::runtime_error(std::string("Invalid adaptive phi parameter: ") + name);
    return parsed;
}

struct Config {
    int first = 75;
    int maximum = 9600;
    double intensity = 0.005;
    double polarization = 0.0025;
    explicit Config(bool enabled = true)
    {
        if (!enabled) return;
        const double lo = Parameter("MBS_PHI_ADAPT_MIN", first);
        const double hi = Parameter("MBS_PHI_ADAPT_MAX", maximum);
        if (lo < 16 || hi > 65536 || hi < 4*lo || lo != std::floor(lo)
            || hi != std::floor(hi))
            throw std::runtime_error("Adaptive phi requires integer 16 <= MIN, 4*MIN <= MAX <= 65536");
        first = (int)lo;
        maximum = (int)hi;
        int n = first;
        while (n < maximum) n *= 2;
        if (n != maximum)
            throw std::runtime_error("Adaptive phi MAX must be MIN times a power of two");
        intensity = Parameter("MBS_PHI_ADAPT_M11_TOL", intensity);
        polarization = Parameter("MBS_PHI_ADAPT_POL_TOL", polarization);
    }
};

inline matrix Mean(const Arr2D &samples, int nPhi, int t, double theta)
{
    matrix out(4, 4);
    out.Fill(0.0);
    const bool forward = PoleMueller::IsForward(theta);
    const bool backward = PoleMueller::IsBackward(theta);
    for (int p = 0; p < nPhi; ++p)
    {
        const double *m = samples.RawCell(p, t);
        const double az = -p * (M_2PI / nPhi);
        const double cs = std::cos(2*az), sn = std::sin(2*az);
        for (int r = 0; r < 4; ++r)
        {
            out[r][0] += m[4*r];
            out[r][3] += m[4*r+3];
            out[r][1] += forward || backward ? m[4*r+1] : m[4*r+1]*cs - m[4*r+2]*sn;
            out[r][2] += forward || backward ? m[4*r+2] : m[4*r+1]*sn + m[4*r+2]*cs;
        }
    }
    out /= nPhi;
    if (forward) PoleMueller::ApplyForward(out);
    if (backward) PoleMueller::ApplyBackward(out);
    return out;
}

inline bool Close(const matrix &a, const matrix &b, const Config &cfg,
                  double &intensityError, double &polarizationError)
{
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
            if (!std::isfinite(a[r][c]) || !std::isfinite(b[r][c]))
                throw std::runtime_error("Adaptive phi encountered a non-finite Mueller value");
    const double scale = std::max(std::fabs(b[0][0]), 1e-300);
    intensityError = std::fabs(a[0][0] - b[0][0]) / scale;
    polarizationError = 0;
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
            polarizationError = std::max(polarizationError,
                std::fabs(a[r][c]/std::max(std::fabs(a[0][0]), 1e-300) - b[r][c]/scale));
    return intensityError <= cfg.intensity && polarizationError <= cfg.polarization;
}

// The prepared orientations are reused at every level: no extra ray tracing.
// Acceptance is local to each size/theta/orientation chunk, never based on a
// forward-dominated global norm. Two consecutive passing comparisons are needed.
inline std::vector<matrix> Compute(HandlerPO &handler, const Light &light,
    const std::vector<PreparedOrientation> &prepared, int start, int count,
    double scale, double waveIndex, const Config &cfg, std::ostream &report,
    size_t sizeIndex, int chunkStart)
{
    struct Restore {
        HandlerPO &handler;
        ScatteringRange sphere;
        bool fft;
        explicit Restore(HandlerPO &h) : handler(h), sphere(h.m_sphere), fft(h.IsFftEnabled()) {}
        ~Restore() { std::swap(handler.m_sphere, sphere); handler.SetFftEnabled(fft); }
    } saved(handler);
    handler.SetFftEnabled(false);
    const int rows = saved.sphere.nZenith + 1;
    std::vector<matrix> previous(rows, matrix(4,4)), result(rows, matrix(4,4));
    std::vector<int> streak(rows, 0), active;
    for (int t = 0; t < rows; ++t) active.push_back(t);
    for (int nPhi = cfg.first; !active.empty(); nPhi *= 2)
    {
        ScatteringRange grid = saved.sphere;
        grid.nAzimuth = nPhi;
        grid.azinuthStep = M_2PI / nPhi;
        grid.isNonUniform = true;
        grid.thetaValues.clear();
        for (int t : active) grid.thetaValues.push_back(saved.sphere.GetZenith(t));
        grid.nZenith = (int)active.size() - 1;
        grid.zenithStart = grid.thetaValues.front();
        grid.zenithEnd = grid.thetaValues.back();
        grid.zenithStep = grid.nZenith > 0 ? (grid.zenithEnd-grid.zenithStart)/grid.nZenith : 0;
        grid.ComputeSphereDirections(light);
        handler.m_sphere = grid;
        std::vector<Arr2D> local;
        local.push_back(Arr2D(nPhi+1, grid.nZenith+1, 4, 4));
        local[0].ClearArr();
        if (!handler.HandleOrientationsToLocalGpuMultiK(prepared, start, count,
                                                       std::vector<double>(1,scale), waveIndex, local))
            throw std::runtime_error("Adaptive phi CUDA diffraction failed; refusing unchecked fallback");
        std::vector<int> next;
        for (size_t j = 0; j < active.size(); ++j)
        {
            const int t = active[j];
            matrix mean = Mean(local[0], nPhi, (int)j, saved.sphere.GetZenith(t));
            double ei = 0, ep = 0;
            const bool pass = nPhi > cfg.first && Close(previous[t], mean, cfg, ei, ep);
            streak[t] = pass ? streak[t]+1 : 0;
            previous[t] = mean;
            if (streak[t] >= 2)
            {
                result[t] = mean;
                report << sizeIndex << ',' << chunkStart << ',' << start << ','
                       << RadToDeg(saved.sphere.GetZenith(t)) << ',' << nPhi << ',' << ei << ',' << ep << ",converged\n";
            }
            else
            {
                next.push_back(t);
                if (nPhi == cfg.maximum)
                    report << sizeIndex << ',' << chunkStart << ',' << start << ','
                           << RadToDeg(saved.sphere.GetZenith(t)) << ',' << nPhi << ',' << ei << ',' << ep << ",limit\n";
            }
        }
        std::cerr << "Adaptive phi size=" << sizeIndex << " chunk=" << chunkStart
                  << " Nphi=" << nPhi << " theta=" << active.size()
                  << " remaining=" << next.size() << '\n';
        if (!next.empty() && nPhi == cfg.maximum)
        {
            report.flush();
            throw std::runtime_error("Adaptive phi did not converge at MAX; increase MBS_PHI_ADAPT_MAX");
        }
        active.swap(next);
    }
    report.flush();
    if (!report) throw std::runtime_error("Cannot write adaptive phi convergence report");
    return result;
}
} // namespace AdaptivePhi
