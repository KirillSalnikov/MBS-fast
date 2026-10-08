#pragma once

#include <complex>
#include <vector>

// Physical return-strip control, not a fitted orientation density. Its frozen
// edge profile has a known spherical mean. Full MBS minus this profile retains
// every traced path and every interference term, even when the strip geometry
// only approximates the real aperture (e.g. a hexagonal base).
namespace AnalyticBackscatter
{
struct Vec
{
    double x, y, z;
    Vec(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
};
struct Face
{
    Vec normal;
    std::vector<Vec> vertices;
};
struct Strip
{
    Vec normalX, normalZ, edge;
    double B, L, H;
};

double FiniteSincIntegral(double a);
std::vector<double> SphericalBessel(int degree, double x);
std::vector<Strip> FindReturnStrips(const std::vector<Face> &faces);
// Coherent sum within one top-entry image family (p round trips normal to
// entry, 2q-1 lateral reflections). Its optical phase cancels in this self term.
double ReturnSelfIntensity(double beta, double index, double B, double L,
                           double H, double wavelength, int p, int q);

class ReturnStrip
{
public:
    ReturnStrip(double index, double B, double L, double H, double wavelength,
                int degree = 96, int returnOrder = 1);
    double Evaluate(double sourceX, double sourceZ, double sourceEdge) const;
    double Mean() const { return mean; }
    double LeadingMean() const { return leadingMean; }
private:
    double n, B, L, H, wave, mean, leadingMean;
    int order;
};

class ReturnControl
{
public:
    ReturnControl(const std::vector<Face> &faces, double index, double wavelength,
                  int returnOrder = 1);
    double Evaluate(const Vec &bodySource) const;
    double Mean() const { return mean; }
    double RefinementError() const { return refinementError; }
    int StripCount() const { return static_cast<int>(strips.size()); }
    // Check the actual strip multiset, rather than assuming user-declared
    // symmetry. A periodic gamma sector and an optional z reflection tile S2.
    bool SupportsDomain(double betaSym, double gammaSym) const;
private:
    std::vector<Strip> strips;
    std::vector<ReturnStrip> models;
    std::vector<int> modelIndex;
    double mean, refinementError;
};
}
