#pragma once
#include "AnalyticBackscatter.h"
#include <array>
#include <string>

struct AnalyticGpuModel;

namespace AnalyticFacetAverage
{
using AnalyticBackscatter::Vec;
using AnalyticBackscatter::Face;
struct Components { double reflection,shadow; };
// Coefficients must be frozen from a separate pilot, never fitted to the
// production samples. The file has theta_deg/reflection_weight/shadow_weight.
std::vector<std::array<double,2>> ReadControlWeights(
    const std::string &path,const std::vector<double> &thetaRadians);
// Unit-vector dot products within eight FP64 epsilons of tangency have
// ambiguous sign after a rigid rotation. Use the same boundary on CPU/CUDA.
constexpr double GrazingTolerance=1.7763568394002505e-15;

// Body-frame directions. Source points TOWARD the source, opposite to the
// incident wave's propagation. theta=0 observes along -source (forward).
// The azimuth reference is projected onto the plane transverse to source.
std::array<Vec,3> SourceFrame(const Vec &source,const Vec &azimuthReference);
Vec Observer(const std::array<Vec,3> &frame,double theta,double azimuth);

// Normalized SO(3) means of physical reflected-facet self intensities and
// projected shadow apertures. Fresnel moments need one smooth 1D quadrature;
// angular spin and all oscillatory orientation phases are integrated by the
// Legendre addition theorem and Rayleigh expansion. Full MBS-C keeps all cross
// terms and clipping/refraction corrections.
class Control
{
public:
    Control(const std::vector<Face>&faces,std::complex<double> index,
            double wavelength,const std::vector<double>&theta,
            const std::string &shadowMode="facets",const std::string &meanCache="",
            const std::string &meanReferenceGrid="");
    // Arbitrary incidence and observation; nonzero finite vectors are
    // normalized and the scattering geometry is derived from their pair.
    Components Evaluate(const Vec&source,const Vec&observer) const;
    // Compatibility overload: rejects an angle inconsistent with directions.
    Components Evaluate(const Vec&source,const Vec&observer,double theta) const;
    Components Mean(int row) const { return means.at(row); }
    double RefinementError() const { return refinementError; }
    bool SupportsDomain(double betaSym,double gammaSym) const;
    AnalyticGpuModel GpuModel() const;
    static std::vector<double> ReflectionMoments(std::complex<double> index,int degree,
                                                double theta);
private:
    struct Patch
    {
        Vec n,h,v;
        std::vector<double> x,y;
        double area,diameter;
    };
    std::vector<Patch> patches;
    std::vector<Face> faces;
    std::vector<Components> means;
    std::complex<double> index;
    double wave,refinementError,forwardShadowMean,circleRadius;
    std::string shadowMode;
    Components EvaluateUnit(const Vec&source,const Vec&observer) const;
    void PrepareMeans(const std::vector<double>&theta);
    void PrepareCachedMeans(const std::vector<double>&theta,const std::string &path);
};
}
