#pragma once
#include "AnalyticBackscatter.h"
#include <vector>

namespace AnalyticAzimuthGaussian
{
// A physical circular Gaussian aperture surrogate for one retained internal
// beam. Its conditional laboratory-spin mean is exact; the native polygon,
// polarization, phases and cross terms remain in the full coherent residual.
struct Beam
{
    double dx,dy,dz,peak,kappa,theta,sine;
};
struct Values { double point,mean; };
// This alpha-independent window preserves the conditional-mean identity.
constexpr double PolarExponentLimit=64;
double ScaledI0(double x);
Beam Aperture(const AnalyticBackscatter::Vec &direction,double projectedArea,
              double jonesSquaredNorm,double wavelength);
Values Evaluate(const Beam &beam,double theta,double azimuth);
bool EvaluateCells(const std::vector<std::vector<Beam>>&poses,
    const std::vector<double>&weights,const std::vector<double>&theta,int phi,
    std::vector<double>&pointCells,std::vector<double>&meanCells,
    std::vector<Values>*samples);
}

bool EvaluateAnalyticAzimuthGaussianGpu(
    const std::vector<std::vector<AnalyticAzimuthGaussian::Beam>>&poses,
    const std::vector<double>&weights,const std::vector<double>&theta,int phi,
    std::vector<double>&pointCells,std::vector<double>&meanCells,
    std::vector<AnalyticAzimuthGaussian::Values>*samples);
