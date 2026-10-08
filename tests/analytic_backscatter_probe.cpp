#include "AnalyticBackscatter.h"
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char **argv)
{
    using namespace AnalyticBackscatter;
    std::cout << std::setprecision(17);
    if (argc==7 && std::string(argv[1])=="hex-control")
    {
        const double pi=std::acos(-1.),height=std::atof(argv[3]),diameter=std::atof(argv[4]);
        std::vector<Face> faces;
        for (double sign : {-1.,1.})
        {
            Face face;face.normal=Vec(0,0,sign);
            for (int i=0; i<6; ++i)
                face.vertices.push_back(Vec(diameter/2*std::cos(i*pi/3),diameter/2*std::sin(i*pi/3),sign*height/2));
            faces.push_back(face);
        }
        for (int i=0; i<6; ++i)
        {
            Face face;face.normal=Vec(std::cos((i+.5)*pi/3),std::sin((i+.5)*pi/3),0);
            for (const auto &pair : {std::make_pair(i,1.),std::make_pair(i+1,1.),
                                     std::make_pair(i+1,-1.),std::make_pair(i,-1.)})
                face.vertices.push_back(Vec(diameter/2*std::cos(pair.first*pi/3),diameter/2*std::sin(pair.first*pi/3),pair.second*height/2));
            faces.push_back(face);
        }
        const auto start=std::chrono::steady_clock::now();
        ReturnControl control(faces,std::atof(argv[2]),std::atof(argv[5]),std::atoi(argv[6]));
        std::cout << control.Mean() << ' ' << control.RefinementError() << ' ' << control.StripCount()
                  << ' ' << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() << '\n';
        double beta,gamma;
        while (std::cin >> beta >> gamma)
        {
            beta*=pi/180;gamma*=pi/180;
            std::cout << control.Evaluate(Vec(-std::sin(beta)*std::cos(gamma),std::sin(beta)*std::sin(gamma),std::cos(beta))) << '\n';
        }
        return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "special")
    {
        for (double x : {0., 1e-8, 1e-5, .001, .1, .99, 1., 1.01, 2., 3., 10., 30., 160., 320., 1000., 1e5})
        {
            std::cout << "sinc " << x << ' ' << FiniteSincIntegral(x) << '\n';
            const auto j=SphericalBessel(160,x);
            for (int l : {0,1,2,16,64,96,160})
                std::cout << "bessel " << x << ' ' << l << ' ' << j[l] << '\n';
        }
        return 0;
    }
    if ((argc == 10 || argc == 11) && std::string(argv[1]) == "strip")
    {
        ReturnStrip strip(std::atof(argv[2]),std::atof(argv[3]),std::atof(argv[4]),
                          std::atof(argv[5]),std::atof(argv[6]),160,argc==11 ? std::atoi(argv[10]) : 1);
        std::cout << strip.Mean() << ' ' << strip.LeadingMean() << ' '
                  << strip.Evaluate(std::atof(argv[7]),std::atof(argv[8]),std::atof(argv[9])) << '\n';
        return 0;
    }
    if (argc==10 && std::string(argv[1])=="self")
    {
        std::cout << ReturnSelfIntensity(std::atof(argv[2]),std::atof(argv[3]),std::atof(argv[4]),
            std::atof(argv[5]),std::atof(argv[6]),std::atof(argv[7]),std::atoi(argv[8]),std::atoi(argv[9])) << '\n';
        return 0;
    }
    std::cerr << "usage: probe special | strip N B L H W SX SZ SE\n";
    return 2;
}
