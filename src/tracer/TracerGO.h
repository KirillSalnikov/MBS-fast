#pragma once

#include "Tracer.h"

class TracerGO : public Tracer
{
public:
	TracerGO(Particle *particle, int reflNum, const std::string &resultFileName);

	void TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange);
	void TraceFixed(const double &beta, const double &gamma);
    void TraceMonteCarlo(const AngleRange &betaRange, const AngleRange &gammaRange,
                         int nOrientations);
    void TraceSobol(int nOrientations, unsigned int seed,
                    double betaMax, double gammaMax, unsigned int offset=0);

protected:
	double CalcNorm(long long orNum);
    void OutputSummary(int orNumber, CalcTimer &timer);
};
