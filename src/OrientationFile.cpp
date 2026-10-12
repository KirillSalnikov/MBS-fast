#include "OrientationFile.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "common/global.h"

namespace
{

std::runtime_error FileError(const std::string &path, int lineNumber,
                             const std::string &message)
{
    std::ostringstream text;
    text << "orientation file '" << path << "'";
    if (lineNumber > 0)
        text << " line " << lineNumber;
    text << ": " << message;
    return std::runtime_error(text.str());
}

} // namespace

std::vector<WeightedEulerOrientationRadians> ReadSO3RuleDegrees(const std::string& path)
{
    std::ifstream input(path.c_str());if(!input)throw FileError(path,0,"cannot open file");
    std::vector<WeightedEulerOrientationRadians> rule;std::string text;int number=0;double sum=0;
    while(std::getline(input,text)){
        ++number;size_t comment=text.find('#');if(comment!=std::string::npos)text.erase(comment);
        std::istringstream row(text);row>>std::ws;if(row.eof())continue;
        double beta,gamma,alpha,weight;std::string extra;
        if(!(row>>beta>>gamma>>alpha>>weight)||(row>>extra))
            throw FileError(path,number,"expected BETA_DEG GAMMA_DEG ALPHA_DEG WEIGHT");
        if(!std::isfinite(beta+gamma+alpha+weight)||beta<0||beta>180||gamma<0||gamma>360||alpha<0||alpha>360||weight<0)
            throw FileError(path,number,"SO(3) angles outside their degree ranges or negative/nonfinite weight");
        rule.push_back({DegToRad(beta),DegToRad(gamma),DegToRad(alpha),weight});sum+=weight;
    }
    if(rule.empty()||std::fabs(sum-1)>1e-10)
        throw FileError(path,0,"SO(3) rule must be nonempty and weights must sum to one");
    return rule;
}

std::vector<EulerOrientationRadians> ReadOrientationFileDegrees(
    const std::string &path)
{
    std::ifstream input(path.c_str());
    if (!input)
        throw FileError(path, 0, "cannot open file");

    std::vector<EulerOrientationRadians> orientations;
    std::string text;
    int lineNumber = 0;
    while (std::getline(input, text))
    {
        ++lineNumber;
        const size_t comment = text.find('#');
        if (comment != std::string::npos)
            text.erase(comment);

        std::istringstream line(text);
        line >> std::ws;
        if (line.eof())
            continue;

        double betaDeg = 0.0;
        double gammaDeg = 0.0;
        if (!(line >> betaDeg >> gammaDeg))
            throw FileError(path, lineNumber,
                            "expected BETA_DEG GAMMA_DEG");
        if (!std::isfinite(betaDeg) || !std::isfinite(gammaDeg))
            throw FileError(path, lineNumber, "angles must be finite");
        if (betaDeg < 0.0 || betaDeg > 180.0)
            throw FileError(path, lineNumber,
                            "beta is outside [0, 180] degrees");

        std::string extra;
        if (line >> extra)
            throw FileError(path, lineNumber,
                            "unexpected token '" + extra + "'");

        EulerOrientationRadians orientation;
        orientation.beta = DegToRad(betaDeg);
        orientation.gamma = DegToRad(gammaDeg);
        orientations.push_back(orientation);
    }

    if (orientations.empty())
        throw FileError(path, 0, "contains no orientations");
    return orientations;
}
