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
