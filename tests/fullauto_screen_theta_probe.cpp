// Exercise experiment selection through the same parser and initialization as fullauto.
#include "../src/FullAuto.cpp"

int main(int argc, const char* argv[]) {
    try {
        using namespace FullAuto;
        Require(argc == 2, "expected private output directory");
        const std::string root = argv[1];
        const std::vector<std::string> base = {"probe", "--method", "po", "--backend", "cpu",
            "--particle", "1", "2", "1", "--ri", "1.3116", "0", "--wavelength-um", ".532",
            "--threads", "1", "--fullauto", ".04", "--output", root, "--allow-experimental-environment"};
        auto check = [&](std::vector<std::string> extra, bool enabled) {
            auto command = base;
            command.insert(command.end(), extra.begin(), extra.end());
            std::vector<const char*> raw;
            for (const auto& arg : command) raw.push_back(arg.c_str());
            ArgPP args;
            RegisterCliOptions(args);
            args.Parse(raw.size(), raw.data());
            auto config = RunConfig::FromCommandLine(args, false, false);
            Engine controller(args, config);
            controller.Init(raw.size(), raw.data());
            Require(controller.screenTheta == enabled, "unexpected theta pilot selection");
        };
        unsetenv("MBS_FULLAUTO_SCREEN_THETA");
        check({}, false);
        setenv("MBS_FULLAUTO_SCREEN_THETA", "1", 1);
        check({}, true);
        check({"--fullauto-theta-audit-count", "64"}, false);
        check({"--fullauto-theta-start", "33"}, false);
        setenv("MBS_FULLAUTO_CARDINAL_THETA", "auto", 1);
        check({}, false);

        const std::vector<double> grid{0, 30, 180};
        Require(FullAutoSeedKnots::Parse("0 30.00000000000001 180", grid).size() == 3,
                "rejected roundoff in a grid angle");
        for (const char* bad : {"", "30 180", "0 30", "0 30 30 180", "0 nan 180",
                                "0 30.01 180", "0 180 x", "0 181"}) {
            bool rejected = false;
            try { FullAutoSeedKnots::Parse(bad, grid); }
            catch (const std::runtime_error&) { rejected = true; }
            Require(rejected, "accepted invalid angular proposal");
        }
        std::cout << "theta pilot CLI precedence and proposal validation passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
