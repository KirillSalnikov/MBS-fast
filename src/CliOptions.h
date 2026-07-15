#pragma once

#include <initializer_list>
#include <string>
#include <vector>

class ArgPP;

struct CliOptionSpec
{
    CliOptionSpec(const char *keyValue,
                  const char *canonicalValue,
                  char valueCountValue,
                  const char *valueSyntaxValue,
                  const char *sectionValue,
                  const char *descriptionValue,
                  bool debugOnlyValue = false,
                  std::initializer_list<const char *> extraAliases = {});

    std::string key;
    std::string canonical;
    char valueCount;
    std::string valueSyntax;
    std::string section;
    std::string description;
    bool debugOnly;
    std::vector<std::string> aliases;
};

const std::vector<CliOptionSpec> &GetCliOptionSpecs();
const CliOptionSpec *FindCliOptionSpec(const std::string &keyOrAlias);
std::string CliCanonicalFlag(const std::string &keyOrAlias);
void RegisterCliOptions(ArgPP &parser);
void PrintGeneratedHelp(bool includeDebugOptions);

