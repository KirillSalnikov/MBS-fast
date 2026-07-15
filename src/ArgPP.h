#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <typeinfo>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <cerrno>
#include <climits>
#include <cmath>
#include <algorithm>
#include <set>

class ArgPP
{
public: // methods

	ArgPP() {}

	// TODO: равнозначные аргументы (либо один, либо другой)
	// TODO: write 'throw' in method's declarations
	void AddRule(const std::string &key, char valueNum = 0,
				 bool isOptional = false, const std::string &dependsOn = "")
	{
		Rule rule(valueNum, isOptional, dependsOn, key);
		AddArgRule(key, rule);
	}

	void AddRule(const std::string &key, const std::string &longKey, char nValues = 0,
				 bool isOptional = false, const std::string &dependsOn = "")
	{
		if (longKey.size() > 1)
		{
			AddRule(key, nValues, isOptional, dependsOn);
			AddAlias(longKey, key);
		}
		else
		{
			std::string msg = "long key must consists of more than 1 chatacter";
			Error(msg);
		}
	}

	void AddAlias(const std::string &alias, const std::string &key)
	{
		if (alias == key)
			return;
		auto target = m_rules.find(key);
		if (target == m_rules.end())
			Error("cannot add alias '" + alias + "': target argument '" + key + "' is not registered");
		Rule aliasRule = target->second;
		aliasRule.storageKey = key;
		AddArgRule(alias, aliasRule);
	}

	void Parse(int argc, const char *argv[])
	{
		m_args.clear();
		m_sourceKeys.clear();
		m_programName = argv[0];
		std::vector<std::string> rawArgs(argv+1, argv+argc);
		FillArgs(rawArgs);
		CheckRequiredArgs();
		CheckDependencies();
	}

	int GetIntValue(const std::string &key, size_t valueIndex = 0) const
	{
		std::string rawValue = FindArgValue(key, valueIndex);
		return ConvertTo<int>(rawValue);
	}

	double GetDoubleValue(const std::string &key, size_t valueIndex = 0) const
	{
		std::string rawValue = FindArgValue(key, valueIndex);
		return ConvertTo<double>(rawValue);
	}

	std::string GetStringValue(const std::string &key, size_t valueIndex = 0) const
	{
		return FindArgValue(key, valueIndex);
	}

	const std::string &GetProgramName() const
	{
		return m_programName;
	}

	bool IsCatched(const std::string &key) const
	{
		auto it = m_args.find(NormalizeKey(key));
		return it != m_args.end();
	}

	unsigned GetArgNumber(const std::string &key) const
	{
		unsigned num = 0;

		auto it = m_args.find(NormalizeKey(key));

		if (it != m_args.end())
		{
			num = it->second.size();
		}

		return num;
	}

	void Reset()
	{
		m_rules.clear();
		m_args.clear();
		m_sourceKeys.clear();
		m_programName.clear();
	}

private: // types

	struct Rule
	{
		Rule(char valueNum, bool isOptional, const std::string &dependsOn,
			 const std::string &storageKey)
			: valueNum(valueNum),
			  isOptional(isOptional),
			  dependsOn(dependsOn),
			  storageKey(storageKey)
		{
		}

		char valueNum;
		bool isOptional;
		std::string dependsOn;
		std::string storageKey;
	};

	typedef std::pair<std::string, Rule> NamedRule;
	typedef std::vector<std::string> Arg;
	typedef std::pair<std::string, Arg> NamedArg;

private: // fields

	std::map<std::string, Rule> m_rules;
	std::map<std::string, Arg> m_args;
	std::map<std::string, std::string> m_sourceKeys;
	std::string m_programName;

private: // methods

	std::string FindArgValue(const std::string &key, size_t valueIndex) const
	{
		Arg arg = FindArg(key);

		if (valueIndex >= arg.size())
		{
			std::string msg = "cannot find value " + std::to_string(valueIndex)
					+ " of argument '" + key + '\'';
			Error(msg);
		}

		return arg[valueIndex];
	}

	void FillArgs(const std::vector<std::string> &rawArgs)
	{
		std::string key;
		size_t valueNum = 0;

		for (size_t i = 0; i < rawArgs.size(); ++i)
		{
			const std::string &rawArg = rawArgs[i];

			if (valueNum == 0) // read key
			{
				std::string sourceKey = RetrieveKey(rawArg);
				if (sourceKey.empty())
				{
					Error("unrecognized argument '" + rawArg + "'.\n"
						  "  Fix: use -X for a one-character flag or --name for a long flag; "
						  "run with --help to see valid names.");
				}
				Rule rule = FindRule(sourceKey);
				key = rule.storageKey;
				if (m_args.find(key) != m_args.end())
				{
					Error("duplicate option '" + FormatFlag(sourceKey) + "'.\n"
						  + "  " + FormatFlag(m_sourceKeys[key])
						  + " already sets the same option.\n"
						  + "  Fix: remove one spelling; list-valued options accept all values after one flag.");
				}
				m_args.insert(NamedArg(key, Arg()));
				m_sourceKeys[key] = sourceKey;
				valueNum = rule.valueNum;
			}
			else // read key value (one of these)
			{
				FillArg(key, rawArg, i, valueNum);
			}
		}

		if (valueNum != 0 && valueNum != '+' && valueNum != '*')
		{
			std::string msg = "not enough values for " + FormatFlag(m_sourceKeys[key])
					+ ": expected " + std::to_string(GetExpectedArgs(key))
					+ " value(s), got " + std::to_string(m_args[key].size())
					+ "\n  Fix: provide the missing value(s); run with --help for the expected syntax.";
			Error(msg);
		}
		// '+' requires at least 1 value
		if (valueNum == '+' && m_args.count(key) && m_args[key].empty())
		{
			Error(FormatFlag(m_sourceKeys[key]) + " requires at least one value.\n"
				  "  Fix: add a value after the flag; run with --help for the expected syntax.");
		}
	}

	void FillArg(const std::string &key, const std::string &rawArg,
				 size_t &i, size_t &valueNum)
	{
		if (valueNum == '+') // more than one args
		{
			if (IsFlagToken(rawArg))
			{
				--i;
				valueNum = 0;

				if (!m_args[key].empty())
				{
					return; // TODO: kostyl'
				}
			}

			AddArgValue(key, rawArg);
		}
		else if (valueNum == '*') // one or more args
		{
			if (IsFlagToken(rawArg))
			{
				--i;
				valueNum = 0;
			}
			else
			{
				AddArgValue(key, rawArg);
			}
		}
		else // fixed number of args
		{
			AddArgValue(key, rawArg);
			--valueNum;
		}
	}

	const Rule &FindRule(const std::string &key)
	{
		auto it = m_rules.find(key);

		if (it == m_rules.end())
		{
			std::string dash = (key.size() == 1) ? "-" : "--";
			std::string msg = "unknown flag '" + dash + key + "'";
			std::vector<std::pair<size_t, std::string> > candidates;
			std::set<std::string> seen;
			for (const auto &registered : m_rules)
			{
				if (registered.first.size() < 2 || !seen.insert(registered.first).second)
					continue;
				size_t distance = EditDistance(key, registered.first);
				size_t limit = std::max((size_t)2, key.size() / 3);
				if (distance <= limit)
					candidates.push_back(std::make_pair(distance, registered.first));
			}
			std::sort(candidates.begin(), candidates.end());
			if (!candidates.empty())
			{
				msg += "\n  Did you mean:";
				for (size_t i = 0; i < std::min((size_t)3, candidates.size()); ++i)
					msg += "\n    " + FormatFlag(candidates[i].second);
			}
			msg += "\n  Fix: replace the flag with a valid name shown above, or run with --help "
				   "for the full list.";
			Error(msg);
		}

		return it->second;
	}

	const Arg &FindArg(const std::string &key) const
	{
		std::string normalized = NormalizeKey(key);
		auto it = m_args.find(normalized);

		if (it == m_args.end())
		{
			std::string msg = "argument with key '" + normalized + "' not found";
			Error(msg);
		}

		return it->second;
	}

	void AddArgValue(const std::string &key, const std::string &value)
	{
		if (NonKey(value))
		{
			m_args[key].push_back(value);

			// TODO: kostyl'
			(void)FindRule(key);
		}
		else
		{
			Error("unexpected '" + value + "' while reading values for "
				  + FormatFlag(m_sourceKeys[key]) + "\n"
				  "  Got so far: --" + key + " " + JoinArg(key) + "\n"
				  "  Fix: check the argument order and provide exactly the number of values "
				  "shown by --help.");
		}
	}

	void CheckDependencies()
	{
		for (auto &arg : m_args)
		{
			const Rule &rule = FindRule(arg.first);
			if (!rule.dependsOn.empty()
				&& m_args.find(NormalizeKey(rule.dependsOn)) == m_args.end())
			{
				Error(FormatFlag(m_sourceKeys[arg.first]) + " requires "
					  + FormatFlag(rule.dependsOn) + ".\n"
					  + "  Fix: add " + FormatFlag(rule.dependsOn)
					  + " anywhere in the command line.");
			}
		}
	}

	std::string JoinArg(const std::string &key) const
	{
		std::string result;
		auto it = m_args.find(key);
		if (it != m_args.end())
			for (size_t i = 0; i < it->second.size(); ++i)
			{
				if (i > 0) result += " ";
				result += it->second[i];
			}
		return result;
	}

	size_t GetExpectedArgs(const std::string &key) const
	{
		auto it = m_rules.find(NormalizeKey(key));
		if (it != m_rules.end())
		{
			char v = it->second.valueNum;
			if (v == '+' || v == '*') return 0;
			return (size_t)v;
		}
		return 0;
	}

	void Error(const std::string &msg) const
	{
		throw std::runtime_error(msg);
	}

	void AddArgRule(const std::string &key, const Rule &rule)
	{
		if (m_rules.find(key) == m_rules.end())
		{
			m_rules.insert(NamedRule(key, rule));
		}
		else
		{
			std::string msg = "argument with key '" + key + "' already exists";
			Error(msg);
		}
	}

	std::string RetrieveKey(const std::string &rawArg)
	{
		using namespace std;
		string key;
		size_t c = 0;
		bool isOk = false;

		if ((rawArg.size() > 1) && (IsKey(rawArg[c])))
		{
			++c;

			if (IsKey(rawArg[c]))
			{
				isOk = (rawArg.size() > 3) && isalpha(rawArg[++c], locale());
			}
			else
			{
				isOk = (rawArg.size() == 2) && isalpha(rawArg[c], locale());
			}
		}

		if (isOk)
		{
			key = rawArg.substr(c);
		}

		return key;
	}

	void CheckRequiredArgs()
	{
		for (const auto &nrule : m_rules)
		{
			Rule rule = nrule.second;
			if (nrule.first != rule.storageKey)
				continue;

			if (!rule.isOptional
					&& m_args.find(nrule.first) == m_args.end())
			{
				Error("required argument " + FormatFlag(nrule.first) + " is missing.\n"
					  "  Fix: add the argument using the syntax shown by --help.");
			}
		}
	}

	std::string NormalizeKey(const std::string &key) const
	{
		auto it = m_rules.find(key);
		return (it == m_rules.end()) ? key : it->second.storageKey;
	}

	std::string FormatFlag(const std::string &key) const
	{
		return (key.size() == 1 ? "-" : "--") + key;
	}

	size_t EditDistance(const std::string &a, const std::string &b) const
	{
		std::vector<size_t> previous(b.size() + 1), current(b.size() + 1);
		for (size_t j = 0; j <= b.size(); ++j)
			previous[j] = j;
		for (size_t i = 1; i <= a.size(); ++i)
		{
			current[0] = i;
			for (size_t j = 1; j <= b.size(); ++j)
			{
				size_t replace = previous[j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1);
				current[j] = std::min(replace,
					std::min(previous[j] + 1, current[j - 1] + 1));
			}
			previous.swap(current);
		}
		return previous[b.size()];
	}

	bool NonKey(const std::string &arg)
	{
		return !IsFlagToken(arg);
	}

	bool IsKey(const char symb)
	{
		return symb == '-';
	}

	bool IsFlagToken(const std::string &arg)
	{
		if (arg.empty() || arg[0] != '-')
			return false;
		return !IsNumberToken(arg);
	}

	bool IsNumberToken(const std::string &arg)
	{
		if (arg.empty())
			return false;
		char *end = nullptr;
		std::strtod(arg.c_str(), &end);
		return end != arg.c_str() && *end == '\0';
	}

	template <class T>
	T ConvertTo(const std::string &rawValue) const
	{
		bool ok = true;
		char *end = nullptr;
		T val = T();
		errno = 0;

		if (typeid(T) == typeid(int))
		{
			long parsed = strtol(rawValue.c_str(), &end, 10);
			if (errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX)
				ok = false;
			val = (T)parsed;
		}
		else if (typeid(T) == typeid(double))
		{
			double parsed = strtod(rawValue.c_str(), &end);
			if (errno == ERANGE || !std::isfinite(parsed))
				ok = false;
			val = (T)parsed;
		}
		else
		{
			ok = false;
		}

		if (end == rawValue.c_str() || (end && strlen(end) != 0) || !ok)
		{
			std::string typeName = (typeid(T) == typeid(int)) ? "integer" : "floating-point number";
			std::string msg = "cannot convert value '" + rawValue
					+ "' to " + typeName
					+ ".\n  Fix: replace it with a finite " + typeName + ".";
			Error(msg);
		}

		return val;
	}
};
