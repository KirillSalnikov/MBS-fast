#include "Tracks.h"

#include <climits>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Beam.h"

int Tracks::FindGroupByTrackId(const IdType &trackId) const
{
	for (size_t i = 0; i < size(); ++i)
	{
		for (int j = 0; j < (*this)[i].size; ++j)
		{
			if ((*this)[i].arr[j] == trackId)
			{
				return (*this)[i].groupID;
			}
		}
	}

	if (size() == 0)
	{
		return 0;
	}

	return -1;
}

void Tracks::ImportTracks(int nFacets, const std::string &filename)
{
	std::ifstream trackFile(filename.c_str());
	if (!trackFile)
		throw std::runtime_error(
			"cannot read trajectory file '" + filename + "'.\n"
			"Fix: provide a readable text file with one facet-index sequence per line.");

	clear();
	TrackGroup ungrouped;
	std::string line;
	int lineNumber = 0;
	while (std::getline(trackFile, line))
	{
		++lineNumber;
		const std::string::size_type comment = line.find('#');
		if (comment != std::string::npos)
			line.erase(comment);

		std::istringstream input(line);
		std::vector<std::string> tokens;
		std::string token;
		while (input >> token)
			tokens.push_back(token);
		if (tokens.empty())
			continue;

		std::vector<int> track;
		bool haveGroup = false;
		size_t groupIndex = 0;
		for (size_t i = 0; i < tokens.size(); ++i)
		{
			if (tokens[i] == ":")
			{
				if (haveGroup || i + 2 != tokens.size())
					throw std::runtime_error(
						"trajectory file '" + filename + "', line "
						+ std::to_string(lineNumber)
						+ ": ':' must be followed by exactly one group index.\n"
						"Fix: use 'FACET ... : GROUP' or remove the malformed ':'.");
				haveGroup = true;
				char *end = nullptr;
				const long parsed = std::strtol(tokens[++i].c_str(), &end, 10);
				if (!end || *end != '\0' || parsed < 0 || parsed >= MAX_GROUP_NUM)
					throw std::runtime_error(
						"trajectory file '" + filename + "', line "
						+ std::to_string(lineNumber)
						+ ": group index must be an integer from 0 to "
						+ std::to_string(MAX_GROUP_NUM - 1) + ".\n"
						"Fix: replace the group index with a value in the supported range.");
				groupIndex = static_cast<size_t>(parsed);
				continue;
			}

			char *end = nullptr;
			const long parsed = std::strtol(tokens[i].c_str(), &end, 10);
			if (!end || *end != '\0' || parsed < 0 || parsed >= nFacets)
				throw std::runtime_error(
					"trajectory file '" + filename + "', line "
					+ std::to_string(lineNumber) + ": facet index '"
					+ tokens[i] + "' is outside [0, "
					+ std::to_string(nFacets - 1) + "].\n"
					"Fix: use zero-based facet indices from the loaded particle.");
			track.push_back(static_cast<int>(parsed));
		}

		if (track.empty())
			throw std::runtime_error(
				"trajectory file '" + filename + "', line "
				+ std::to_string(lineNumber) + ": trajectory has no facet indices.\n"
				"Fix: place at least one zero-based facet index before the optional ': GROUP'.");

		IdType trackID = 0;
		for (int facet : track)
		{
			trackID += facet + 1;
			trackID *= nFacets + 1;
		}

		TrackGroup *group = &ungrouped;
		if (haveGroup)
		{
			if (groupIndex >= size())
				resize(groupIndex + 1);
			group = &(*this)[groupIndex];
			group->groupID = static_cast<int>(groupIndex);
		}
		if (group->size >= MAX_GROUP_NUM)
			throw std::runtime_error(
				"trajectory file '" + filename + "', line "
				+ std::to_string(lineNumber) + ": group contains more than "
				+ std::to_string(MAX_GROUP_NUM) + " trajectories.\n"
				"Fix: split the trajectories across additional groups.");
		group->tracks.push_back(track);
		group->arr[group->size++] = trackID;
	}

	if (!trackFile.eof())
		throw std::runtime_error(
			"failed while reading trajectory file '" + filename + "'.\n"
			"Fix: check file permissions and remove overlong or binary records.");
	if (empty() && ungrouped.size == 0)
		throw std::runtime_error(
			"trajectory file '" + filename + "' contains no trajectories.\n"
			"Fix: add one or more lines of zero-based facet indices.");

	for (int i = 0; i < ungrouped.size; ++i)
	{
		if (size() >= MAX_GROUP_NUM)
			throw std::runtime_error(
				"trajectory file '" + filename + "' creates more than "
				+ std::to_string(MAX_GROUP_NUM) + " groups.\n"
				"Fix: combine trajectories with the ': GROUP' syntax.");
		TrackGroup group;
		group.arr[group.size++] = ungrouped.arr[i];
		group.tracks.push_back(ungrouped.tracks[i]);
		group.groupID = static_cast<int>(size());
		push_back(group);
	}
}

void Tracks::RecoverTrack(const Beam &beam, int facetNum,
						  std::vector<int> &track)
{
	const int coef = facetNum + 1;
	const size_t oldSize = track.size();
	const size_t trackLength = static_cast<size_t>(beam.nActs + 1);
	track.resize(oldSize + trackLength);

#ifndef _DEBUG
	// Most practical trajectories fit in at most two BigInteger blocks. Decode
	// those with native integer arithmetic instead of repeated arbitrary-
	// precision quotient/remainder operations in the diffraction hot path.
#if defined(__SIZEOF_INT128__)
	__extension__ typedef unsigned __int128 NativeTrackId;
	const BigInteger::Index nativeBlockLimit = 2;
#else
	using NativeTrackId = unsigned long;
	const BigInteger::Index nativeBlockLimit = 1;
#endif
	if (beam.id.getSign() != BigInteger::negative
		&& beam.id.getLength() <= nativeBlockLimit)
	{
		NativeTrackId tmpId = beam.id.isZero()
			? NativeTrackId(0) : NativeTrackId(beam.id.getBlock(0));
#if defined(__SIZEOF_INT128__)
		if (beam.id.getLength() == 2)
		{
			tmpId |= NativeTrackId(beam.id.getBlock(1))
				<< (sizeof(unsigned long) * CHAR_BIT);
		}
#endif
		tmpId /= static_cast<NativeTrackId>(coef);
		for (int i = 0; i <= beam.nActs; ++i)
		{
			const int tmp = static_cast<int>(
				tmpId % static_cast<NativeTrackId>(coef));
			tmpId /= static_cast<NativeTrackId>(coef);
			if (tmp < 1 || tmp > facetNum)
				throw std::runtime_error(
					"beam trajectory ID cannot be decoded for the loaded particle.\n"
					"Fix: unset MBS_DISABLE_TRACK_IDS when using absorption or "
					"--trajectories; if it is already unset, report this command as a bug.");
			track[oldSize + trackLength - 1 - static_cast<size_t>(i)] = tmp - 1;
		}
		return;
	}
#endif

	auto tmpId = beam.id/coef;

	for (int i = 0; i <= beam.nActs; ++i)
	{
#ifdef _DEBUG // DEB
		int tmp = (tmpId%coef);
#else
		int tmp = (tmpId%coef).toInt();
#endif
		tmpId -= tmp;
		tmpId /= coef;
		if (tmp < 1 || tmp > facetNum)
			throw std::runtime_error(
				"beam trajectory ID cannot be decoded for the loaded particle.\n"
				"Fix: unset MBS_DISABLE_TRACK_IDS when using absorption or "
				"--trajectories; if it is already unset, report this command as a bug.");
		track[oldSize + trackLength - 1 - static_cast<size_t>(i)] = tmp - 1;
	}
}
