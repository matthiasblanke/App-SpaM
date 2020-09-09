/**
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * author: Matthias Blanke
 * mail  : matthias.blanke@biologie.uni-goettingen.de
 */

#include "MatchManager.h"
#include <algorithm>

MatchManager::MatchManager() {

}

bool MatchManager::insert_match(Match &match) {
	if (readIDsToReadPosToMatches.find(match.readID) == readIDsToReadPosToMatches.end()) {
		readIDsToReadPosToMatches.insert(matchMap::value_type (match.readID, std::unordered_map<pos_t, std::vector<Match>>()));

		readIDsToReadPosToMatches[match.readID].insert(std::unordered_map<pos_t, std::vector<Match>>::value_type (match.readPos, std::vector<Match>()));
	}
	else {
		if (readIDsToReadPosToMatches[match.readID].find(match.readPos) == readIDsToReadPosToMatches[match.readID].end()) {
			readIDsToReadPosToMatches[match.readID].insert(std::unordered_map<pos_t, std::vector<Match>>::value_type (match.readPos, std::vector<Match>()));
		}
	}
	readIDsToReadPosToMatches[match.readID][match.readPos].push_back(match);

	return true;
}

bool MatchManager::clear_readIDsToReadPosToMatches() {
	readIDsToReadPosToMatches.clear();
	return true;
}

/** Sort matches for every readID and readPos by their score. */
bool MatchManager::sort_matches() {
	for (auto readIDToReadPosToMatches : readIDsToReadPosToMatches) {
		for (auto readPosToMatches : readIDToReadPosToMatches.second) {
			std::sort(readPosToMatches.second.begin(), readPosToMatches.second.end());
		}
	}
	return true;
}

std::vector<Match> MatchManager::get_highest_scoring_matches() {
	std::vector<Match> matches;

	for (auto readIDToReadPosToMatches : readIDsToReadPosToMatches) {		// Loop through all reads
		for (auto readPosToMatches : readIDToReadPosToMatches.second) {		// Loop through all positions
			int max_score = std::numeric_limits<double>::min();
			Match max_match = readPosToMatches.second[0];

			for (auto match : readPosToMatches.second) {					// Find best match for this read and position
				if (match.score > max_score) {
					max_score = match.score;
					max_match = match;
				}
			}
			matches.push_back(max_match);
		}
	}

	return matches;
}

std::vector<Match> MatchManager::get_high_scoring_matches() {
	std::vector<Match> matches;
	sort_matches();

	for (auto readIDToReadPosToMatches : readIDsToReadPosToMatches) {		// Loop through all reads
		for (auto readPosToMatches : readIDToReadPosToMatches.second) {		// Loop through all positions

			std::unordered_map<seq_id_t, Match> referencesMatchMap;

			for (auto match : readPosToMatches.second) {					// Find best match for this read and position
				if (referencesMatchMap.find(match.genomeID) == referencesMatchMap.end()) {	// If no match for this genome found so far, insert the current one
					referencesMatchMap.insert(std::unordered_map<seq_id_t, Match>::value_type (match.genomeID, match));
				}
			}

			for (auto referencesToMatches : referencesMatchMap) {
				matches.push_back(referencesToMatches.second);
			}
		}
	}

	return matches;
}
