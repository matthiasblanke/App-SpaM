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
#ifndef FSWM_MATCHMANAGER_H
#define FSWM_MATCHMANAGER_H

#include <unordered_map>
#include <vector>
#include "Word.h"
#include "Match.h"

typedef std::unordered_map<seq_id_t, std::unordered_map<pos_t, std::vector<Match> > > matchMap;

class MatchManager {
	private:
		/** For every read and every position in such read a vector with all matches is stored */
		matchMap readIDsToReadPosToMatches;

	public:
		MatchManager();

		bool insert_match(Match &match);				// Insert match at appropriate position into matchMap
		bool clear_readIDsToReadPosToMatches();			// Clear map that holds all matches
		bool sort_matches();							// Sort all match vectors in matchMap by their score
		std::vector<Match> get_highest_scoring_matches();	// Return the highest scoring match for each read at every position
		std::vector<Match> get_high_scoring_matches();	// Return the highest scoring match for each read at every position to each reference
};

#endif
