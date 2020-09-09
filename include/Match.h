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
#ifndef FSWM_MATCHES_H_
#define FSWM_MATCHES_H_

#include "Word.h"

struct Match {
	int score;
	int mismatches;
	seq_id_t readID;
	seq_id_t genomeID;
	pos_t readPos;
	pos_t genomePos;

	Match(int score, int mismatches, seq_id_t readID, seq_id_t genomeID, pos_t readPos, pos_t genomePos);

	bool operator<(const Match& otherMatch) const {
		return score < otherMatch.score;
	}
	bool operator>(const Match& otherMatch) const {
		return score > otherMatch.score;
	}
	bool operator==(const Match& otherMatch) const {
		return score == otherMatch.score;
	}

	friend std::ostream& operator<<(std::ostream &strm, const Match &match) {
		return strm << "Match(score=" << match.score << ", mismatches=" << match.mismatches << ", readID=" << match.readID
			  << ", readPos=" << match.readPos << ", genomeID=" << match.genomeID << ", genomePos=" << match.genomePos << ")";
	}
};

#endif
