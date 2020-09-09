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
#ifndef FSWM_WORD_H_
#define FSWM_WORD_H_

#include <string>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include "GlobalParameters.h"

struct Word {
		word_t matches;
		word_t dontCares;
		seq_id_t seqID;
		pos_t seqPos;
		minimizer_t minimizer;

		// Constructor
		Word(seq_id_t seqID, pos_t seqPos, word_t matches, word_t dontCares);

		// Words are compared by evaluating their hash of the matching positions
		bool operator<(const Word& otherWord) const {
			return matches < otherWord.matches;
		}

		bool operator>(const Word& otherWord) const {
			return matches > otherWord.matches;
		}

		bool operator==(const Word& otherWord) const {
			return (matches == otherWord.matches);
		}

		// Print the hash and bit values of the word. 
		friend std::ostream& operator<<(std::ostream &strm, const Word &word);
};

#endif
