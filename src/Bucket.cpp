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

#include <iostream>
#include <fstream>
#include "Bucket.h"

Bucket::Bucket(minimizer_t minimizer) {
	this->minimizer = minimizer;
	bucketSize = 0;
	words.reserve(10000);
}

/**
 * Create groups of words based on same hash of matching positions. 
 */
bool Bucket::create_wordGroups() {
	this->sort_words();
	int currentGroupSize = 1;
	int currentWord_idx = 0;
	word_t currentMatchesHash;
	if (words.size() > 0) {
		currentMatchesHash = words[currentWord_idx].matches;
	}

	for (int idx = 1; idx < words.size(); idx++) {
		if (words[idx].matches == currentMatchesHash) {
			currentGroupSize++;
		}
		else {
			wordGroups.push_back(std::pair<int,int> (currentWord_idx, currentGroupSize));
			currentGroupSize = 1;
			currentWord_idx = idx;
			currentMatchesHash = words[currentWord_idx].matches;
		}
	}
	return true;
}
