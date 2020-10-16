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
#include "Seed.h"
#include "GlobalParameters.h"

/**
 * Initialize Seed with length attributes and generate the pattern.
 * @param weight Weight of the pattern (number of match positions)
 * @param dontCare Number of don't care positions in pattern. weight + don'tcare = total length
 */
Seed::Seed(int32_t weight, int32_t dontCare) {
	this->weight = weight;
	this->dontCare = dontCare;
	this->length = weight + dontCare;
	this->spaces = length - weight;
}

/**
 * Generates a pattern suitable for the dataset.
 */
bool Seed::generate_pattern(std::string &patternStr) {
	for(int i = 0; i < patternStr.length(); i++) {
		if (patternStr[uint(i)] == '1') {
			matchPos.push_back(i);
		} else {
			dontCarePos.push_back(i);
		}
	}

	if (fswm_params::g_verbose) {
		std::cout << "\tGenerated pattern with weight " << weight << " and length " << length << "." << std::endl;
		std::cout << "\t\tPattern: " << patternStr << std::endl;
	}

	return true;
}
