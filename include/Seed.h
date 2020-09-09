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
#ifndef FSWM_SEED_H_
#define FSWM_SEED_H_

#include <vector>
#include "Pattern.h"

class Seed{
	private:
		int32_t length;
		int32_t weight;
		int32_t dontCare;
		int32_t spaces;

		std::vector<int> matchPos;
		std::vector<int> dontCarePos;

	public:
		Seed(int32_t weight, int32_t dontCare);
		bool generate_pattern(std::string &patternStr);

		int32_t get_length();
		int32_t get_weight();
		int32_t get_dontCare();
		int32_t get_spaces();

		std::vector<int> get_matchPos() const;
		std::vector<int> get_dontCarePos() const;
};


/**
 * Return pattern length.
 */
inline int32_t Seed::get_length() {
	return length;
}

/**
 * Return pattern weight (number of match positions).
 */
inline int32_t Seed::get_weight() {
	return weight;
}

/**
 * Return the number of dont care positions
 */
inline int32_t Seed::get_spaces() {
	return spaces;
}

/**
 * Return the number of positions where matches occur in the pattern.
 */
inline std::vector<int> Seed::get_matchPos() const {
	return matchPos;
}

/**
 * Return the number of positions where dont care's occur in the pattern.
 */
inline std::vector<int> Seed::get_dontCarePos() const {
	return dontCarePos;
}

#endif
