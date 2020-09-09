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
#include "Word.h"

/**
 * Create a (spaced) word based on hashes for matches and dont cares and calculate minimizer.
 */
Word::Word(seq_id_t seqID, pos_t seqPos, word_t matches, word_t dontCares) {
	this->seqID = seqID;
	this->seqPos = seqPos;
	this->matches = matches;
	this->dontCares = dontCares;
	this->minimizer = this->matches & 0xF;
}

/**
 * Write word with all parameters and bit representation to stream.
 */
std::ostream& operator<<(std::ostream &strm, const Word &word) {
	std::string matchesWord = "";
	word_t temp = word.matches;
	for (int i = 0; i < fswm_params::g_weight; i++) {
		int nextChar = temp & 0x0000003;
		switch (nextChar) {
			case 0: matchesWord.push_back('A'); break;
			case 1: matchesWord.push_back('C'); break;
			case 2: matchesWord.push_back('G'); break;
			case 3: matchesWord.push_back('T'); break;
		}
		temp = temp >> 2;
	}
	std::reverse(matchesWord.begin(), matchesWord.end());

	std::string dontCaresWord = "";
	temp = word.dontCares;
	for (int i = 0; i < fswm_params::g_spaces; i++) {
		int next_char = temp & 0x0000003;
		switch (next_char) {
			case 0: dontCaresWord.push_back('A'); break;
			case 1: dontCaresWord.push_back('C'); break;
			case 2: dontCaresWord.push_back('G'); break;
			case 3: dontCaresWord.push_back('T'); break;
		}
		temp = temp >> 2;
	}
	std::reverse(dontCaresWord.begin(), dontCaresWord.end());

	return strm << "Word(seqID=" << word.seqID << ", pos=" << word.seqPos << ", weight=" << fswm_params::g_weight
	  << ", spaces=" << fswm_params::g_spaces << ")" << std::endl
	  << "Matches:" << std::endl << matchesWord << std::endl
	  << "DontCares:" << std::endl << dontCaresWord << std::endl;
}
