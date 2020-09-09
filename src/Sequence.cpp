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

#include "Sequence.h"
#include "Crc32.h"

/**
 * Create a single sequence and its reverse complement.
 * @param header Name of sequence (after '>' in fasta format).
 * @param seqLine Reference to string of nucleotide sequence itself.
 */
Sequence::Sequence(std::string &header, std::string &seqLine, seq_id_t seqID) {
	this->header = header;
	this->seqID = seqID;
	seq.reserve(seqLine.size());
	seqRev.reserve(seqLine.size());

	// Create normal sequence and its reverse
	for (std::string::iterator it = seqLine.begin(); it != seqLine.end(); it++) {
		if (!std::isspace(*it)) {
			char c = std::toupper(*it);
			switch(c) {
				case 'A': seq.push_back(0x00); break;
				case 'C': seq.push_back(0x01); break;
				case 'G': seq.push_back(0x02); break;
				case 'T': seq.push_back(0x03); break;
				case 'U': seq.push_back(0x03); break;
				default:  break;
			}
		}
	}

	seqRev.reserve(seq.size());
	for (std::string::reverse_iterator rit = seqLine.rbegin(); rit != seqLine.rend(); rit++) {
		if (!std::isspace(*rit)) {
			char c = std::toupper(*rit);
			switch(c) {
				case 'A': seqRev.push_back(0x03); break;
				case 'C': seqRev.push_back(0x02); break;
				case 'G': seqRev.push_back(0x01); break;
				case 'T': seqRev.push_back(0x00); break;
				case 'U': seqRev.push_back(0x00); break;
				default: break;
			}
		}
	}
}

/**
 * Go through sequence and fill buckets with spaced words in sequence.
 */
void Sequence::fill_buckets(std::vector<Seed> &seeds, BucketManager &bucketManager) {

	const size_t NumBytes = 8;

	for (auto &seed : seeds) {
		std::vector<int> matchPos = seed.get_matchPos();
		std::vector<int> dontCarePos = seed.get_dontCarePos();

	// Go through all spaced words in sequence and save to bucket
	uint32_t go_until = std::max(int(seq.size() - fswm_params::g_weight - fswm_params::g_spaces + 1), 0);
	for (uint32_t i = 0; i < go_until; i++) {
		// Create spaced word at position i and give word to BucketManager for further processing
		word_t matches = 0;
		for (auto const &pos : matchPos) {
			matches = matches << 2;
			matches += seq[i+pos];
		}
	
		// Create spaced word at position i and give word to BucketManager for further processing
		word_t dontCares = 0;
		for (auto const &pos : dontCarePos) {
			dontCares = dontCares << 2;
			dontCares += seq[i+pos];
		}

		if (fswm_params::g_sampling) {
			auto res = crc32_fast(&matches, NumBytes);
			if (res < fswm_params::g_minHashLowerLimit) {
				Word newWord = Word(seqID, i, matches, dontCares);
				bucketManager.insert_word(newWord);
			}
		}
		else {
			Word newWord = Word(seqID, i, matches, dontCares);
			bucketManager.insert_word(newWord);
		}
	}

	for (uint32_t i = 0; i < go_until; i++) {
		// Create spaced word at position i and give word to BucketManager for further processing
		word_t matches = 0;

		for (auto const &pos : matchPos) {
			matches = matches << 2;
			matches += seqRev[i+pos];
		}
	
		// Create spaced word at position i and give word to BucketManager for further processing
		word_t dontCares = 0;
		for (auto const &pos : dontCarePos) {
			dontCares = dontCares << 2;
			dontCares += seqRev[i+pos];
		}

		if (fswm_params::g_sampling) {
			auto res = crc32_fast(&matches, NumBytes);
			if (res < fswm_params::g_minHashLowerLimit) {
				Word newWord = Word(seqID, i, matches, dontCares);
				bucketManager.insert_word(newWord);
			}
		}
		else {
			Word newWord = Word(seqID, i, matches, dontCares);
			bucketManager.insert_word(newWord);
		}
		
	}
	}
}
