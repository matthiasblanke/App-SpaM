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
#ifndef FSWM_SEQUENCE_H_
#define FSWM_SEQUENCE_H_

#include <vector>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <bitset>
#include "Word.h"
#include "Seed.h"
#include "Bucket.h"
#include "BucketManager.h"

class Sequence {

	private:
		std::vector<char> seq;
		std::vector<char> seqRev;
		std::string header;
		seq_id_t seqID;

	public:
		Sequence(std::string &header, std::string &seqLine, seq_id_t seqID);

		void fill_buckets(std::vector<Seed> &seeds, BucketManager &bucket_manager);

		std::string get_header() const;
		seq_id_t get_seqID() const;
};

inline std::string Sequence::get_header() const {
	return header;
}

inline seq_id_t Sequence::get_seqID() const {
	return seqID;
}

#endif
