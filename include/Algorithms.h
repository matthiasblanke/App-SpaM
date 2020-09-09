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
#ifndef FSWM_ALGORITHMS_H_
#define FSWM_ALGORITHMS_H_

#include <vector>
#include "Word.h"
#include "BucketManager.h"

class Algorithms {		
	public:
		// Count spaced k-mers between genomes and reads
		static bool count_kmers(BucketManager &genomeBucketManager, BucketManager &readBucketManager, Scoring &fswm_distances);

		// Complete checks the quadratic number of matches between corresponding buckets
		static bool fswm_complete(BucketManager &genomeBucketManager, BucketManager &readBucketManager, Scoring &fswm_distances);

		// Calculates distances between all spaced words within one bucket
		static bool match_reads_against_phyloDB(BucketManager &readBucketManager);
};

#endif
