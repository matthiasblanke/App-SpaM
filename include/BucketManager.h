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
#ifndef FSWM_BucketManager_H_
#define FSWM_BucketManager_H_

#include <unordered_map>
#include "Bucket.h"
#include "Word.h"
#include "Scoring.h"

class BucketManager {

	private:
		int32_t bucketCount;
		std::unordered_map<minimizer_t,Bucket> minimizersToBuckets;
		std::vector<minimizer_t> minimizers;

	public:
		// Constructors
		BucketManager();

		// Functions
		bool insert_word(Word &word);
		bool sort_words_in_buckets();
		bool create_wordGroups();

		// Debug functions
		bool print_bucket_information() const;

		// Get & Set
		int get_bucketCount() const;
		std::vector<minimizer_t> get_minimizers();
		Bucket get_bucket(minimizer_t minimizer);
};

inline bool BucketManager::insert_word(Word &word) {
	minimizersToBuckets.find(word.minimizer)->second.add_word(word);
	return true;
}

inline int BucketManager::get_bucketCount() const {
	return bucketCount;
}

inline std::vector<minimizer_t> BucketManager::get_minimizers() {
	return minimizers;
}

inline Bucket BucketManager::get_bucket(minimizer_t minimizer) {
	return minimizersToBuckets.find(minimizer)->second;
}

#endif
