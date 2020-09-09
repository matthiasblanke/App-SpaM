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
#ifndef FSWM_BUCKET_H_
#define FSWM_BUCKET_H_

#include <vector>
#include "Word.h"

class Bucket {
	private:
		std::vector<Word> words;
		minimizer_t minimizer;
		uint32_t bucketSize;

		// Word groups are groups of words with the same hash of matching positions.
		// They are represented by pairs of integers, the first encoding the start
		// position of the group in the sorted words vector, the second encoding the
		// number of elements in this group
		std::vector<std::pair<int, int>> wordGroups;

	public:
		// Constructors
		Bucket(minimizer_t minimizer);

		// Functions
		bool add_word(Word &newWord);
		bool sort_words();
		bool words_sorted() const;
		bool create_wordGroups();

		// Get & Set
		std::vector<Word>& get_words();
		minimizer_t get_minimizer() const;
		uint32_t get_bucketSize() const;
		std::vector<std::pair<int, int>>& get_wordGroups();

		// Compare buckets only based on their associated minimizer
		bool operator>(const Bucket& otherBucket) const {
			return minimizer > otherBucket.minimizer;
		}

		bool operator<(const Bucket& otherBucket) const {
			return minimizer < otherBucket.minimizer;
		}

		bool operator==(const Bucket& otherBucket) const {
			return minimizer == otherBucket.minimizer;
		}
};

inline bool Bucket::add_word(Word &newWord) {
	words.push_back(newWord);
	bucketSize++;
	return true;
}

inline bool Bucket::sort_words() {
	std::sort(words.begin(), words.end());
	return true;
}

inline bool Bucket::words_sorted() const {
	return std::is_sorted(words.begin(), words.end());
}

inline std::vector<Word>& Bucket::get_words() {
	return words;
}

inline uint32_t Bucket::get_bucketSize() const {
	return bucketSize;
}

inline minimizer_t Bucket::get_minimizer() const {
	return minimizer;
}

inline std::vector<std::pair<int, int>>& Bucket::get_wordGroups() {
	return wordGroups;
}

#endif
