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

#include <numeric>
#include "BucketManager.h"

BucketManager::BucketManager() {
	std::vector<minimizer_t> minimizersList(16);
	std::iota (std::begin(minimizersList), std::end(minimizersList),0);

	for (auto const minimizer : minimizersList) {
		minimizers.push_back(minimizer);
		minimizersToBuckets.insert(std::unordered_map<minimizer_t, Bucket>::value_type (minimizer, Bucket(minimizer)));
	}
}

bool BucketManager::sort_words_in_buckets() {
 	for (auto &minimizerToBucket : minimizersToBuckets) {
		minimizerToBucket.second.sort_words();
	}
	return true;
}

bool BucketManager::create_wordGroups() {
 	for (auto &minimizerToBucket : minimizersToBuckets) {
		minimizerToBucket.second.create_wordGroups();
	}
	return true;
}

bool BucketManager::print_bucket_information() const {
	for (auto & bucket: minimizersToBuckets) {
		std::cout << bucket.first << ": " << bucket.second.get_bucketSize() << std::endl;
		std::cout << "Is sorted: " << bucket.second.words_sorted() << std::endl << std::endl;
	}
	return true;
}
