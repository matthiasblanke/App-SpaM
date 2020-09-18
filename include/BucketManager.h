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
#include <vector>
#include "Scoring.h"
#include "GlobalParameters.h"

class BucketManager {

	private:
		bool isQuery;

	public:
		std::unordered_map<minimizer_t,std::vector<std::pair<seq_id_t, pos_t>>> minimizersToBuckets;

		// Constructors
		BucketManager();
		BucketManager(bool isQuery);

		// Functions
		void insert_word(minimizer_t minimizer, seq_id_t seq_id, pos_t seq_pos);
};

inline void BucketManager::insert_word(minimizer_t minimizer, seq_id_t seq_id, pos_t seq_pos) {

	if (minimizersToBuckets.find(minimizer) == minimizersToBuckets.end()){
		minimizersToBuckets[minimizer] = std::vector<std::pair<seq_id_t, pos_t>>();
	}
	minimizersToBuckets[minimizer].push_back(std::pair<seq_id_t, pos_t>(seq_id, seq_pos));
}

#endif
