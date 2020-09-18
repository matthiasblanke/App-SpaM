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

/**
 * Functionality:
 * Reads the sequences in partitions instead of all at once. 
 * Delivers the BucketManager of partitions one after another 
 * via get_next_partition_BucketManager.
 *
 * Example:
 * 	ReadManager	readManager(readBlockSize, readsfname);
 * 	for (int currentPartition = 0; currentPartition < readManager.get_partitions(); currentPartition++) {
 * 		BucketManager bucketManagerReads = readManager.get_next_partition_BucketManager(seed);
 *	}
 *
 */
#ifndef FSWM_READMANAGER_H_
#define FSWM_READMANAGER_H_

#include <string>
#include <vector>
#include "GlobalParameters.h"
#include "Seed.h"
#include "BucketManager.h"

class ReadManager {
	private:
    	uint32_t readCount;
    	BucketManager bucketManagerReads;

	public:
		ReadManager(std::string readsfname, std::vector<Seed> &seeds);
		BucketManager get_BucketManager();

		uint32_t get_readCount() const;
};

#endif