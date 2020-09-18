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

#include <math.h>
#include "ReadManager.h"
#include "SeqIO.h"

ReadManager::ReadManager(std::string readsfname, std::vector<Seed> &seeds) {
	if (fswm_params::g_verbose) { std::cout << "-> Reading reads from file: " << readsfname << std::endl; }

	bucketManagerReads = BucketManager(true);

	SeqIO::read_sequences(readsfname, false, bucketManagerReads);
}

/**
 * Return BucketManager for the next partition of genomes from input sequences.
 */
BucketManager ReadManager::get_BucketManager() {
	return (bucketManagerReads);
}

uint32_t ReadManager::get_readCount() const {
	return readCount;
}