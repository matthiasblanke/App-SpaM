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
#include "GenomeManager.h"
#include "SeqIO.h"
#include "GlobalParameters.h"

GenomeManager::GenomeManager(std::string genomesfname, std::vector<Seed> &seeds) {
	if (fswm_params::g_verbose) { std::cout << "-> Reading genomes from file: " << genomesfname << std::endl; }

	bucketManagerGenomes = BucketManager(false);

	SeqIO::read_sequences(genomesfname, false, bucketManagerGenomes);
}

/**
 * Return BucketManager for the next partition of genomes from input sequences.
 */
BucketManager& GenomeManager::get_BucketManager() {
	return (bucketManagerGenomes);
}

uint32_t GenomeManager::get_genomeCount() const {
	return genomeCount;
}
