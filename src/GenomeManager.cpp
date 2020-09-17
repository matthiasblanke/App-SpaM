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
#include "Word.h"
#include "SeqIO.h"
#include "GlobalParameters.h"

GenomeManager::GenomeManager(std::string genomesfname, std::vector<Seed> &seeds) {
	if (fswm_params::g_verbose) { std::cout << "-> Reading genomes from file: " << genomesfname << std::endl; }

	bucketManagerGenomes = BucketManager();

	std::vector<Sequence> genomes;
	SeqIO::read_sequences(genomesfname, genomes);
	if (fswm_params::g_verbose) { std::cout << "\t" << genomes.size() << " genomes found and read."<< std::endl; }
	fswm_internal::g_numberGenomes = genomes.size();

	this->genomeCount = genomes.size();

	if (fswm_params::g_verbose) { std::cout << "\t" << "Creating spaced words for genomes." << std::endl; }
	for (auto &genome : genomes) {
		genome.fill_buckets(seeds, bucketManagerGenomes);
		
		std::string header = genome.get_header();
		seq_id_t seqID = genome.get_seqID();

		if (fswm_internal::seqIDsToNames.find(seqID) != fswm_internal::seqIDsToNames.end() or fswm_internal::namesToSeqIDs.find(header) != fswm_internal::namesToSeqIDs.end()) {
			std::cerr << "Multiple sequences in the genomes seem to have the same name. Please fix." << std::endl;
			exit(EXIT_FAILURE);
		}
		fswm_internal::seqIDsToNames[seqID] = header;
		fswm_internal::namesToSeqIDs[header] = seqID;
		fswm_internal::genomeIDsToNames[seqID] = header;
		fswm_internal::namesToGenomeIDs[header] = seqID;
	}

	//genomes.clear();
	//genomes.shrink_to_fit();
}

/**
 * Return BucketManager for the next partition of genomes from input sequences.
 */
BucketManager GenomeManager::get_BucketManager() {
	return (bucketManagerGenomes);
}

std::vector<Sequence>& GenomeManager::get_genomes() {
	return genomes;
}

uint32_t GenomeManager::get_genomeCount() const {
	return genomeCount;
}
