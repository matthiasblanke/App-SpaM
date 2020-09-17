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

ReadManager::ReadManager(std::string readsfname) {
	if (fswm_params::g_verbose) { std::cout << "-> Reading reads from file: " << readsfname << std::endl; }

	SeqIO::read_sequences(readsfname, reads);
	if (fswm_params::g_verbose) { std::cout << "\t" << reads.size() << " reads found and read."<< std::endl; }

	partitions = ceil((double) reads.size() / fswm_params::g_readBlockSize);
	if (fswm_params::g_verbose) { std::cout << "\tDividing into " << partitions << " partitions" << std::endl; }

	this->readCount = reads.size();
	this->currentSeq = 0;
	this->currentPartition = 0;
}

/**
 * Return BucketManager for the next partition of reads from input read sequences.
 */
void ReadManager::get_next_partition_BucketManager(std::vector<Seed> &seeds, BucketManager &bucketManagerReads) {

	if (fswm_params::g_verbose) { std::cout << "\t-> Creating spaced words for read partition " << currentPartition << std::endl; }

	fswm_internal::readIDsToNames.clear();
	fswm_internal::namesToReadIDs.clear();

	for (int i = 0; i < fswm_params::g_readBlockSize and currentSeq < reads.size(); i++) {
		reads[currentSeq].fill_buckets(seeds, bucketManagerReads);
		fswm_internal::readIDsToNames[reads[currentSeq].get_seqID()] = reads[currentSeq].get_header();
		fswm_internal::namesToReadIDs[reads[currentSeq].get_header()] = reads[currentSeq].get_seqID();		
		currentSeq++;
	}

	currentPartition++;
}

uint32_t ReadManager::get_partitions() const {
	return partitions;
}

std::vector<Sequence>& ReadManager::get_reads() {
	return reads;
}

uint32_t ReadManager::get_readCount() const {
	return readCount;
}