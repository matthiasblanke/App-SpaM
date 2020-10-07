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
#ifndef FSWM_SCORING_H_
#define FSWM_SCORING_H_

#include <unordered_map>
#include <string>
#include <vector>
#include "Word.h"

// Unordered map from sequence IDs to integer counts used e.g. for mismatches and number of spaced words
typedef std::unordered_map<seq_id_t,count_t> seqIDtoCount_t;

// Unordered map from sequence IDs to seqIDtoCount_t for "matrix"-functionality
typedef std::unordered_map<seq_id_t,seqIDtoCount_t> countMap_t;

// Unordered map from sequence IDs to scoring_t used e.g. for fswm-scores
typedef std::unordered_map<seq_id_t,scoring_t> seqIDtoScoring_t;

// Unordered map from sequence IDs to seqIDtoScoring_t for "matrix"-functionality
typedef std::unordered_map<seq_id_t,seqIDtoScoring_t> scoringMap_t;

class Scoring {
	public:
		// For each assigned read (first seqID) it records the seqID of the assigned genome or internal leave (second seqID)
		std::vector<std::pair<seq_id_t, int>> readAssignment;

		// Maps for seqIDs of reads to map of seqIDs of genomes to scoring/mismatchCount/spacedWordMatchCount
		countMap_t kmerCountsMap;
		scoringMap_t scoringMap;
		countMap_t mismatchCount;
		countMap_t spacedWordMatchCount;

		Scoring();

		/**
		 * Calculate jk-corrected distances between fswm based on mismatch counts.
		 */
		void calculate_fswm_distances();

		/**
		 * Assign reads to reference genome based on shortest jk-corrected distance.
		 */
		void read_assignment();

		/**
		 * Assign reads to reference tree of genome.
		 */
		void phylogenetic_placement();

		/*
		Assign reads to reference tree of genome (only phylo-kmer based)
		*/
		void phylo_kmer_placement();
		void phylo_kmer_placement_best();
		void phylo_kmer_placement_all_best();
		void phylo_kmer_placement_path();

		/**
		 * Assign reads to reference genome based on most k-mer matches.
		 */
		void assign_reads_kmer_counting();

		/**
		 * Write read assignments to file.
		 * @param filename Filename to be written to.
		 */
		void write_assignment_to_file();

		// Debug functions
		/**
		 * Write jk-corrected distances between all reads and genomes to file.
		 * The file will be located in the specified results-folder (defaul: results) and named scoring.txt
		 */
		void write_scoring_to_file();

		/**
		 * Write jk-corrected distances between all reads and genomes to file in tab-delimited table.
		 */
		void write_scoring_to_file_as_table();
};

#endif
