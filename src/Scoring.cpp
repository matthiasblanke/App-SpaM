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

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <math.h>
#include "Scoring.h"
#include "Tree.h"
#include <vector>


Scoring::Scoring() {

}

/** Calculate jk-corrected distances between fswm based on mismatch counts. */
void Scoring::calculate_fswm_distances() {
	double substFreq = 0;	// Calculated plain substitution frequency
	double jk = 0;			// Jukes-Cantor corrected substitution frequency

	scoringMap_t::iterator it1 = scoringMap.begin();
	countMap_t::iterator it2 = mismatchCount.begin();
	countMap_t::iterator it3 = spacedWordMatchCount.begin();

	while (it1 != scoringMap.end()) { 	// Iterate through reads

		seqIDtoScoring_t::iterator it11 = it1->second.begin();
		seqIDtoCount_t::iterator it21 = it2->second.begin();
		seqIDtoCount_t::iterator it31 = it3->second.begin();

		while (it11 != it1->second.end()) {		// Iterate through genomes
			if (it31->second <= 0) { 			// Reads with not matches get default distance
				it11->second = fswm_params::g_defaultDistance;
			}
			else {
				substFreq = (double) it21->second / (it31->second * fswm_params::g_spaces);
				jk = -0.75 * log(1.0 - ((4.0/3.0) * substFreq));

				it11->second = jk;
			}
			it11++;
			it21++;
			it31++;
		}
		it1++;
		it2++;
		it3++;
	}
}

/** Assign reads to reference tree of genome. */
void Scoring::phylogenetic_placement() {
	Tree tree(fswm_params::g_reftreefname);				// Read and create reference tree
	int min_j;											// Currently minimum assigned genome. -1 for unassigned.
	countMap_t::iterator countMap_it = spacedWordMatchCount.begin();

	std::unordered_map<seq_id_t, bool> readAssignmentTracker;  // Track which reads were assigned and which not (only reads that have at least one entry are assigned)
	for (const auto &read : fswm_internal::readIDsToNames) {
		readAssignmentTracker[read.first] = false;             // Later on: assign reads that have not been assigned by algorithm to root of tree
	}

   	for (scoringMap_t::iterator scoringMap_it = scoringMap.begin(); scoringMap_it != scoringMap.end(); scoringMap_it++) {		// Iterate through reads
		min_j = -1;

		if (fswm_params::g_assignmentMode == "SPAMCOUNT") {
			min_j = tree.get_node_best_count(countMap_it);
		}
		else if (fswm_params::g_assignmentMode == "MINDIST") {
			min_j = tree.get_node_best_score(scoringMap_it);
		}
		else if (fswm_params::g_assignmentMode == "LCACOUNT") {
			min_j = tree.get_LCA_best_count(countMap_it);
		}
		else if (fswm_params::g_assignmentMode == "LCADIST") {
			min_j = tree.get_LCA_best_score(scoringMap_it);
		}
		else if (fswm_params::g_assignmentMode == "SPAMX") {
			min_j = tree.get_LCA_best_count_exp(countMap_it, fswm_params::g_spam_X);
		}
		readAssignment.push_back(std::pair<seq_id_t, int> (scoringMap_it->first, min_j));  // assign read to some internal leave, determined based on assignment mode
		readAssignmentTracker[scoringMap_it->first] = true;

		countMap_it++;
   	}

   	// In the rare case, that no spaced words are found: 
   	// Assign all reads that were not assigned so far to root
   	for (const auto &read : readAssignmentTracker) {
   		if (!read.second) {
   			readAssignment.push_back(std::pair<seq_id_t, int> (read.first, tree.get_rootID()));
   		}
   	}

   	if (fswm_params::g_assignmentMode != "APPLES") {
   		tree.write_jplace_placement_data(readAssignment, this->scoringMap);
   	}

   	// tree.write_newick(fswm_params::g_outfoldername + "tree.nwk");
}

/** Write jk-corrected distances between all reads and genomes to file. */
void Scoring::write_scoring_to_file() {
	std::ofstream results;
	results.open(fswm_params::g_outfoldername + "scoring_list.txt", std::ios_base::app);

	scoringMap_t::iterator it1 = scoringMap.begin();
	while (it1 != scoringMap.end()) {
		seqIDtoScoring_t::iterator it11 = it1->second.begin();
		while (it11 != it1->second.end()) {
			results << it1->first << "\t" << it11->first << "\t" << it11->second << std::endl;
			it11++;
		}
		it1++;
	}
	results.close();
}

/** Write jk-corrected distances between reads and genomes to table. */
void Scoring::write_scoring_to_file_as_table() {
	std::ofstream results;
	results.open(fswm_params::g_outfoldername + "scoring_table.txt", std::ios_base::app);

	std::unordered_map<seq_id_t, std::string> seqIDtoNamesMap;
	seqIDtoNamesMap = fswm_internal::readIDsToNames;

	for (auto read : seqIDtoNamesMap) {		// For all reads: write distances to all genomes to file
		results << read.second;

		if (scoringMap.find(read.first) != scoringMap.end()) {					// If read has distances to any genome
			for (auto genome : fswm_internal::genomeIDsToNames) {	// Write those distances to file and use
				if (scoringMap[read.first].find(genome.first) != scoringMap[read.first].end()) {
					results << "\t" << scoringMap[read.first][genome.first];
				}
				else {
					results << "\t" << fswm_params::g_defaultDistance;
				}
			}
		}
		else {
			for (auto genome : fswm_internal::genomeIDsToNames) {	// Write those distances to file and use
				results << "\t" << fswm_params::g_defaultDistance;
			}
		}
		results << std::endl;
	}

	results.close();
}
