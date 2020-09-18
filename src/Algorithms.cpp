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

#include <iostream>
#include <fstream>
#include "Algorithms.h"
#include "Scoring.h"
#include "SubstitutionMatrix.h"
#include "Seed.h"

/**
 * Calculate fswm distance between reads and genomes considering all spaced words.
 */
bool Algorithms::fswm_complete(BucketManager &genomeBucketManager, 
				BucketManager &readBucketManager, Scoring &fswm_distances) {
	SubstitutionMatrix substMat;
	std::ofstream histogramFile;
	if (fswm_params::g_writeHistogram) { histogramFile.open(fswm_params::g_outfoldername + "histogram.txt", std::ios_base::app); }

	// Loop through minimizers and compare each bucket on its own
	std::cout << genomeBucketManager.minimizersToBuckets.size() << std::endl;
	for (std::pair<minimizer_t, std::vector<std::pair<seq_id_t, pos_t>>> minimizerToBucket : genomeBucketManager.minimizersToBuckets) {
		
		// If no spaced words for this minimizer in reads, continue
		minimizer_t minimizer = minimizerToBucket.first;
		if (readBucketManager.minimizersToBuckets.find(minimizer) == readBucketManager.minimizersToBuckets.end()) {
			continue;
		}

		// Loop through buckets and compare spaced words
		std::vector<std::pair<seq_id_t, pos_t>> wordsGenomes = minimizerToBucket.second;
		std::vector<std::pair<seq_id_t, pos_t>> wordsReads = readBucketManager.minimizersToBuckets[minimizer];

		if (fswm_params::g_verbose) {
			std::cout << "\tBucket: " << minimizer << std::endl;
			std::cout << "\tEntries Genomes: " << wordsGenomes.size() << std::endl;
			std::cout << "\tEntries Reads: " << wordsReads.size() << std::endl;
		}

		int count = 0;
		// Loop through all word groups
		for (auto const &wordGenome : wordsGenomes) {
			for (auto const &wordRead : wordsReads) {

				// For each match calculate spaced word score
				for (auto &seed : fswm_internal::seeds) {
					std::vector<int> dontCarePos = seed.get_dontCarePos();
					std::cout << "Test" << std::endl;
					int score = 0;
					int mismatches = 0;
					for (auto const &pos : dontCarePos) {
						score += substMat.chiaromonte
									[fswm_internal::referenceSequences[wordGenome.first][wordGenome.second + pos]]
									[fswm_internal::querySequences[wordRead.first][wordRead.second + pos]];
						mismatches += substMat.mismatch
									[fswm_internal::referenceSequences[wordGenome.first][wordGenome.second + pos]]
									[fswm_internal::querySequences[wordRead.first][wordRead.second + pos]];
            		}
            		std::cout << "Test" << std::endl;

					if (fswm_params::g_writeHistogram) {
						int readSeqID = wordRead.first;
						int genomeSeqID = wordGenome.first;
						histogramFile << readSeqID << "\t" << genomeSeqID << "\t" << score << std::endl;
					}

					if (score > fswm_params::g_filteringThreshold) {
						count++;
						int readSeqID = wordRead.first;
						int genomeSeqID = wordGenome.first;
						// std::cout << readSeqID << "\t" << genomeSeqID << "\t" << wordsReads[wordRead_it->first + readCounter].seqPos << std::endl;
						if (fswm_distances.scoringMap.find(readSeqID) == fswm_distances.scoringMap.end()) {
							fswm_distances.scoringMap[readSeqID] = std::unordered_map<seq_id_t, scoring_t>();
							fswm_distances.mismatchCount[readSeqID] = std::unordered_map<seq_id_t, count_t>();
							fswm_distances.spacedWordMatchCount[readSeqID] = std::unordered_map<seq_id_t, count_t>();
						}
						fswm_distances.scoringMap[readSeqID][genomeSeqID] += score;
						fswm_distances.mismatchCount[readSeqID][genomeSeqID] += mismatches;
						fswm_distances.spacedWordMatchCount[readSeqID][genomeSeqID] += 1;
					}
				}
			}
		}
		std::cout << "\t\t# matches: " << count << std::endl;
	}

	if (fswm_params::g_writeHistogram) { histogramFile.close(); }

	return true;
}
