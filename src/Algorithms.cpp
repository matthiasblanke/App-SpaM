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
#include "Match.h"
#include "MatchManager.h"

/**
 * Calculate fswm distance between reads and genomes considering all spaced words.
 */
bool Algorithms::fswm_complete(BucketManager &genomeBucketManager, Scoring &fswm_distances) {
	SubstitutionMatrix substMat;
	std::ofstream histogramFile;
	if (fswm_params::g_writeHistogram) { histogramFile.open(fswm_params::g_outfoldername + "histogram.txt", std::ios_base::app); }

	// Loop through minimizers and compare each bucket on its own
	for (auto const minimizer : genomeBucketManager.get_minimizers()) {
		//Get buckets
		Bucket bucketGenomes = genomeBucketManager.get_bucket(minimizer);

		// Loop through buckets and compare spaced words
		std::vector<Word> wordsGenomes = bucketGenomes.get_words();

		// Get vector of word groups. First int is starting position, second int length of group
		std::vector<std::pair<uint,uint>> wordGroupGenomes = bucketGenomes.get_wordGroups();

		std::vector<std::pair<uint,uint>>::const_iterator wordGenome_it1 = wordGroupGenomes.cbegin();
		std::vector<std::pair<uint,uint>>::const_iterator wordGenome_it2 = wordGroupGenomes.cbegin();

		if (fswm_params::g_verbose) {
			std::cout << "\tBucket: " << bucketGenomes.get_minimizer() << std::endl;
			std::cout << "\t\tBucket size genomes: " << bucketGenomes.get_bucketSize() << std::endl;
		}

		word_t dontCaresGenome1;
		word_t dontCaresGenome2;
		int count = 0;
		// Loop through all word groups
		while (wordGenome_it1 != wordGroupGenomes.end()) {
			if (wordsGenomes[wordGenome_it1->first] < wordsGenomes[wordGenome_it2->first]) {
				wordGenome_it1++;
			}
			else if (wordsGenomes[wordGenome_it1->first] > wordsGenomes[wordGenome_it2->first]) {
				wordGenome_it2++;
			}
			else {
				// Loop through words with same spaced k-mer
				for (int genomeCounter1 = 0; genomeCounter1 < wordGenome_it1->second; genomeCounter1++) {
					for (int genomeCounter2 = 0; genomeCounter2 < wordGenome_it2->second; genomeCounter2++) {

						// For each match calculate spaced word score
						dontCaresGenome1 = wordsGenomes[wordGenome_it1->first + genomeCounter2].dontCares;
						dontCaresGenome2 = wordsGenomes[wordGenome_it2->first + genomeCounter1].dontCares;

						int score = 0;
						int mismatches = 0;

						for (int i = 0; i < fswm_params::g_spaces; i++) {
							score += substMat.chiaromonte[(dontCaresGenome1 & 0x03)][(dontCaresGenome2 & 0x03)];
							mismatches += substMat.mismatch[(dontCaresGenome1 & 0x03)][(dontCaresGenome2 & 0x03)];
							dontCaresGenome2 = dontCaresGenome2 >> 2;
							dontCaresGenome1 = dontCaresGenome1 >> 2;
						}

						if (fswm_params::g_writeHistogram) {
							int readSeqID = wordsGenomes[wordGenome_it2->first + genomeCounter1].seqID;
							int genomeSeqID = wordsGenomes[wordGenome_it1->first + genomeCounter2].seqID;
							histogramFile << readSeqID << "\t" << genomeSeqID << "\t" << score << std::endl;
						}

						if (score > fswm_params::g_filteringThreshold) {
							count++;
							int readSeqID = wordsGenomes[wordGenome_it2->first + genomeCounter1].seqID;
							int genomeSeqID = wordsGenomes[wordGenome_it1->first + genomeCounter2].seqID;
							// std::cout << readSeqID << "\t" << genomeSeqID << "\t" << wordsGenomes[wordRead_it->first + readCounter].seqPos << std::endl;
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
				wordGenome_it1++;
				wordGenome_it2++;
			}
		}
		if (fswm_params::g_verbose) { std::cout << "\t\t# matches: " << count << std::endl; }
	}

	if (fswm_params::g_writeHistogram) { histogramFile.close(); }

	return true;
}
