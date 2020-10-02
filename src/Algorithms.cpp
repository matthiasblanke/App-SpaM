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
bool Algorithms::fswm_complete(BucketManager &genomeBucketManager, BucketManager &readBucketManager, Scoring &fswm_distances) {
	SubstitutionMatrix substMat;
	std::ofstream histogramFile;
	if (fswm_params::g_writeHistogram) { histogramFile.open(fswm_params::g_outfoldername + "histogram.txt", std::ios_base::app); }

	// Loop through minimizers and compare each bucket on its own
	for (auto const minimizer : genomeBucketManager.get_minimizers()) {
		//Get buckets
		Bucket bucketGenomes = genomeBucketManager.get_bucket(minimizer);
		Bucket bucketReads = readBucketManager.get_bucket(minimizer);

		// Loop through buckets and compare spaced words
		std::vector<Word> wordsGenomes = bucketGenomes.get_words();
		std::vector<Word> wordsReads = bucketReads.get_words();

		// Get vector of word groups. First int is starting position, second int length of group
		std::vector<std::pair<int,int>> wordGroupGenomes = bucketGenomes.get_wordGroups();
		std::vector<std::pair<int,int>> wordGroupReads = bucketReads.get_wordGroups();

		std::vector<std::pair<int,int>>::const_iterator wordGenome_it = wordGroupGenomes.cbegin();
		std::vector<std::pair<int,int>>::const_iterator wordRead_it = wordGroupReads.cbegin();

		if (fswm_params::g_verbose) {
			std::cout << "\tBucket: " << bucketGenomes.get_minimizer() << std::endl;
			std::cout << "\t\tBucket size genomes: " << bucketGenomes.get_bucketSize() << std::endl;
			std::cout << "\t\tBucket size reads: " << bucketReads.get_bucketSize() << std::endl;
		}

		word_t dontCaresGenome;
		word_t dontCaresRead;
		int count = 0;
		// Loop through all word groups
		while (wordRead_it != wordGroupReads.end() and wordGenome_it != wordGroupGenomes.end()) {
			if (wordsGenomes[wordGenome_it->first] < wordsReads[wordRead_it->first]) {
				wordGenome_it++;
			}
			else if (wordsGenomes[wordGenome_it->first] > wordsReads[wordRead_it->first]) {
				wordRead_it++;
			}
			else {
				// Loop through words with same spaced k-mer
				for (int readCounter = 0; readCounter < wordRead_it->second; readCounter++) {
					for (int genomeCounter = 0; genomeCounter < wordGenome_it->second; genomeCounter++) {

						// For each match calculate spaced word score
						dontCaresGenome = wordsGenomes[wordGenome_it->first + genomeCounter].dontCares;
						dontCaresRead = wordsReads[wordRead_it->first + readCounter].dontCares;

						int score = 0;
						int mismatches = 0;
						int transitions = 0;
						int transversions = 0;

						for (int i = 0; i < fswm_params::g_spaces; i++) {
							score += substMat.chiaromonte[(dontCaresGenome & 0x03)][(dontCaresRead & 0x03)];
							mismatches += substMat.mismatch[(dontCaresGenome & 0x03)][(dontCaresRead & 0x03)];
							transitions += substMat.transition[(dontCaresGenome & 0x03)][(dontCaresRead & 0x03)];
							transversions += substMat.transversion[(dontCaresGenome & 0x03)][(dontCaresRead & 0x03)];
							dontCaresRead = dontCaresRead >> 2;
							dontCaresGenome = dontCaresGenome >> 2;
						}

						if (fswm_params::g_writeHistogram) {
							int readSeqID = wordsReads[wordRead_it->first + readCounter].seqID;
							int genomeSeqID = wordsGenomes[wordGenome_it->first + genomeCounter].seqID;
							histogramFile << readSeqID << "\t" << genomeSeqID << "\t" << score << std::endl;
						}

						if (score > fswm_params::g_filteringThreshold) {
							count++;
							int readSeqID = wordsReads[wordRead_it->first + readCounter].seqID;
							int genomeSeqID = wordsGenomes[wordGenome_it->first + genomeCounter].seqID;
							// std::cout << readSeqID << "\t" << genomeSeqID << "\t" << wordsReads[wordRead_it->first + readCounter].seqPos << std::endl;
							if (fswm_distances.scoringMap.find(readSeqID) == fswm_distances.scoringMap.end()) {
								fswm_distances.scoringMap[readSeqID] = std::unordered_map<seq_id_t, scoring_t>();
								fswm_distances.mismatchCount[readSeqID] = std::unordered_map<seq_id_t, count_t>();
								fswm_distances.transitionCount[readSeqID] = std::unordered_map<seq_id_t, count_t>();
								fswm_distances.transversionCount[readSeqID] = std::unordered_map<seq_id_t, count_t>();
								fswm_distances.spacedWordMatchCount[readSeqID] = std::unordered_map<seq_id_t, count_t>();
							}
							fswm_distances.scoringMap[readSeqID][genomeSeqID] += score;
							fswm_distances.mismatchCount[readSeqID][genomeSeqID] += mismatches;
							fswm_distances.transitionCount[readSeqID][genomeSeqID] += transitions;
							fswm_distances.transversionCount[readSeqID][genomeSeqID] += transversions;
							fswm_distances.spacedWordMatchCount[readSeqID][genomeSeqID] += 1;
						}
					}
				}
				wordGenome_it++;
				wordRead_it++;
			}
		}
		if (fswm_params::g_verbose) { std::cout << "\t\t# matches: " << count << std::endl; }
	}

	if (fswm_params::g_writeHistogram) { histogramFile.close(); }

	return true;
}
