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

#include "PairwiseDistances.h"
#include "BucketManager.h"
#include "GenomeManager.h"
#include "ReadManager.h"
#include "Sequence.h"
#include "Algorithms.h"
#include "Tree.h"
#include "GlobalParameters.h"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Algorithms.h"
#include "Scoring.h"
#include "SubstitutionMatrix.h"
#include "Match.h"
#include "MatchManager.h"

void PairwiseDistances::pairwise_distances() {
	// Initialize pattern
	Pattern pattern = Pattern(fswm_params::g_numPatterns, fswm_params::g_weight + fswm_params::g_spaces,
			fswm_params::g_weight, 0); // Create pattern with seed

	pattern.Silent();
	pattern.ImproveSecure();
	pattern.Improve(10);

	//pattern.Print();
	std::vector<std::string> patterns = pattern.GetPattern();

	if (fswm_params::g_verbose) { std::cout << "-> Pattern size : " << patterns.size() << std::endl; }

	std::vector<Seed> seeds;
	for (int i = 0; i < fswm_params::g_numPatterns; i++) {
		Seed seed(fswm_params::g_weight, fswm_params::g_spaces);
		seed.generate_pattern(patterns[i]);
		seeds.push_back(seed);
	}

	std::cout << "-> Reading sequences." << std::endl;
	// Read reads
	GenomeManager genomeManager1(fswm_params::g_genomesfname, seeds, "_1");

	// Read genomes, create spaced words and organize BucketManagers
	GenomeManager genomeManager2(fswm_params::g_genomesfname, seeds, "_2");

	if (fswm_params::g_writeIDs) { GlobalParameters::write_read_ids_to_file(); };
	if (fswm_params::g_writeIDs) { GlobalParameters::write_seq_ids_to_file(); };

	// Create empty output files
	PairwiseDistances::create_output_files();

	if (fswm_params::g_writeScoring) {
		std::ofstream results;
		results.open(fswm_params::g_outfoldername + "scoring_table.txt");

		for (auto genome : fswm_internal::genomeIDsToNames) {		// Write columns names (references to file)
			results << "\t" << genome.second;
		}
		results << std::endl;
		results.close();

		results.open(fswm_params::g_outfoldername + "scoring_list.txt");
		results.close();
	}

	// Compare buckets of reads and genomes
	std::cout << "-> Comparing sets of genomes." << std::endl;

	Scoring fswm_distances = Scoring();

	BucketManager bucketManagerGenomes1 = genomeManager1.get_BucketManager();

	BucketManager bucketManagerGenomes2 = genomeManager2.get_BucketManager();
	
	Algorithms::fswm_complete(bucketManagerGenomes1, bucketManagerGenomes2, fswm_distances);

	fswm_distances.calculate_fswm_distances();

	if (fswm_params::g_writeScoring) {
		fswm_distances.write_scoring_to_file();
		fswm_distances.write_scoring_to_file_as_table();
	}
}

void PairwiseDistances::create_output_files() {
	if (fswm_params::g_writeHistogram) {
		std::ofstream histogramFile;
		histogramFile.open(fswm_params::g_outfoldername + "histogram.txt");
		histogramFile.close();
	}
}
