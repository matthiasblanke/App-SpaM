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

#include "Placement.h"
#include "BucketManager.h"
#include "GenomeManager.h"
#include "ReadManager.h"
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

void Placement::phylogenetic_placement() {
	// Initialize pattern
	Pattern pattern = Pattern(NULL, NULL, fswm_params::g_numPatterns, fswm_params::g_weight + fswm_params::g_spaces,
			fswm_params::g_weight, 10000, 10000, 10000, 0.75, 0.25, 0); // Create pattern with seed

	pattern.Silent();
	pattern.ImproveSecure();
	pattern.Improve(10);

	//pattern.Print();
	std::vector<std::string> patterns = pattern.GetPattern();

	std::cout << "pattern size : " << patterns.size() << std::endl;

	for (int i = 0; i < fswm_params::g_numPatterns; i++) {
		Seed seed(fswm_params::g_weight, fswm_params::g_spaces);
		seed.generate_pattern(patterns[i]);
		fswm_internal::seeds.push_back(seed);
	}

	// Read reads
	ReadManager	readManager(fswm_params::g_readsfname, fswm_internal::seeds);

	// Read genomes, create spaced words and organize BucketManagers
	GenomeManager genomeManager(fswm_params::g_genomesfname, fswm_internal::seeds);

	// Create empty output files
	Placement::create_output_files();

	// Compare buckets of reads and genomes
	std::cout << "-> Compare reads and genomes." << std::endl;
	// #pragma omp parallel for

	Scoring fswm_distances = Scoring();

	BucketManager bucketManagerGenomes = genomeManager.get_BucketManager();
	BucketManager bucketManagerReads = readManager.get_BucketManager();
	Algorithms::fswm_complete(bucketManagerGenomes, bucketManagerReads, fswm_distances);

	{
		std::cout << "-> Calculate distances." << std::endl;
		fswm_distances.calculate_fswm_distances();
		std::cout << "-> Placing reads in reference tree." << std::endl;
		fswm_distances.phylogenetic_placement();

		if (fswm_params::g_writeScoring or fswm_params::g_assignmentMode == "APPLES") {
			fswm_distances.write_scoring_to_file();
			fswm_distances.write_scoring_to_file_as_table();
		}
	}

	if (fswm_params::g_assignmentMode == "APPLES") {
		int res = std::system(("run_apples.py -t " + fswm_params::g_reftreefname + " -d " + fswm_params::g_outfoldername +
					"scoring_table.txt -o " + fswm_params::g_outfoldername + fswm_params::g_outjplacename).c_str());
		if (res != 0) {
			std::cout << "Apples did not run properly." << std::endl;
			exit (EXIT_FAILURE);
		}
	}
}

void Placement::create_output_files() {
	std::ofstream assignmentJPlaceFile;
	assignmentJPlaceFile.open(fswm_params::g_outfoldername + fswm_params::g_outjplacename);
	assignmentJPlaceFile.close();

	if (fswm_params::g_writeHistogram) {
		std::ofstream histogramFile;
		histogramFile.open(fswm_params::g_outfoldername + "histogram.txt");
		histogramFile.close();
	}
}
