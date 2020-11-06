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
#include <string>
#include "BucketManager.h"
#include "Seed.h"
#include "ReadManager.h"
#include "GenomeManager.h"
#include "Algorithms.h"
#include "Scoring.h"
#include "GlobalParameters.h"
#include "Placement.h"
#include <omp.h>

int main(int argc, char *argv[]) {
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << " Alignment-free phylogenetic placement algorithm" << std::endl;
	std::cout << "           based on spaced word matches         " << std::endl;
	std::cout << "------------------------------------------------" << std::endl << std::endl;

	// Parse command line options and check for correctness
	GlobalParameters::parse_parameters(argc,  argv);
	GlobalParameters::check_parameters();
	if (fswm_params::g_verbose) {GlobalParameters::print_to_console(); };

	omp_set_dynamic(0);
	omp_set_num_threads(fswm_params::g_threads);

	Placement::phylogenetic_placement();

	if (fswm_params::g_writeParameter) { GlobalParameters::save_parameters(); };

	std::cout << std::endl << "-> Placement finished. Output files are in the folder: "
			  << fswm_params::g_outfoldername << std::endl;
}
