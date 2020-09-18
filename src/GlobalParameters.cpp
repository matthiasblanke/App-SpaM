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
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>
#include "GlobalParameters.h"

// Initialize global parameters to and set default values

// IO parameters
std::string fswm_params::g_genomesfname = "";
std::string fswm_params::g_reftreefname = "";
std::string fswm_params::g_readsfname = "";
std::string fswm_params::g_outjplacename = "appspam_placement_results.jplace";
std::string fswm_params::g_outfoldername = "./";
std::string fswm_params::g_paramfname = "";

// General parameters
uint16_t fswm_params::g_weight = 12;
uint16_t fswm_params::g_spaces = 32;
std::string fswm_params::g_assignmentMode = "LCACOUNT";
std::string fswm_params::phylo_placement_mode = "ALLBEST";
bool fswm_params::g_verbose = true;
int fswm_params::g_filteringThreshold = GlobalParameters::calculate_filteringThreshold();
int fswm_params::g_filteringThresholdMultiplicator = 0;
bool fswm_params::g_sampling = false;
int fswm_params::g_minHashLowerLimit = 10000;

// Additional options
uint16_t fswm_params::g_threads = 1;
uint32_t fswm_params::g_readBlockSize = 1000;
bool fswm_params::g_writeHistogram = false;
bool fswm_params::g_writeScoring = false;
int fswm_params::g_allowance = 0;
double fswm_params::default_distance_new_leaves = 0.001;
int fswm_params::g_numPatterns = 10;
double fswm_params::g_defaultDistance = 10;

// Initialize global internal mappings between sequence IDs and names.
std::unordered_map<seq_id_t, std::string> fswm_internal::seqIDsToNames = std::unordered_map<seq_id_t, std::string>();
std::unordered_map<std::string, seq_id_t> fswm_internal::namesToSeqIDs = std::unordered_map<std::string, seq_id_t>();

// Only mappings between reference sequences and ids
std::unordered_map<seq_id_t, std::string> fswm_internal::genomeIDsToNames = std::unordered_map<seq_id_t, std::string>();
std::unordered_map<std::string, seq_id_t> fswm_internal::namesToGenomeIDs = std::unordered_map<std::string, seq_id_t>();

// Mappings between query sequences and ids
std::unordered_map<seq_id_t, std::string> fswm_internal::readIDsToNames = std::unordered_map<seq_id_t, std::string>();
std::unordered_map<std::string, seq_id_t> fswm_internal::namesToReadIDs = std::unordered_map<std::string, seq_id_t>();
std::unordered_map<seq_id_t, seq_id_t> fswm_internal::IDsToPlacementIDs = std::unordered_map<seq_id_t, seq_id_t>();
std::unordered_map<seq_id_t, seq_id_t> fswm_internal::placementIDsToIDs = std::unordered_map<seq_id_t, seq_id_t>();

std::unordered_map<seq_id_t, std::vector<char>> fswm_internal::referenceSequences = std::unordered_map<seq_id_t, std::vector<char>>();
std::unordered_map<seq_id_t, std::vector<char>> fswm_internal::querySequences = std::unordered_map<seq_id_t, std::vector<char>>();

std::vector<Seed> fswm_internal::seeds;

int fswm_internal::g_numberGenomes = 0;

bool GlobalParameters::save_parameters() {
	std::ofstream foutstream(fswm_params::g_outfoldername + "fswm_parameters.txt");
	foutstream << "  Parameters : {" << std::endl;
	foutstream << "\treference : " << fswm_params::g_genomesfname << "," << std::endl;
	foutstream << "\ttree : " << fswm_params::g_reftreefname << "," << std::endl;
	foutstream << "\tquery : " << fswm_params::g_readsfname << "," << std::endl;
	foutstream << "\tout_jplace : " << fswm_params::g_outjplacename << "," << std::endl;
	foutstream << "\tweight : " << fswm_params::g_weight << "," << std::endl;
	foutstream << "\tspaces : " << fswm_params::g_spaces << "," << std::endl;
	foutstream << "\tallowance :" << fswm_params::g_allowance << "," << std::endl;
	foutstream << "\tassignment_mode : " << fswm_params::g_assignmentMode << "," << std::endl;
	foutstream << "\tread_block_size : " << fswm_params::g_readBlockSize << "," << std::endl;
	foutstream << "  }" << std::endl << "}" << std::endl;
	foutstream.close();

	return true;
}

bool GlobalParameters::load_parameters(std::string filename = "") {
	if (filename == "") {
		return false;
	}

	std::ifstream finstream(filename);
	if (!finstream.is_open()) {
		return false;
	}

	std::size_t found = filename.find_last_of("/");
	std::string param_folder = filename.substr(0,found) + "/";

	std::string line;
	std::string delimiter = " : ";
	while (std::getline(finstream, line)) {
		if (line.find(delimiter) != std::string::npos) {
			std::string key = line.substr(0, line.find(delimiter));
			std::string value = line.substr(line.find(delimiter)+delimiter.size(), line.size() - line.find(delimiter) - delimiter.size() - 1);
			if (key.find("weight") != std::string::npos) {
				fswm_params::g_weight = std::stoi(value);
			}
			if (key.find("spaces") != std::string::npos) {
				fswm_params::g_spaces = std::stoi(value);
				fswm_params::g_filteringThreshold = calculate_filteringThreshold();
			}
			if (key.find("allowance") != std::string::npos) {
				fswm_params::g_allowance = std::stoi(value); 					//stoi: string to int
			}
			if (key.find("assignment_mode") != std::string::npos) {
				fswm_params::g_assignmentMode = value;
			}
			if (key.find("threads") != std::string::npos) {
				fswm_params::g_threads = std::stoi(value);
			}
			if (key.find("read_block_size") != std::string::npos) {
				fswm_params::g_readBlockSize = std::stoi(value);
			}
			if (key.find("verbose") != std::string::npos) {
				fswm_params::g_verbose = std::stoi(value);
			}
			if (key.find("threshold") != std::string::npos) {
				fswm_params::g_filteringThresholdMultiplicator = std::stoi(value);
				calculate_filteringThreshold();
			}
			if (key.find("reference") != std::string::npos) {
				if (value.rfind("/", 0) == 0) {
					fswm_params::g_genomesfname = value;
				}
				else {
					fswm_params::g_genomesfname = param_folder + value;
				}
			}
			if (key.find("query") != std::string::npos) {
				if (value.rfind("/", 0) == 0) {
					fswm_params::g_readsfname = value;
				}
				else {
					fswm_params::g_readsfname = param_folder + value;
				}
			}
			if (key.find("tree") != std::string::npos) {
				if (value.rfind("/", 0) == 0) {
					fswm_params::g_reftreefname = value;
				}
				else {
					fswm_params::g_reftreefname = param_folder + value;
				}
			}
			if (key.find("out_jplace") != std::string::npos) {
				fswm_params::g_outjplacename = value;
				if (fswm_params::g_outjplacename.find('/') != std::string::npos) {
    				fswm_params::g_outfoldername = fswm_params::g_outjplacename.substr(0,fswm_params::g_outjplacename.find_last_of("/")) + "/";
					fswm_params::g_outjplacename = fswm_params::g_outjplacename.substr(fswm_params::g_outjplacename.find_last_of("/") + 1);
				}
			}
		}
	}

	finstream.close();
	return true;
}

/** Parse option parameters from parameter file or command line. */
bool GlobalParameters::parse_parameters(int argc, char *argv[]) {
	int option_param;
	std::string possible_params = "l:s:t:q:o:k:d:hm:g:a:b:vi";
	bool usingParameterfile = false;

    int index = -1;
	static const struct option long_options[] =
    {
    	{ "load", required_argument,      	 	nullptr, 'l' },
        { "reference", required_argument,       nullptr, 's' },
        { "tree", required_argument,			nullptr, 't' },
        { "query", required_argument,  			nullptr, 'q' },
        { "out_jplace", required_argument, 		nullptr, 'o' },
        { "weight", required_argument, 			nullptr, 'k' },
        { "dontCare", required_argument,       	nullptr, 'd' },
        { "threads", required_argument, 		nullptr, 1   },
        { "help", no_argument, 					nullptr, 'h' },
        { "assignment_mode", required_argument, nullptr, 'g' },
		{ "allowance", required_argument,		nullptr, 'a' },
        { "read_block_size", required_argument, nullptr, 'b' },
        { "verbose", no_argument, 				nullptr, 'v' },
        { "histogram", no_argument, 			nullptr, 2   },
        { "scoring", no_argument, 				nullptr, 3   },
        { "sampling", no_argument, 				nullptr, 4   },
        { "threshold", required_argument, 		nullptr, 5   },
        0
    };

	// Scan for parameter file first and load parameters from file
	while ((option_param = getopt_long(argc, argv, possible_params.c_str(), long_options, &index)) != -1) {
		switch (option_param) {
			case 'l':
				usingParameterfile = true;
				fswm_params::g_paramfname = optarg;
				if (!GlobalParameters::load_parameters(fswm_params::g_paramfname)) {
					std::cout << "Please supply an existing parameter file." << std::endl;
					exit (EXIT_FAILURE);
				}
				break;
		}
	}
	optind = 1;

	// If no parameters are given
	if (argc == 1) {
		print_help();
		exit (EXIT_FAILURE);
	}

    while ((option_param = getopt_long(argc, argv, possible_params.c_str(), long_options, &index)) != -1) {
		switch (option_param) {
			case 's':
				fswm_params::g_genomesfname = optarg;
				break;
			case 'q':
				fswm_params::g_readsfname = optarg;
				break;
			case 'o':
				fswm_params::g_outjplacename = optarg;
				if (fswm_params::g_outjplacename.find('/') != std::string::npos) {
    				fswm_params::g_outfoldername = fswm_params::g_outjplacename.substr(0,fswm_params::g_outjplacename.find_last_of("/")) + "/";
					fswm_params::g_outjplacename = fswm_params::g_outjplacename.substr(fswm_params::g_outjplacename.find_last_of("/") + 1);
				}
				break;
			case 't':
				fswm_params::g_reftreefname = optarg;
				break;
			case 'k':
				fswm_params::g_weight = atoi(optarg);
				break;
			case 'd':
				fswm_params::g_spaces = atoi(optarg);
				fswm_params::g_filteringThreshold = calculate_filteringThreshold();
				break;
			case 1:
				fswm_params::g_threads = atoi(optarg);
				break;
			case 'h':
				print_help();
				exit (EXIT_FAILURE);
				break;
			case 'b':
				fswm_params::g_readBlockSize = atoi(optarg);
				break;
			case 'v':
				fswm_params::g_verbose = true;
				break;
			case 2:
				fswm_params::g_writeHistogram = true;
				break;
			case 3:
				fswm_params::g_writeScoring = true;
				break;
			case 'g':
				fswm_params::g_assignmentMode = optarg;
				break;
			case 4:
				fswm_params::g_sampling = true;
				break;
			case 5:
				fswm_params::g_filteringThresholdMultiplicator = atoi(optarg);
				calculate_filteringThreshold();
				break;
			case 'a':
				fswm_params::g_allowance = atoi(optarg);				//atoi: character string to int
				break;
			case '?':
				print_help();
				exit (EXIT_FAILURE);
      	}
	}

	return true;
}

bool GlobalParameters::check_parameters() {
	if(fswm_params::g_weight < 2 || fswm_params::g_weight > 32) {
		std::cerr << "Weight (-k) must be between 2 and 32"<< std::endl;
		print_to_console();
		exit (EXIT_FAILURE);
	}
	if(fswm_params::g_spaces < 2 || fswm_params::g_spaces > 32) {
		std::cerr << "SPACES (-d) must be between 2 and 32"<< std::endl;
		print_to_console();
		exit (EXIT_FAILURE);
	}
	/*if (fswm_params::g_assignmentMode != "BEST" and fswm_params::g_assignmentMode != "LCA" and fswm_params::g_assignmentMode != "RANDOM" and fswm_params::g_assignmentMode != "PHYLOKMERS" and fswm_params::g_assignmentMode != "APPLES") {
		std::cerr << "AssignmentMode must be \"BEST\" or \"LCA\" or \"RANDOM\" or \"PHYLOKMERS\" or \"APPLES\"."<< std::endl;
		print_to_console();
		exit (EXIT_FAILURE);
	}*/
	if (fswm_params::g_assignmentMode.rfind("PHYLOKMERS", 0) == 0) {
		std::stringstream test;
		test.str (fswm_params::g_assignmentMode);
		std::string segment;
		std::vector<std::string> seglist;

		while(std::getline(test, segment, 'X')) {
   			seglist.push_back(segment);
		}
		if (seglist.size() != 2) {
			std::cerr << "PHYLOKMERS placement must be followed by mode." << std::endl;
			print_to_console();
			exit (EXIT_FAILURE);
		}
		fswm_params::g_assignmentMode = seglist[0];
		fswm_params::phylo_placement_mode = seglist[1];
	}
	if (fswm_params::g_allowance < 0 ) {
		std::cerr << "Allowance must be an integer >= 0." << std::endl;
		print_to_console();
		exit (EXIT_FAILURE);
	}
	if (fswm_params::g_readBlockSize < 1 or fswm_params::g_readBlockSize > 200000) {
		std::cerr << "Choose a block size between 1 and 200000."<< std::endl;
		print_to_console();
		exit (EXIT_FAILURE);
	}
	if(fswm_params::g_threads < 1) {
		std::cerr << "Threads (-t) must be an integer larger than 0."<< std::endl;
		print_to_console();
		exit (EXIT_FAILURE);
	}
	std::ifstream f(fswm_params::g_genomesfname.c_str());
	if (!f.good()) {
		std::cout << "Please supply an existing file for the genomes." << std::endl;
		print_to_console();
		print_help();
		exit (EXIT_FAILURE);
	}

	std::ifstream g(fswm_params::g_readsfname.c_str());
	if (!g.good()) {
		std::cout << "Please supply an existing file for the reads." << std::endl;
		print_to_console();
		print_help();
		exit (EXIT_FAILURE);
	}

	std::ifstream h(fswm_params::g_reftreefname.c_str());
	if (!h.good() and !(fswm_params::g_reftreefname == "not set")) {
		std::cout << "Please supply an existing file for reference tree." << std::endl;
		print_to_console();
		print_help();
		exit (EXIT_FAILURE);
	}
	return true;
}

bool GlobalParameters::print_to_console() {
	std::cout << "Parameters:" << std::endl;
	std::cout << "\tweight  : " << fswm_params::g_weight << std::endl;
	std::cout << "\tspaces  : " << fswm_params::g_spaces << std::endl;
	std::cout << "\tallowance : " << fswm_params::g_allowance << std::endl;
	std::cout << "\tthreads : " << fswm_params::g_threads << std::endl;
	std::cout << "\tassignment : " << fswm_params::g_assignmentMode << std::endl;
	std::cout << "\tread_block_size  : " << fswm_params::g_readBlockSize << std::endl;
	std::cout << "\tVerbose  : " << fswm_params::g_verbose << std::endl;
	std::cout << "\treference  : " << fswm_params::g_genomesfname << std::endl;
	std::cout << "\tquery  : " << fswm_params::g_readsfname << std::endl;
	std::cout << "\ttree  : " << fswm_params::g_reftreefname << std::endl;
	std::cout << "\tout_jplace  : " << fswm_params::g_outjplacename << std::endl;
	std::cout << "\tout_folder  : " << fswm_params::g_outfoldername << std::endl;
	return true;
}

void GlobalParameters::write_genome_ids_to_file() {
	std::ofstream genomeIDsStream(fswm_params::g_outfoldername + "genomeIDsToNames.txt");
	for (auto const &entry : fswm_internal::genomeIDsToNames) {
		genomeIDsStream << entry.first << "\t" << entry.second << std::endl;
	}
	genomeIDsStream.close();
}

void GlobalParameters::write_seq_ids_to_file() {
	std::ofstream seqIDsOutStream(fswm_params::g_outfoldername + "namesToSeqIDs.txt");
	for (auto const &entry : fswm_internal::namesToSeqIDs) {
		seqIDsOutStream << entry.first << "\t" << entry.second << std::endl;
	}
	seqIDsOutStream.close();
}

void GlobalParameters::write_read_ids_to_file() {
	std::ofstream seqIDsOutStream(fswm_params::g_outfoldername + "readsToSeqIDs.txt");
	for (auto const &entry : fswm_internal::readIDsToNames) {
		seqIDsOutStream << entry.first << "\t" << entry.second << std::endl;
	}
	seqIDsOutStream.close();
}

int GlobalParameters::calculate_filteringThreshold() {
	fswm_params::g_filteringThreshold = fswm_params::g_spaces * fswm_params::g_filteringThresholdMultiplicator;
	return fswm_params::g_spaces * fswm_params::g_filteringThresholdMultiplicator;
}

/**
 * Print helpfile and exit.
 */
void GlobalParameters::print_help() {
	std::cout << R""""(
------------------------------------------------------------

Execute appspam with:
	./appspam [parameters]

------------------------------------------------------------

A typical call might look like:
	./appspam -l appspam_parameters.txt
	./appspam -h
	./appspam -s references.fasta -q query.fasta -t tree.nwk
	./appspam -l fswm_parameters.txt -d 10 -w 8

The following parameters are necessary:
	-s 	Reference sequences.
		Full path of fasta file with references.
		(Normally marker genes or fully assembled genomes,
		i.e. every sequence is corresponding to one reference.
		If unassembled references are used, the corresponding flag
		and delimiter parameters need to be used.)
	-q 	Query sequences.
		Full path of fasta file with query sequences.
	-t	Reference tree.
		File of reference tree in newick format.
		(Rooted, bifurcating tree in newick format.
		All leaves must have identical names to reference genomes.)

The following parameters are optional.
	--outputFolder
	-o 	Output folder				./results/
		Path to folder for all output files

	--weight
	-w	Weight of pattern.			12
		Number of match positions in pattern. Pattern will be optimized
		automatically.

	--dontCare
	-d	Don't care positions. 		16
		Number of don't care positions in pattern.
		Match positions (-k) plus don't care positions (-d) equal the
		overall pattern length.

	--threads
	 	Number of threads.			1
		The calculations are not parallelized at the moment.

	--mode
	-m	Matching-mode 					0
		Switch between different modes of fswm, i.e. which spaced words are
		compared to each other. Possible values:
		all:		 	"Comare all spaced words against all others"
		first:	 		"Compare only each spaced word only against a random other spaced word match"
		best:		    "Compare each spaced word only against highest scoring match"

	--samling
		Samples the spaced word matches.

	--readBlockSize
	-b	Read block size				100000
		Number of reads that are processed in one chunk.
		Every chunk is read in and processed consecutively.

Following additional flags exist:
	-h	Print out help and exit.
	-v 	Turn on verbose mode with additional information printed to std_out.

)"""";
}
