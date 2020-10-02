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
#ifndef FSWM_GLOBAL_H_
#define FSWM_GLOBAL_H_

#include <unordered_map>
#include <vector>
#include <string>

// Spaced words and dont care positions are encoded in 64bit (2 bits per base).
// Thus the number of match and don't care positions cannot exceed 32 characters.
typedef uint64_t word_t;

// Minimizers are encoded with 16 characters as maximum.
typedef uint32_t minimizer_t;

// Position of spaced word occurences in sequences.
typedef uint32_t pos_t;

// Each input sequence has its own internal id.
typedef uint32_t seq_id_t;

typedef double scoring_t;
typedef int32_t count_t;

// Parameters
namespace fswm_params {
	// Weight of spaced words (k)
	extern uint16_t g_weight;

	// Number of don't care positions of spaced word
	extern uint16_t g_spaces;

	// Defines methods with which reads are assigned or placed in tree
	extern std::string g_assignmentMode;

	// Number of threads used
	extern uint16_t g_threads;

	// Number of reads that are processed at once
	extern uint32_t g_readBlockSize;

	// Toggles verbose mode with additional comments to std::cout
	extern bool g_verbose;

	// Toggles if a histogram of spaced word matches is written to file
	extern bool g_writeHistogram;

	// Toggles if scoring list and table are written to files
	extern bool g_writeScoring;

	// Specifies minimum score (filtering threshold) that determines if a spaced word match is considered homologous.
	extern int g_filteringThreshold;

	// Specifies multiplicator to calculate filtering threshold.
	extern int g_filteringThresholdMultiplicator;

	// Switch turns sampling on or off.
	extern bool g_sampling;

	// When sampling is on, only consider spaced words with crc32(word) < g_minHashLowerLimit.
	extern int g_minHashLowerLimit;

	extern std::string g_delimiter;
	extern bool g_draftGenomes;

	//
	extern int g_numPatterns;
	extern double g_defaultDistance;

	// Full file names of input files
	extern std::string g_genomesfname;
	extern std::string g_reftreefname;
	extern std::string g_readsfname;
	extern std::string g_outjplacename;
	extern std::string g_outfoldername;
	extern std::string g_paramfname;

	// Set default distance of new branch for phylogenetic distance
	extern double default_distance_new_leaves;
}

// Mappings of internally used ids to names of input sequences.
namespace fswm_internal {
	extern std::unordered_map<seq_id_t, std::string> seqIDsToNames;
	extern std::unordered_map<std::string, seq_id_t> namesToSeqIDs;

	extern std::unordered_map<seq_id_t, std::string> genomeIDsToNames;
	extern std::unordered_map<std::string, seq_id_t> namesToGenomeIDs;

	extern std::unordered_map<seq_id_t, std::string> readIDsToNames;
	extern std::unordered_map<std::string, seq_id_t> namesToReadIDs;

	extern std::unordered_map<seq_id_t, seq_id_t> placementIDsToIDs;
	extern std::unordered_map<seq_id_t, seq_id_t> IDsToPlacementIDs;
	extern int g_numberGenomes;
}

class GlobalParameters {
	public:
		// Save parameters to default txt file.
		static bool save_parameters();

		// Load parameters from txt file given as argument.
		static bool load_parameters(std::string filename);

		// Check parameters for correctness.
		static bool check_parameters();

		// Parse command line parameters.
		static bool parse_parameters(int argc, char *argv[]);

		// Print set of parametes to std::out.
		static bool print_to_console();

		// Print help.txt as manual and exit.
		static void print_help();

		// Write mapping of sequence names to their internal ids to file.
		static void write_genome_ids_to_file();
		static void write_seq_ids_to_file();
		static void write_read_ids_to_file();

		// Calculates threshold for spaced word filter.
		static int calculate_filteringThreshold();
};

#endif
