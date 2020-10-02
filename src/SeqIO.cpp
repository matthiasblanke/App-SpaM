/**
 * APPSPAM is a tool for alignment-free phylogenetic placement.
 *  Copyright (C) 2020 Matthias Blanke
 *
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
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include "SeqIO.h"

seq_id_t SeqIO::seqID_counter = -1;

// Read sequence from fasta file 'filename'
void SeqIO::read_sequences(std::string fastafname, std::vector<Sequence> &sequences, bool genomes) {
    std::ifstream fastafstream(fastafname);
    std::string line;
    std::string read;
    std::string header;
    std::string last_header = header;

    if (genomes) {
        std::getline(fastafstream, line, '>');
        while (!fastafstream.eof()) {
            std::getline(fastafstream, header);
            header = header.substr(0, header.find(' '));
            std::getline(fastafstream, line, '>');

            if (!fswm_params::g_draftGenomes) {
                SeqIO::seqID_counter++;
                if (fswm_internal::seqIDsToNames.find(SeqIO::seqID_counter) != fswm_internal::seqIDsToNames.end() or fswm_internal::namesToSeqIDs.find(header) != fswm_internal::namesToSeqIDs.end()) {
                    std::cerr << "Multiple sequences in the genomes seem to have the same name. Please fix: " << header << std::endl;
                    exit(EXIT_FAILURE);
                }
                else {
                    fswm_internal::seqIDsToNames[SeqIO::seqID_counter] = header;
                    fswm_internal::namesToSeqIDs[header] = SeqIO::seqID_counter;
                    fswm_internal::genomeIDsToNames[SeqIO::seqID_counter] = header;
                    fswm_internal::namesToGenomeIDs[header] = SeqIO::seqID_counter;
                    sequences.push_back(Sequence(header, line, SeqIO::seqID_counter));
                    fswm_internal::g_numberGenomes += 1;
                }
            }
            else {
                header = header.substr(0, header.find(fswm_params::g_delimiter));
                if (fswm_internal::namesToSeqIDs.find(header) != fswm_internal::namesToSeqIDs.end()) {
                    sequences.push_back(Sequence(header, line, SeqIO::seqID_counter));
                }
                else {
                    SeqIO::seqID_counter++;
                    fswm_internal::seqIDsToNames[SeqIO::seqID_counter] = header;
                    fswm_internal::namesToSeqIDs[header] = SeqIO::seqID_counter;
                    fswm_internal::genomeIDsToNames[SeqIO::seqID_counter] = header;
                    fswm_internal::namesToGenomeIDs[header] = SeqIO::seqID_counter;
                    sequences.push_back(Sequence(header, line, SeqIO::seqID_counter));
                    fswm_internal::g_numberGenomes += 1;
                }
            }

            
        }
    }
    else {
        std::getline(fastafstream, line, '>');
        while (!fastafstream.eof()) {
            std::getline(fastafstream, header);
            SeqIO::seqID_counter++;
            header = header.substr(0, header.find(' '));
            std::getline(fastafstream, line, '>');
            sequences.push_back(Sequence(header, line, SeqIO::seqID_counter));
        }
    }

    fastafstream.close();
}
