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
#include <fstream>
#include <iostream>
#include "SeqIO.h"

seq_id_t SeqIO::seqID_counter = -1;

// Read sequence from fasta file 'filename'
void SeqIO::read_sequences(std::string fastafname, bool isQuery, BucketManager &bucketManager) {
    std::ifstream fastafstream(fastafname);
    std::string line;
    std::string header;
    std::string last_header = header;

    std::getline(fastafstream, line, '>');
    while (!fastafstream.eof()) {
        std::getline(fastafstream, header);
        SeqIO::seqID_counter++;
        header = header.substr(0, header.find(' '));
        last_header = header;
        std::getline(fastafstream, line, '>');
        parse_sequence(header, line, SeqIO::seqID_counter, isQuery, bucketManager);
    }
    fastafstream.close();
}

void SeqIO::parse_sequence(std::string header, std::string seqLine, seq_id_t seqID, bool isQuery, BucketManager &bucketManager) {
    std::vector<char> seq;
    seq.reserve(seqLine.size()*2);

    // Create normal sequence and its reverse
    for (std::string::iterator it = seqLine.begin(); it != seqLine.end(); it++) {
        if (!std::isspace(*it)) {
            char c = std::toupper(*it);
            switch(c) {
                case 'A': seq.push_back(0x00); break;
                case 'C': seq.push_back(0x01); break;
                case 'G': seq.push_back(0x02); break;
                case 'T': seq.push_back(0x03); break;
                case 'U': seq.push_back(0x03); break;
                default:  break;
            }
        }
    }

    for (std::string::reverse_iterator rit = seqLine.rbegin(); rit != seqLine.rend(); rit++) {
        if (!std::isspace(*rit)) {
            char c = std::toupper(*rit);
            switch(c) {
                case 'A': seq.push_back(0x03); break;
                case 'C': seq.push_back(0x02); break;
                case 'G': seq.push_back(0x01); break;
                case 'T': seq.push_back(0x00); break;
                case 'U': seq.push_back(0x00); break;
                default: break;
            }
        }
    }

    /*
    if (isQuery) {
        fswm_internal::querySequences[seqID] = seq;
    }
    else {
        fswm_internal::referenceSequences[seqID] = seq;
    }*/
    
    const size_t NumBytes = 8;

    for (auto &seed : fswm_internal::seeds) {
        std::vector<int> matchPos = seed.get_matchPos();
        std::vector<int> dontCarePos = seed.get_dontCarePos();

        // Go through all spaced words in sequence and save to bucket
        uint32_t go_until = std::max(int(seq.size() - fswm_params::g_weight - fswm_params::g_spaces + 1), 0);
        for (uint32_t i = 0; i < go_until; i++) {
            // Create spaced word at position i and give word to BucketManager for further processing
            word_t matches = 0;
            word_t dontCares = 0;
            for (auto const &pos : matchPos) {
                matches = matches << 2;
                matches += seq[i+pos];
            }

            for (auto const &pos : dontCarePos) {
                dontCares = dontCares << 2;
                dontCares += seq[i+pos];
            }
    
            if (fswm_params::g_sampling) {
                auto res = crc32_fast(&matches, NumBytes);
                if (res < fswm_params::g_minHashLowerLimit) {
                    bucketManager.insert_word(matches, dontCares, seqID);
                }
            }
            else {
                bucketManager.insert_word(matches, dontCares, seqID);
            }
        }
    }

    if (fswm_internal::seqIDsToNames.find(seqID) != fswm_internal::seqIDsToNames.end() or fswm_internal::namesToSeqIDs.find(header) != fswm_internal::namesToSeqIDs.end()) {
        std::cerr << "Multiple sequences in the genomes seem to have the same name. Please fix." << std::endl;
        exit(EXIT_FAILURE);
    }
    fswm_internal::seqIDsToNames[seqID] = header;
    fswm_internal::namesToSeqIDs[header] = seqID;
    fswm_internal::genomeIDsToNames[seqID] = header;
    fswm_internal::namesToGenomeIDs[header] = seqID;
}