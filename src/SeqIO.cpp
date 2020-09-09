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
void SeqIO::read_sequences(std::string fastafname, std::vector<Sequence> &sequences) {
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
        sequences.push_back(Sequence(header, line, SeqIO::seqID_counter));
    }
    fastafstream.close();
}
