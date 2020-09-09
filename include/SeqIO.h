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

#ifndef FSWM_SEQIO_H_
#define FSWM_SEQIO_H_

#include "Sequence.h"

class SeqIO {
    public:
    	static seq_id_t seqID_counter;
        static void read_sequences(std::string fastafname, std::vector<Sequence> &sequences);
       	static void reset_seqID_counter();
};

inline void SeqIO::reset_seqID_counter() {
	seqID_counter = -1;
}

#endif
