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
#ifndef FSWM_SUBSTITUTIONMATRIX_H_
#define FSWM_SUBSTITUTIONMATRIX_H_

struct SubstitutionMatrix {
	int chiaromonte[4][4] = {{  91, -114,  -31, -123},
                             {-114,  100, -125,  -31},
                             { -31, -125,  100, -114},
                             {-123,  -31, -114,   91}};

    int binary[4][4] = 		{{ 1, -1, -1, -1},
                             {-1,  1, -1, -1},
                             {-1, -1,  1, -1},
                             {-1, -1, -1,  1}};

	int mismatch[4][4]    = {{   0,    1,    1,    1},
                             {   1,    0,    1,    1},
                             {   1,    1,    0,    1},
                             {   1,    1,    1,    0}};

    int transition[4][4]  = {{   0,    0,    1,    0},
                             {   0,    0,    0,    1},
                             {   1,    0,    0,    0},
                             {   0,    1,    0,    0}};

    int transversion[4][4] ={{   0,    1,    0,    1},
                             {   1,    0,    1,    0},
                             {   0,    1,    0,    1},
                             {   1,    0,    1,    0}};
};

#endif
