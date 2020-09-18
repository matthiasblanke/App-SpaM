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
#ifndef FSWM_NODE_H_
#define FSWM_NODE_H_

#include <string>
#include <vector>
#include <iostream>
#include "GlobalParameters.h"

struct Node {
	std::string name;				// Node name should be identical to reference sequence name for leaves
	seq_id_t ID;					// ID is identical to reference sequence IDs for leaves and unique for internal nodes
	Node* parent;
	std::vector<Node*> children;
	scoring_t similarityScore;			// Similarity score represents the similarity of node to current read
	count_t weight; 					// Weight for calculating similarity scores of inner nodes
	count_t leaves_below;
	double distance;				// Distance from father node to this node

	Node(std::string name);

	bool add_child(Node *child);
	bool remove_child(Node *child);

	friend std::ostream& operator<<(std::ostream &strm, const Node &node); 
};

#endif
