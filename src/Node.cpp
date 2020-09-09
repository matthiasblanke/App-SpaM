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
#include "Node.h"
#include "SeqIO.h"

Node::Node(std::string name) {
	this->parent = nullptr;
	this->distance = 0;
	this->similarityScore = -1;
	this->weight = -1;
	this->leaves_below = -1;
	this->name = name;

	if (fswm_internal::namesToSeqIDs.find(name) != fswm_internal::namesToSeqIDs.end()) { // Use existing ID for leaves
		this->ID = fswm_internal::namesToSeqIDs[name];
	}
	else {		// Create new sequence ID for internal nodes that reads can be assigned to
		SeqIO::seqID_counter++;
		this->ID = SeqIO::seqID_counter;
		fswm_internal::seqIDsToNames[this->ID] = name;
		fswm_internal::namesToSeqIDs[this->name] = this->ID;
	}
}

bool Node::add_child(Node *child) {
	children.push_back(child);
	child->parent = this;
    return true;
}

bool Node::remove_child(Node *child) {
	for (int i = 0; i < children.size(); i++) {
		if (children[i] == child) {
			children.erase(children.begin() + i);
			return true;
		}
	}
    return false;
}

std::ostream& operator<<(std::ostream &strm, const Node &node) {
	strm << "Name:" << node.name << "\tDist:" << node.distance << "\tSimS" << node.similarityScore << "\tWeight:" << node.weight << "\tNodesBelow:" << node.leaves_below << "\tID:" << node.ID;
	return strm;
}
