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
#include <fstream>
#include <math.h>
#include <limits>
#include <unordered_set>
#include <string>
#include "Tree.h"
#include "SubstitutionMatrix.h"

/**
 * Create tree from newick file.
 */
Tree::Tree(std::string filename) {
	root = new Node("internal_1");
	internalNodeCounter = 1;
	is_rooted = true;

	std::ifstream newickFile(filename);
	std::string line;

	if (newickFile.is_open()) {
		newickFile >> line;
		parse_newick_tree(line);
	}
	else {
		std::cout << "Tree file does not exist or is not correctly formatted." << std::endl;
		exit (EXIT_FAILURE);
	}

	newickFile.close();

	dfs_iterator = dfs_iterator_recurse(root);
	bfs_iterator = bfs_iterator_recurse(root);
	leave_iterator = leave_iterator_recurse(root);
}

bool Tree::parse_newick_tree(std::string treeStr) {
	std::string nodeName = "";
	std::string distance = "";
	Node *currentNode_pt = root;
	Node *currentDistance_pt = root;

	for (std::string::iterator it = treeStr.begin() + 1; it != treeStr.end(); it++) {
		if (!std::isspace(*it)) {
			switch(*it) {
				case '(': {
					Node* temp = new Node("internal_" + std::to_string(++internalNodeCounter));
					currentNode_pt->add_child(temp);
					currentNode_pt = temp;
					break;
				}
				case ')': {
					currentDistance_pt = currentNode_pt;
					currentNode_pt = currentNode_pt->parent;
					break;
				}
				case ',': {
					break;
				}
				case ':': {
					distance = "";
					it++;
					while ((*it) != ')' and (*it) != '(' and (*it) != ',' and (*it) != ':' and (*it) != ';') {
						distance += (*it);
						it++;
					}
					it--;
					currentDistance_pt->distance = std::stod(distance);
					break;
				}
				case ';': {
					break;
				}
				default: {
					nodeName = "";
					bool createNewNode = true;
					if (*(it-1) == ')') {	// internal node; thus won't create new node and rename existing
						createNewNode = false;
					}
					while ((*it) != ')' and (*it) != '(' and (*it) != ',' and (*it) != ':' and (*it) != ';') {
						nodeName += (*it);
						it++;
					}
					it--;

					if (createNewNode) {
						Node* temp = new Node(nodeName);	// create new leave node and name it
						if (fswm_internal::namesToGenomeIDs.find(nodeName) == fswm_internal::namesToGenomeIDs.end()) {
							std::cerr << "The following sequence name is in the tree, but not in the references: " << nodeName  << std::endl;
							exit (EXIT_FAILURE);
						}
						currentNode_pt->add_child(temp);
						currentDistance_pt = temp;
					}
					else {
						currentDistance_pt->name = nodeName;  // update name of existing internal node
					}
				}
			}
		}
	}
	// Check if tree is unrooted an root arbitrarily at lowest level
	if (root->children.size() > 2) {
		if (fswm_params::g_verbose) { std::cout << "\tThe input tree is unrooted, please use a rooted tree."
			"\n\tThe tree will be rooted at the implicit trifurcating root now." << std::endl; }

		// Save children that will be below new subtree and remove from children
		Node* child2 = root->children[1];
		Node* child3 = root->children[2];
		root->children.erase(root->children.begin()+1,root->children.begin()+3);

		// Create new internal node
		Node* temp = new Node("internal_" + std::to_string(++internalNodeCounter));
		root->add_child(temp);
		currentNode_pt = temp;

		// Add old children of root as children of new internal node
		currentNode_pt->add_child(child2);
		currentNode_pt->add_child(child3);

		// Fix internal node names
		fix_internalNodeNumbers();

		is_rooted = false;
	}
	return true;
}

void Tree::fix_internalNodeNumbers() {
	int i = 1;
	for (auto const& node : dfs_iterator) {
		if (!node->children.empty()) {
			node->name = "internal_" + std::to_string(i++);
		}
	}
}

/** Recursively build vector of node pointers for depth first seach iterators. */
std::vector<Node*> Tree::bfs_iterator_recurse(Node* currentNode) {
	std::vector<Node*> bfs_nodes;

	bfs_nodes.push_back(currentNode);

	for (auto const &child : currentNode->children) {
		for (auto const &node_pt : bfs_iterator_recurse(child)) {
			bfs_nodes.push_back(node_pt);
		}
	}

	return bfs_nodes;
}

std::vector<Node*> Tree::dfs_iterator_recurse(Node* currentNode) {
	std::vector<Node*> dfs_nodes;

	for (auto const &child : currentNode->children) {
		for (auto const &node_pt : dfs_iterator_recurse(child)) {
			dfs_nodes.push_back(node_pt);
		}
	}

	dfs_nodes.push_back(currentNode);

	return dfs_nodes;
}

std::vector<Node*> Tree::leave_iterator_recurse(Node* currentNode) {
	std::vector<Node*> dfs_leaves;

	for (auto const &node : dfs_iterator_recurse(currentNode)) {
		if (node->children.empty()) {
			dfs_leaves.push_back(node);
		}
	}

	return dfs_leaves;
}

/** Return leave with most filtered matching k-mers. */
seq_id_t Tree::get_node_best_count(countMap_t::iterator &countMap_it) {
	seq_id_t bestID = get_rootID();
	count_t maxCount = std::numeric_limits<count_t>::min();
	for (auto const &seqIDtoCount : countMap_it->second) {
		if (seqIDtoCount.second > maxCount) {
			maxCount = seqIDtoCount.second;
			bestID = seqIDtoCount.first;
		}
	}
	return bestID;
}

/** Return leave with highest similarity score. */
seq_id_t Tree::get_node_best_score(scoringMap_t::iterator &scoringMap_it) {
	seq_id_t bestID = get_rootID();
	double bestScore = std::numeric_limits<double>::max();
	for (auto const &seqIDtoScoring : scoringMap_it->second) {
		if (seqIDtoScoring.second < bestScore) {
			bestScore = seqIDtoScoring.second;
			bestID = seqIDtoScoring.first;
		}
	}
	return bestID;
}

/** Return LCA of top n nodes with most filtered matching k-mers. */
seq_id_t Tree::get_LCA_best_count(countMap_t::iterator &it) {
	if (it->second.size() == 0) {
		return get_rootID();
	}
	else if (it->second.size() == 1) {
		return it->second.begin()->first;
	}

	count_t first = std::numeric_limits<count_t>::min();
	seq_id_t first_id = 0;
	count_t second = std::numeric_limits<count_t>::min();
	seq_id_t second_id = 0;

	for (auto const &seqIDtoScoring : it->second) {
		if (seqIDtoScoring.second > first) {
			second = first;
			second_id = first_id;
			first = seqIDtoScoring.second;
			first_id = seqIDtoScoring.first;
		}
		else if (seqIDtoScoring.second > second) {
			second = seqIDtoScoring.second;
			second_id = seqIDtoScoring.first;
		}
	}
	return find_LCA(std::vector<seq_id_t> {first_id, second_id})->ID;
}

/** Return LCA of top n nodes with most filtered matching k-mers. */
seq_id_t Tree::get_LCA_best_count_exp(countMap_t::iterator &it, double div) {
	if (it->second.size() == 0) {
		return get_rootID();
	}
	else if (it->second.size() == 1) {
		return it->second.begin()->first;
	}

	count_t first = std::numeric_limits<count_t>::min();
	seq_id_t first_id = 0;
	count_t second = std::numeric_limits<count_t>::min();
	seq_id_t second_id = 0;

	for (auto const &seqIDtoScoring : it->second) {
		if (seqIDtoScoring.second > first) {
			second = first;
			second_id = first_id;
			first = seqIDtoScoring.second;
			first_id = seqIDtoScoring.first;
		}
		else if (seqIDtoScoring.second > second) {
			second = seqIDtoScoring.second;
			second_id = seqIDtoScoring.first;
		}
	}

	// If highest count is much higher than second highest count, return id of first instead of LCA
	if ((first - second) > (first + second)/div) {
		return first_id;
	}

	return find_LCA(std::vector<seq_id_t> {first_id, second_id})->ID;
}

/** Return LCA of top n nodes with smallest similarity scores. */
seq_id_t Tree::get_LCA_best_score(scoringMap_t::iterator &it) {
	std::vector<std::pair<scoring_t, seq_id_t>> minimal_scores;

	for (auto seqIDtoScoring : it->second) {
		minimal_scores.push_back(std::pair<scoring_t, seq_id_t> {seqIDtoScoring.second, seqIDtoScoring.first});
	}

	std::sort(minimal_scores.begin(), minimal_scores.end());

	if (minimal_scores.size() == 1) {
		return minimal_scores[0].second;
	}
	else if (minimal_scores.size() > 1) {
		return find_LCA(std::vector<seq_id_t> {minimal_scores[0].second, minimal_scores[1].second})->ID;
	}

	return get_rootID();
}

/**
 * Return pointer to node of LCA of all leaves given in parameter vector
 */
Node* Tree::find_LCA(std::vector<seq_id_t> leaves) {
	std::vector<std::vector<Node*>> paths;

	for (int i = 0; i < leaves.size(); i++) {		// Find all paths from leaves to root...
		std::vector<Node*> nodes;
		Node* currentNode = find_node(leaves[i]);
		nodes.push_back(currentNode);

		while (currentNode != root) {
			nodes.push_back(currentNode->parent);
			currentNode = currentNode->parent;
		}

		std::reverse(std::begin(nodes), std::end(nodes));
		paths.push_back(nodes);                  	// ...and save them in vector.
	}

	for (int i = 0; i < paths[0].size(); i++) {     // Go through all paths in parallel...
		Node* currentNode = paths[0][i];
		for (auto const &path : paths) {
			if (path[i] != currentNode) {
				return path[i-1];					//...and return last node that is the same in all paths.
			}
		}
	}

	return paths[0][paths[0].size()-1];
}

/**
 * Return pointer to node of LCA of all leaves given in set
 */
Node* Tree::find_LCA(std::unordered_set<seq_id_t> leaves) {
	std::vector<std::vector<Node*>> paths;

	for (auto const &leave : leaves) {		// Find all paths from leaves to root...
		std::vector<Node*> nodes;
		Node* currentNode = find_node(leave);
		nodes.push_back(currentNode);

		while (currentNode != root) {
			nodes.push_back(currentNode->parent);
			currentNode = currentNode->parent;
		}

		std::reverse(std::begin(nodes), std::end(nodes));
		paths.push_back(nodes);                  	// ...and save them in vector.
	}

	for (int i = 0; i < paths[0].size(); i++) {     // Go through all paths in parallel...
		Node* currentNode = paths[0][i];
		for (auto const &path : paths) {
			if (path[i] != currentNode) {
				return path[i-1];					//...and return last node that is the same in all paths.
			}
		}
	}

	return paths[0][paths[0].size()-1];
}

/** Return pointer to node with ID */
Node* Tree::find_node(seq_id_t seqID) {
	for (auto node : dfs_iterator) {
		if (node->ID == seqID) {
			return node;
		}
	}
	std::cerr << "Could not find node with the given node ID in tree: " << seqID  << std::endl;
	exit (EXIT_FAILURE);
	return nullptr;
}

/** Reset leave similarity scores to 0. */
void Tree::reset_similarityScores() {
	for (auto node : dfs_iterator) {
		node->similarityScore = -1;
	}
}

/**
 * Set leave similarity scores based on scoring map of genomes for one read as created by Scoring.
 */
void Tree::set_similarityScores(scoringMap_t::iterator &it) {
	for (auto const &leave : leave_iterator) {
		if (it->second.find(leave->ID) != it->second.end()) {
			leave->similarityScore = it->second[leave->ID];
		}
		else {
			leave->similarityScore = 10;
		}
	}
}

/** Reset node weights to -1. */
void Tree::reset_weights() {
	for (auto &node : dfs_iterator) {
		node->weight = -1;
	}
}

/**
 * Set leave weights based on scoring map of genomes for one read as created by Scoring.
 */
void Tree::set_weights_to_counts(countMap_t::iterator &it) {
	for (auto const &leave : leave_iterator) {
		if (it->second.find(leave->ID) != it->second.end()) {
			leave->weight = it->second[leave->ID];
		}
		else {
			leave->weight = 0;
		}
	}
}

/**
 * Fill all internal node similarity scores with lowest one from children.
 */
void Tree::fill_internals_min_score() {
	for (auto const &node : dfs_iterator) {
		if (node->similarityScore < 0) {
			double min_similarity_score = std::numeric_limits<double>::max();
			for (auto const child : node->children) {
				if (child->similarityScore < min_similarity_score) {
					min_similarity_score = child->similarityScore;
				}
			}
			node->similarityScore = min_similarity_score;
		}
	}
}

/**
 * Fill all internal node similarity scores with sum of all children.
 */
void Tree::fill_internals_sum_score() {
	for (auto const &node : dfs_iterator) {
		if (node->similarityScore < 0) {
			double sum_similarity_score = 0;
			for (auto const child : node->children) {
				sum_similarity_score += child->similarityScore;
			}
			node->similarityScore = sum_similarity_score;
		}
	}
}

/**
 * Fill all internal node similarity scores with avg of all children.
 */
void Tree::fill_internals_avg_score() {
	for (auto const &node : dfs_iterator) {
		if (node->similarityScore < 0) {
			double avg_similarity_score = 0;
			int weight = 0;
			for (auto const child : node->children) {
				weight += child->weight;
			}
			for (auto const child : node->children) {
				avg_similarity_score += child->similarityScore * child->weight / weight;
			}
			node->similarityScore = avg_similarity_score;
		}
	}
}

/** Fill all internal node similarity scores according to spaced word counts */
void Tree::fill_internals_max_count() {
	for (auto const &node : dfs_iterator) {
		if (node->weight < 0) {
			int max_weight = std::numeric_limits<int>::min();
			for (auto const child : node->children) {
				if (child->weight > max_weight) {
					max_weight = child->weight;
				}
			}
			node->weight = max_weight;
		}
	}
}

/** Fill all internal node weights with sum of all children. */
void Tree::fill_internals_sum_count() {
	for (auto const &node : dfs_iterator_recurse(root)) {
		if (node->children.size() > 0) {
			count_t sum_weight = 0;
			for (auto const child : node->children) {
				sum_weight += child->weight;
			}
			node->weight = sum_weight;
		}
	}
}

/** Fill nodes below fields in all inner nodes. */
void Tree::fill_leaves_below() {
	for (auto const node : dfs_iterator) {
		if (node->children.empty()) {
			node->leaves_below = 1;
		}
		else {
			count_t leaves_below = 0;
			for (auto const child : node->children) {
				leaves_below += child->leaves_below;
			}
			node->leaves_below = leaves_below;
		}
	}
}

/** Check if &child is a node or same as &parent and if so return true, otherwise return false. */
bool Tree::is_child_of(seq_id_t child_id, seq_id_t parent_id) {
	Node* child = find_node(child_id);
	Node* parent = find_node(parent_id);
	if (child == parent) {
		return true;
	}
	if (child == root) {
		return false;
	}
	if (parent == root) {
		return true;
	}

	while (child->parent != root) {
		if (child->parent == parent) {
			return true;
		}
		child = child->parent;
	}
	return false;
}

/**  Write tree to file 'filename' in newick format. */
bool Tree::write_newick(std::string filename) {
	std::ofstream* outputTreeStream = new std::ofstream(filename);
	*outputTreeStream << get_newick_str(false);
	outputTreeStream->close();

	return true;
}

/** Return string of tree in newick format. */
std::string Tree::get_newick_str(bool write_edge_nums = true) {
	std::stringstream outputTreeStream;
	//outputTreeStream << "(";
	get_newick_str_recurse(outputTreeStream, root, 0, write_edge_nums);
	outputTreeStream << ";";

	return outputTreeStream.str();;
}

/** Helper function for get_newick_str(). */
int Tree::get_newick_str_recurse(std::stringstream &outputTreeStream, Node *node, int count, bool write_edge_nums = true) {
	if (!node->children.empty()) {
		outputTreeStream << "(";
		for (auto const child : node->children) {
			count = get_newick_str_recurse(outputTreeStream, child, count, write_edge_nums);
			if (!(child == *(node->children.end() - 1))) {
				outputTreeStream << ",";
			}
		}
		outputTreeStream << ")" << node->name << ":" << node->distance;
		if (write_edge_nums) {
			outputTreeStream << "{" << count << "}";
		}
		fswm_internal::IDsToPlacementIDs[node->ID] = count;
		fswm_internal::placementIDsToIDs[count] = node->ID;
		count++;
	}
	else {
		outputTreeStream << node->name << ":" << node->distance;
		if (write_edge_nums) {
			outputTreeStream << "{" << count << "}";
		}
		fswm_internal::IDsToPlacementIDs[node->ID] = count;
		fswm_internal::placementIDsToIDs[count] = node->ID;
		count++;
	}
	return count;
}

/** Return ID of root. */
seq_id_t Tree::get_rootID() {
	return root->ID;
}

/** Write metainformation of jplace file, such as version, fields, metadata, tree. */
void Tree::write_jplace_data_beginning() {
	std::ofstream jPlaceFile;
	jPlaceFile.open(fswm_params::g_outfoldername + fswm_params::g_outjplacename);
	jPlaceFile << "{\n\t\"version\":3,\n\t"
		"\"fields\":[\"edge_num\",\"distal_length\",\"pendant_length\",\"like_weight_ratio\",\"likelihood\"],\n"
		"\t\"metadata\":{\n"
			"\t\t\"software\"\t:\t\"App-SpaM\",\n"
			"\t\t\"More info\"\t:\t\"https://github.com/matthiasblanke/APP-SpaM\",\n\n"
			"\t\t\"reference_fasta\"\t:\t\"" + fswm_params::g_genomesfname + "\",\n"
			"\t\t\"tree_newick\"\t:\t\"" + fswm_params::g_reftreefname + "\",\n"
			"\t\t\"query_fasta\"\t:\t\"" + fswm_params::g_readsfname + "\",\n"
			"\t\t\"number of patterns\"\t:\t" + std::to_string(fswm_params::g_numPatterns) + ",\n"
			"\t\t\"weight\"\t:\t" + std::to_string(fswm_params::g_weight) + ",\n"
			"\t\t\"dont cares\"\t:\t" + std::to_string(fswm_params::g_spaces) + ",\n"
			"\t\t\"mode\"\t:\t\"" + fswm_params::g_assignmentMode + "\",\n"
			"\t\t\"filtering threshold\"\t:\t" + std::to_string(fswm_params::g_filteringThreshold) + ",\n"
			"\t\t\"sampling\"\t:\t" + std::to_string(fswm_params::g_sampling) + ",\n"
			"\t\t\"minHashLowerLimit\"\t:\t" + std::to_string(fswm_params::g_minHashLowerLimit) + ",\n"
			"\t\t\"unassembled\"\t:\t" + std::to_string(fswm_params::g_draftGenomes) + ",\n"
			"\t\t\"delimiter\"\t:\t\"" + fswm_params::g_delimiter + "\"\n"
			"\t},\n\t"
		"\"tree\":\"" + get_newick_str() + "\",\n"
		"\t\"placements\":\n"
		"\t[\n";

	jPlaceFile.close();
}

/** Write closing brackets after placement data to jplace file. */
void Tree::write_jplace_data_end() {
	std::ofstream jPlaceFile;
	jPlaceFile.open(fswm_params::g_outfoldername + fswm_params::g_outjplacename, std::ios_base::app);
	jPlaceFile << "\t]\n"
		"}";
	jPlaceFile.close();
}

/** For each assigned read, write placement data to jplace file. */
void Tree::write_jplace_placement_data(std::vector<std::pair<seq_id_t, int>> &readAssignment, scoringMap_t &scoringMap) {
	std::ofstream jPlaceFile;
	jPlaceFile.open(fswm_params::g_outfoldername + fswm_params::g_outjplacename, std::ios_base::app);

	int i = 0;

	double distal_length = 0;
	double pendant_length = fswm_params::default_distance_new_leaves;
	double dist_refs = 0;
	double dist_current_edge = 0;

	for (auto const& read : readAssignment) {
		if (fswm_internal::jplace_tracking) {
			;
		}
		else {
			jPlaceFile << ",";
		}
		
		// Determine distal and pendant branch lengths
		if (fswm_params::g_assignmentMode == "MINDIST" or fswm_params::g_assignmentMode == "SPAMCOUNT") {
			dist_refs = scoringMap[read.first][read.second];
			dist_current_edge = find_node(read.second)->distance;
			if (dist_refs < 2*dist_current_edge) {
				distal_length = dist_refs / 2;
				pendant_length = dist_refs / 2;
			}
			else {
				distal_length = dist_current_edge;
				pendant_length = dist_refs - dist_current_edge;
			}
		}
		if (fswm_params::g_assignmentMode == "LCACOUNT" or fswm_params::g_assignmentMode == "LCADIST") {
			distal_length = find_node(read.second)->distance/2; 
		}

		jPlaceFile << "\t\t{\n"
					  "\t\t\t\"p\":\n"
		              "\t\t\t[[" + std::to_string(fswm_internal::IDsToPlacementIDs[read.second]) + "," + 
		              	std::to_string(distal_length) + "," + std::to_string(pendant_length) + ",1,1]],\n"
		              "\t\t\t\"nm\":\n"
		              "\t\t\t[[\"" + fswm_internal::readIDsToNames[read.first] + "\", 1]]\n"
		              "\t\t}";
		i++;
		jPlaceFile << "\n";
		fswm_internal::jplace_tracking = false;
	}

	jPlaceFile.close();
}
