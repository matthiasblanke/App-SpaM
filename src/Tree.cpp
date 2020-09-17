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
						currentNode_pt->add_child(temp);
						currentDistance_pt = temp;
					}
					else {
						;// currentDistance_pt->name = nodeName;  // update name of existing internal node
					}
				}
			}
		}
	}
	// Check if tree is unrooted an root arbitrarily at lowest level
	if (root->children.size() > 2) {
		std::cout << "\t### Be aware that the tree is unrooted and will be rooted arbitrarily. ###" << std::endl;

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
 * Descent tree from root to leaves. Return node, when difference in scores
 * between both subtrees is small.
 */
seq_id_t Tree::get_node_root_descent_best_count(countMap_t::iterator &it) {
	reset_weights();
	set_weights_to_counts(it);
	fill_internals_max_count();

	int rel_diff = 10;
	Node *currentNode_pt = root;
	while (currentNode_pt->children.size() == 2) {
		int score1 = currentNode_pt->children[0]->weight;
		int score2 = currentNode_pt->children[1]->weight;
		if (abs(score1 - score2) < ((score1+score2)/rel_diff)) {
			return currentNode_pt->ID;
		}
		else {
			if (score1 < score2) {
				currentNode_pt = currentNode_pt->children[1];
			}
			else {
				currentNode_pt = currentNode_pt->children[0];
			}
		}
	}
	return currentNode_pt->ID;
}

/**
 * Descent tree from root to leaves. Return node, when difference in scores
 * between both subtrees is small.
 */
seq_id_t Tree::get_node_root_descent_best_score(scoringMap_t::iterator &it) {
	reset_similarityScores();
	set_similarityScores(it);
	fill_internals_min_score();

	int rel_diff = 10;
	Node *currentNode_pt = root;
	while (currentNode_pt->children.size() == 2) {
		double score1 = currentNode_pt->children[0]->similarityScore;
		double score2 = currentNode_pt->children[1]->similarityScore;
		if (abs(score1 - score2) < ((score1+score2)/rel_diff)) {
			return currentNode_pt->ID;
		}
		else {
			if (score1 < score2) {
				currentNode_pt = currentNode_pt->children[0];
			}
			else {
				currentNode_pt = currentNode_pt->children[1];
			}
		}
	}
	return currentNode_pt->ID;
}

seq_id_t Tree::get_node_best_weighted(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count) {
	reset_similarityScores();
	set_similarityScores(it_score);
	reset_weights();
	set_weights_to_counts(it_count);
	fill_internals_weighted_min_score();

	Node *currentNode_pt = root;
	while (currentNode_pt->children.size() == 2) {
		double score1 = currentNode_pt->children[0]->similarityScore;
		double score2 = currentNode_pt->children[1]->similarityScore;
		if (score1 < score2) {
			currentNode_pt = currentNode_pt->children[0];
		}
		else {
			currentNode_pt = currentNode_pt->children[1];
		}
	}

	return get_rootID();
}

seq_id_t Tree::get_LCA_best_weighted(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count) {
	reset_similarityScores();
	set_similarityScores(it_score);
	reset_weights();
	set_weights_to_counts(it_count);
	fill_internals_weighted_min_score();

	std::vector<std::pair<scoring_t, seq_id_t>> minimal_scores;

	for (auto const &leave : leave_iterator) {
		minimal_scores.push_back(std::pair<scoring_t, seq_id_t> {leave->similarityScore, leave->ID});
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
 * Descent tree from root to leaves. Return node, when difference in scores
 * between both subtrees is small.
 */
seq_id_t Tree::get_node_root_descent_best_weighted_score(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count) {
	reset_similarityScores();
	set_similarityScores(it_score);
	reset_weights();
	set_weights_to_counts(it_count);
	fill_internals_weighted_min_score();

	int rel_diff = 10;
	Node *currentNode_pt = root;
	while (currentNode_pt->children.size() == 2) {
		double score1 = currentNode_pt->children[0]->similarityScore;
		double score2 = currentNode_pt->children[1]->similarityScore;
		if (abs(score1 - score2) < ((score1+score2)/rel_diff)) {
			return currentNode_pt->ID;
		}
		else {
			if (score1 < score2) {
				currentNode_pt = currentNode_pt->children[0];
			}
			else {
				currentNode_pt = currentNode_pt->children[1];
			}
		}
	}
	return currentNode_pt->ID;
}

void Tree::score_all(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count, bool first) {
	std::ofstream weighting_out;
	weighting_out.open(fswm_params::g_outfoldername + "weighting.txt", std::ios_base::app);
	weighting_out << "ReadID\tNodeID\tNodeName\tqual_score\tavg_count\tsum_count\tavg_score\tsum_score\tleaves_below" << std::endl;

	std::vector<std::pair<double, seq_id_t>> quality_values;

	fill_leaves_below();

	reset_weights();
	set_weights_to_counts(it_count);
	fill_internals_sum_count();

	reset_similarityScores();
	set_similarityScores(it_score);
	fill_internals_sum_score();

	// Fill in node information
	int i = 0;
	for (auto const &node : dfs_iterator) {

		double avg_count = (double) node->weight / node->leaves_below;
		double avg_score = abs(node->similarityScore / node->leaves_below);
		double quality_score = (custom_weighting_function_count(node->leaves_below)) * (avg_count + 0.001);
		quality_score = quality_score * (1 / (abs(avg_score) + 0.001));

		quality_values.push_back(std::pair<double, seq_id_t>(quality_score, node->ID));
		weighting_out << it_score->first << "\t" << node->ID << "\t" << node->name << "\t" << quality_score << "\t"
					  << avg_count << "\t" << node->weight << "\t" << avg_score << "\t" << node->similarityScore << "\t" << node->leaves_below << std::endl;
	}
	std::cout << std::endl;

	std::sort(quality_values.rbegin(), quality_values.rend());

	double sum = 0;
	for (int i = 0; i < 2; i++) {
		sum += quality_values[i].first;
	}

	std::vector<std::pair<seq_id_t, double>> placements;
	for (int i = 0; i < 2; i++) {
		placements.push_back(std::pair<seq_id_t, double> (quality_values[i].second, quality_values[i].first / sum));
	}

	for (auto const &entry : quality_values) {
		std::cout << entry.first << "\t" << entry.second << std::endl;
	}

	/*
	double q1 = quality_values[int(1 * quality_values.size() / 10)].first;
	double q3 = quality_values[int(3 * quality_values.size() / 10)].first;
	double iqr = q1 - q3;
	double threshold = quality_values[0].first / 2;
	std::cout << "Q1: " << q1 << "\tQ3: " << q3 << "\tiqr: " << iqr << "\tt: " << threshold << std::endl;

	double sum = 0;
	for (auto const &qual_val : quality_values) {
		if (qual_val.first > threshold) {
			sum += qual_val.first;
		}
	}

	std::vector<std::pair<seq_id_t, double>> placements;
	for (auto const &qual_val : quality_values) {
		if (qual_val.first > threshold) {
			placements.push_back(std::pair<seq_id_t, double> (qual_val.second, qual_val.first / sum));
		}
	}

	if (placements.size() == 0) {
		placements.push_back(std::pair<seq_id_t, double> (quality_values[0].second, 1));
	}*/

	write_multiple_jplace(placements, first, it_score->first);
}

void Tree::score_all_avg_sim_score(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count, bool first) {
	std::vector<std::pair<double, seq_id_t>> quality_values;

	fill_leaves_below();
	reset_similarityScores();
	set_similarityScores(it_score);
	fill_internals_sum_score();

	// Fill in node information
	int i = 0;
	for (auto const &node : dfs_iterator) {

		double average_score = node->similarityScore / node->leaves_below;
		double quality_score = (custom_weighting_function_score(node->leaves_below)) * (average_score + 0.01);

		quality_values.push_back(std::pair<double, seq_id_t>(quality_score, node->ID));
	}

	std::sort(quality_values.begin(), quality_values.end());

	std::vector<std::pair<seq_id_t, double>> placements;

	double sum = 0;
	for (int i = 0; i < 3; i++) {
		sum += quality_values[i].first;
	}

	for (int i = 0; i < 3; i++) {
		placements.push_back(std::pair<seq_id_t, double> (quality_values[i].second, quality_values[i].first / sum));
	}

	write_multiple_jplace(placements, first, it_score->first);
}

/** Write metainformation of jplace file, such as version, fields, metadata, tree. */
void Tree::write_multiple_jplace(std::vector<std::pair<seq_id_t, double>> placements, bool first, seq_id_t seqID) {
	std::ofstream jPlaceFile;
	jPlaceFile.open(fswm_params::g_outfoldername + fswm_params::g_outjplacename, std::ios_base::app);

	if (!first) {
		jPlaceFile << "\t\t,\n";
	}
	jPlaceFile << "\t\t{\n"
				  "\t\t\t\"p\":\n\t\t\t[\n";

	first = true;
	for (auto const &placement : placements) {
		if (!first) {
			jPlaceFile << "\t\t\t,\n";
		}
		jPlaceFile << "\t\t\t [" + std::to_string(fswm_internal::IDsToPlacementIDs[placement.first]) + "," + std::to_string(find_node(placement.first)->distance/2) + "," +
		std::to_string(fswm_params::default_distance_new_leaves) + "," << placement.second << ",1]\n";
		first = false;
	}
	if (placements.size() == 0) {
		jPlaceFile << "\t\t\t [" + std::to_string(fswm_internal::IDsToPlacementIDs[root->ID]) + "," + std::to_string(find_node(root->ID)->distance/2) + "," +
		std::to_string(fswm_params::default_distance_new_leaves) + ",1,1]\n";
		first = false;
	}
	jPlaceFile << "\n\t\t\t],\n";
	jPlaceFile << "\t\t\t\"nm\":\n"
		              "\t\t\t[[\"" + fswm_internal::readIDsToNames[seqID] + "\", 1]]\n"
		              "\t\t}\n";
	jPlaceFile << "\n";

	jPlaceFile.close();
}

void Tree::score_all_avg_sim_count(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count, bool first) {
	std::ofstream weighting_out;
	weighting_out.open(fswm_params::g_outfoldername + "weighting.txt", std::ios_base::app);
	weighting_out << "ReadID\tNodeID\tNodeName\tqual_score\tavg_score\tcount\tleaves_below" << std::endl;

	std::vector<std::pair<double, seq_id_t>> quality_values;

	fill_leaves_below();
	reset_weights();
	set_weights_to_counts(it_count);
	fill_internals_sum_count();

	// Fill in node information
	int i = 0;
	for (auto const &node : dfs_iterator) {

		double average_score = (double) node->weight / node->leaves_below;
		double quality_score = (custom_weighting_function_count(node->leaves_below)) * (average_score + 0.01);

		quality_values.push_back(std::pair<double, seq_id_t>(quality_score, node->ID));
		weighting_out << it_score->first << "\t" << node->ID << "\t" << node->name << "\t" << quality_score << "\t"
					  << average_score << "\t" << node->weight << "\t" << node->leaves_below << std::endl;
 		// std::cout << "AVG: " << average_score << "\tGAMMA: " << quality_score << "\tNB: " << node->leaves_below << std::endl;
	}
	std::cout << std::endl;

	std::sort(quality_values.rbegin(), quality_values.rend());

	for (auto const &entry : quality_values) {
		std::cout << entry.first << "\t" << entry.second << std::endl;
	}
	std::cout << std::endl;

	/*double q1 = quality_values[int(1 * quality_values.size() / 4)].first;
	double q3 = quality_values[int(3 * quality_values.size() / 4)].first;
	double iqr = q1 - q3;
	double threshold = q1 + 7 * iqr;
	std::cout << "Q1: " << q1 << "\tQ3: " << q3 << "\tiqr: " << iqr << "\tt: " << threshold << std::endl;

	double sum = 0;
	for (auto const &qual_val : quality_values) {
		if (qual_val.first > threshold) {
			sum += qual_val.first;
		}
	}

	std::vector<std::pair<seq_id_t, double>> placements;
	for (auto const &qual_val : quality_values) {
		if (qual_val.first > threshold) {
			placements.push_back(std::pair<seq_id_t, double> (qual_val.second, qual_val.first / sum));
		}
	}

	if (placements.size() == 0) {
		placements.push_back(std::pair<seq_id_t, double> (quality_values[0].second, 1));
	}
	*/

	double sum = 0;
	for (int i = 0; i < 4; i++) {
		sum += quality_values[i].first;
	}

	std::vector<std::pair<seq_id_t, double>> placements;
	for (int i = 0; i < 4; i++) {
		placements.push_back(std::pair<seq_id_t, double> (quality_values[i].second, quality_values[i].first / sum));
	}

	for (auto const &entry : placements) {
		std::cout << entry.first << "\t" << entry.second << std::endl;
	}
	std::cout << std::endl;
	std::cout << "-------------" << std::endl;

	std::ofstream jPlaceFile;
	jPlaceFile.open(fswm_params::g_outfoldername + fswm_params::g_outjplacename, std::ios_base::app);

	if (!first) {
		jPlaceFile << "\t\t,\n";
	}
	jPlaceFile << "\t\t{\n"
				  "\t\t\t\"p\":\n\t\t\t[\n";

	first = true;
	for (auto const &placement : placements) {
		if (!first) {
			jPlaceFile << "\t\t\t,\n";
		}
		jPlaceFile << "\t\t\t [" + std::to_string(fswm_internal::IDsToPlacementIDs[placement.first]) + "," + std::to_string(find_node(placement.first)->distance/2) + "," +
		std::to_string(fswm_params::default_distance_new_leaves) + "," << placement.second << ",1]\n";
		first = false;
	}
	jPlaceFile << "\n\t\t\t],\n";
	jPlaceFile << "\t\t\t\"nm\":\n"
		              "\t\t\t[[\"" + fswm_internal::readIDsToNames[it_score->first] + "\", 1]]\n"
		              "\t\t}\n";
	jPlaceFile << "\n";

	jPlaceFile.close();

	weighting_out.close();
}

double Tree::gamma_k2(double x) {
	return 1 - exp(-x / 2);
}

double Tree::custom_weighting_function_score(int x) {
	switch(x) {
		case 1: return 1; break;
		case 2: return 0.97; break;
		case 3: return 0.93; break;
		case 4: return 0.9; break;
		case 5: return 0.87; break;
		case 6: return 0.83; break;
		case 7: return 0.8; break;
		default: return 0.8; break;
	}
}

double Tree::custom_weighting_function_count(int x) {
	switch(x) {
		case 1: return 0.9; break;
		case 2: return 0.95; break;
		case 3: return 1; break;
		case 4: return 1.05; break;
		case 5: return 1.1; break;
		case 6: return 1.15; break;
		case 7: return 1.2; break;
		case 8: return 1.25; break;
		default: return 1.3; break;
	}
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
 * Set leave weights based on scoring map of genomes for one read as created by Scoring.
 */
void Tree::set_weights_to_lca_counts(std::unordered_map<seq_id_t,int> &map) {
	for (auto const &node : dfs_iterator) {
		if (map.find(node->ID) != map.end()) {
			node->weight = map[node->ID];
		}
		else {
			node->weight = 0;
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

/**
 * Fill all internal node similarity scores
 */
void Tree::fill_internals_weighted_min_score() {
	for (auto const &node : dfs_iterator) {
		if (node->children.empty()) {
			node->similarityScore = node->similarityScore * 1 / node->weight;
		}
		else if (node->similarityScore < 0) {
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

/** Test if given number of seqIDs is a monophyletic group in tree with an allowance of n nodes. */
bool Tree::test_monophyly(std::unordered_set<seq_id_t> &seqID_set) {
	Node* lca = find_LCA(seqID_set);	// Find LCA of seqIDs

	int count_strange_leaves = 0; 		// number of leaves that are not in seqID_set

	for (auto const& node : leave_iterator_recurse(lca)) {		// If some node below LCA...
		if (seqID_set.find(node->ID) == seqID_set.end()) { 	// ...is not in set of seqIDs...
			count_strange_leaves++;
			if (count_strange_leaves > fswm_params::g_allowance) {
				return false;		// ...the group is not monophyletic.
			}
		}
	}

	return true;	// All nodes below LCA are in group, i.e. the group is monophyletic.
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
	jPlaceFile.open(fswm_params::g_outfoldername + fswm_params::g_outjplacename, std::ios_base::app);
	jPlaceFile << "{\n\t\"version\":3,\n\t"
		"\"fields\":[\"edge_num\",\"distal_length\",\"pendant_length\",\"like_weight_ratio\",\"likelihood\"],\n"
		"\t\"metadata\":{\"invocation\":\"\"},\n\t"
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
void Tree::write_jplace_placement_data(std::vector<std::pair<seq_id_t, int>> readAssignment) {
	std::ofstream jPlaceFile;
	jPlaceFile.open(fswm_params::g_outfoldername + fswm_params::g_outjplacename, std::ios_base::app);

	int i = 0;
	bool first = true;

	for (auto const& read : readAssignment) {
		if (first) {
			;
		}
		else {
			jPlaceFile << ",";
		}
		jPlaceFile << "\t\t{\n"
					  "\t\t\t\"p\":\n"
		              "\t\t\t[[" + std::to_string(fswm_internal::IDsToPlacementIDs[read.second]) + "," + std::to_string(find_node(read.second)->distance/2) + "," + std::to_string(fswm_params::default_distance_new_leaves) + ",1,1]],\n"
		              "\t\t\t\"nm\":\n"
		              "\t\t\t[[\"" + fswm_internal::readIDsToNames[read.first] + "\", 1]]\n"
		              "\t\t}";
		i++;
		jPlaceFile << "\n";
		first = false;
	}

	jPlaceFile.close();
}
