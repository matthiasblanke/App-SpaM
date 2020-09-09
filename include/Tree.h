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
#ifndef FSWM_TREE_H_
#define FSWM_TREE_H_

#include <unordered_set>
#include <sstream>
#include "Node.h"
#include "Scoring.h"
#include "BucketManager.h"
#include "Algorithms.h"

class Tree {
	private:
		Node *root;
		int internalNodeCounter;
		bool is_rooted;

		bool parse_newick_tree(std::string treeStr);
		std::vector<Node*> bfs_iterator_recurse(Node* currentNode);
		std::vector<Node*> dfs_iterator_recurse(Node* currentNode);
		std::vector<Node*> leave_iterator_recurse(Node* currentNode);
		int get_newick_str_recurse(std::stringstream &outputTreeStream, Node *node, int count, bool write_edge_nums);

	public:
		Tree(std::string filename);

		bool write_newick(std::string filename);
		std::vector<Node*> dfs_iterator;
		std::vector<Node*> leave_iterator;
		std::vector<Node*> bfs_iterator;

		// Assignment mode methods
		seq_id_t get_node_best_count(countMap_t::iterator &countMap_it);
		seq_id_t get_node_best_score(scoringMap_t::iterator &scoringMap_it);
		seq_id_t get_node_best_weighted(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count);

		seq_id_t get_LCA_best_count(countMap_t::iterator &it);
		seq_id_t get_LCA_best_score(scoringMap_t::iterator &it);
		seq_id_t get_LCA_best_weighted(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count);

		seq_id_t get_node_root_descent_best_score(scoringMap_t::iterator &it);
		seq_id_t get_node_root_descent_best_count(countMap_t::iterator &it);
		seq_id_t get_node_root_descent_best_weighted_score(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count);
		seq_id_t get_node_root_descent_best_path();

		void score_all(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count, bool first);
		void score_all_avg_sim_score(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count, bool first);
		void score_all_avg_sim_count(scoringMap_t::iterator &it_score, countMap_t::iterator &it_count, bool first);

		double gamma_k2(double x);
		double custom_weighting_function_score(int x);
		double custom_weighting_function_count(int x);

		// Helper functions for assignment mode methods
		void fill_internals_min_score();
		void fill_internals_sum_score();
		void fill_internals_avg_score();
		void fill_internals_max_count();
		void fill_internals_sum_count();
		void fill_internals_sum_count_phylokmers();
		void fill_internals_weighted_min_score();

		void fill_leaves_below();

		// Fill node weights and scores
		void set_similarityScores(scoringMap_t::iterator &it);
		void reset_similarityScores();
		void set_weights_to_counts(countMap_t::iterator &it);
		void set_weights_to_lca_counts(std::unordered_map<seq_id_t,int> &map);
		void reset_weights();

		seq_id_t get_rootID();
		Node* find_LCA(std::vector<seq_id_t> leaves);
		Node* find_LCA(std::unordered_set<seq_id_t> leaves);
		Node* find_node(seq_id_t seqID);
		bool is_child_of(seq_id_t child_id, seq_id_t parent_id);

		// Phylo kmer methods
		void build_phylokmer_db(BucketManager &bucketManager);
		void build_kmercount_db(BucketManager &bucketManager);
		bool test_monophyly(std::unordered_set<seq_id_t> &seqID_set);

		// JPlace writing
		void write_jplace_data_beginning();
		void write_jplace_data_end();
		void write_jplace_placement_data(std::vector<std::pair<seq_id_t, int>> readAssignment);
		void write_multiple_jplace(std::vector<std::pair<seq_id_t, double>> placements, bool first, seq_id_t seqID);
		std::string get_newick_str(bool write_edge_nums);
		void write_phyloDB_toFile();

		void fix_internalNodeNumbers();
};

#endif
