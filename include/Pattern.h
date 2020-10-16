/**
 * This programm calculates the variance of a set of pattern with the same length and weight.
 * It is possible to improve your patternset, estimate values for p, q and S, l_hom, li, lj
 *   from a multiple alignment file in fasta format, and also read patterns from a file.
 *
 * pattern object header
 *
 * For theory please have a look at:
 *
 * B. Morgenstern, B. Zhu, S. Horwege, C.-A Leimeister (2015)
 * Estimating evolutionary distances between genomic sequences from spaced-word matches
 * Algorithms for Molecular Biology 10, 5. (http://www.almob.org/content/10/1/5/abstract)
 *
 *
 * @author: Lars Hahn - 03.06.2015, Georg-August-Universitaet Goettingen
 * @version: 1.0 06/2015
 */
#ifndef FSWM_PATTERN_H_
#define FSWM_PATTERN_H_


#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <string.h>

class Pattern{
	public:
		Pattern();
		Pattern(long size, long length, int weight, uint seed);
		~Pattern();

		void ReinitPattern();

		std::vector<std::string> GetPattern();
		std::vector<std::string> GetBestPattern();
		std::string GetPattern(uint number);
		std::string GetBestPattern(uint number);

		double Variance();
		double GetVariance();
		double GetBestVariance();
		double GetNormVariance();
		double GetBestNormVariance();
		int GetWeight();
		int GetSize();
		int GetLength();
		int GetWorstPatMaxVal();
		int GetWorstPatMaxPat();
		bool UniqPattern(uint number);
		void Improve(int limit);
		void ImproveLoop(int limit);
		void ImproveMaxValue(int limit);
		void ImproveMaxValuePattern(int limit);
		void ImproveSecure();

		void Quiet();
		void Silent();
		void Print();
		void ChangePatternRandom(uint number);

	protected:
		void InitMatrix();

		std::vector<std::string> CreateRandomPattern();
		std::vector<std::string> PatternCopy(std::vector<std::string>old_pattern);

		double CalcVariance();
		int ShiftPos(int p1, int p2, int s);
		int WorstPattern_max_val();
		int WorstPattern_max_pat();
		void DoImprove(int limit, bool max_val, bool max_pat, bool loop);

		int LengthSeq(std::vector<std::string> seq);

		double Gauss();
		double MaxNumberPattern(int p_weight, int p_length);
		double Faculty(int value);
		void SecureMessage(std::string errmsg, int pos);


	private:
		std::vector<std::vector<double> > q_values;
		std::vector<std::vector<std::string> > seq_matrix;
		std::vector<std::string> pattern_set;
		std::vector<std::string> best_pattern;
		std::vector<std::vector<double> > var_sum;			/*Contains for each pattern pair the share of the complete variance...*/
		std::vector<int> seq_leng;
		double variance;
		double best_variance;
		long size;
		long length;
		int weight;
		double p;
		double q;
		bool improve;
		bool quiet;
		bool silent;
		bool secure;
};
#endif
