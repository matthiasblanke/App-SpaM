/**
 * This programm calculates the variance of a set of pattern with the same length and weight.
 *
 * pattern object file
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
#include "Pattern.h"


/*---Variables---------------------------------------------------------------*/

std::default_random_engine generator(std::random_device{}());

/*---Constructor-& Init------------------------------------------------------*/
/**
 * Default constructor, sets the default vaulues, pattern will be generated automatically.
 */
Pattern::Pattern() {
	this->size = 10;
	this->length = 14;
	this->weight = 8;
	Pattern(size, length, weight, 0);
}

/**
 * Long constructor, sets the values; resets automatically, if there are problems.
 *
 * @param size
 * 		The amount of patterns; pattern number.
 *
 * @param length
 *		The pattern length for each pattern of the pattern set.
 *
 * @param weigth
 * 		The weight (match positions; '1') for each pattern of the pattern set.
 */
Pattern::Pattern(long size, long length, int weight, uint seed) {
	this->size = size;
	this->length = length;
	this->weight = weight;
	this->p = 0.75;
	this->q = 0.25;
	this->variance = 0;
	this->best_variance = 0;
	this->quiet = false;
	this->silent = false;
	this->secure = false;
	generator.seed(seed);
	ReinitPattern();
}

/**
 * Default destructor, deletes all vectors and matrices in the object.
 */
Pattern::~Pattern() {
	for (uint i = 0; i < seq_matrix.size(); i++) {
		seq_matrix[i].clear();
	}
	seq_matrix.clear();
	seq_leng.clear();

	for (uint i = 0; i < var_sum.size(); i++) {
		var_sum[i].clear();
	}
	var_sum.clear();

	pattern_set.clear();
	best_pattern.clear();
}

/**
 * Creates for the submitted or default values a set of pattern and calculates the first variance.
 * If possible estimates p, q, all sequences lengths and all combinations of homologous sequence positions
 *
 *Also reset Pattern if needed, not necessary to create new object
 */
void Pattern::ReinitPattern() {
	std::ifstream check_align;
	std::vector<std::string> pattern_tmp;
	std::string tmp;
	char tokens[4] = {'.',' ',',',';'};			/*These tokens are allowed to seperate patterns*/

	pattern_set.clear();
	best_pattern.clear();

	if (size <= 0) {						/*Coping with all ...*/
		SecureMessage("patsize", -1);
		this->size = 10;
	}

	if (weight <= 0 || size <= 0 || length <= 1) {		/*... possible wrong ...*/
		SecureMessage("nkl",-1);
		this->size = 10;
		this->weight = 8;
		this->length = 14;
		this->pattern_set.clear();
		this->best_pattern.clear();
	}

	if (weight > length) {					/*..by the user submitted conditions for weight, length etc.*/
		SecureMessage("weight_pat", -1);
		this->weight = 8;
		this->length = 14;
	}
	if (size >= int(MaxNumberPattern(weight-2, length-2))) {	/*Match positions at the start and end do not alterate, therefore -2*/
		size = int(MaxNumberPattern(weight-2, length-2));
		SecureMessage("max_number_pattern",size);	/*We can create all patterns directly*/
		this->improve = false;
		pattern_set.clear();
	}
	else {
		this->improve = true;

	}
	if (improve && (length < 4 || weight == length || weight < 3)) {
		SecureMessage("noimprove",-1);			/*The improvement switches positions with '1' and '0'*/
		this->improve = false;				/*At least 2 positions + end '1' and start '0' are nedeed*/
	}

	if (pattern_set.size() == 0) {				/*worst case, if not wished, generating autopattern*/
		this->pattern_set = CreateRandomPattern();
		for(uint i = 0; i < size; i++) {
			ChangePatternRandom(i);			/*Creating just uniq Pattern!*/
		}
		this->best_pattern = PatternCopy(pattern_set);
		//std::cout << "\n ... Done!\n" << std::endl;
	}

	InitMatrix();

	this->variance = CalcVariance();
	this->best_variance = variance;
	this->best_pattern = PatternCopy(pattern_set);

	return;
}

/**
 * By if it may used more than once, this method is saved.
 */

void Pattern::InitMatrix() {
	std::vector<double> tmp;
	for(int i = 0; i < size; i++) {
		tmp.push_back(0.0);
	}
	for(int i = 0; i < size; i++) {
		var_sum.push_back(tmp);
	}
}


/*---Get-&-SetFunc-----------------------------------------------------------*/
/**
 * Returns complete pattern set
 * @return complete pattern set
 */
std::vector<std::string> Pattern::GetPattern() {
	return pattern_set;
}

/**
 * Returns complete best pattern set
 * @return complete best pattern set
 */

std::vector<std::string> Pattern::GetBestPattern() {
	return best_pattern;
}

/**
 * Returns a specific pattern of the pattern set
 * @return pattern_set[number], specific string
 */
std::string Pattern::GetPattern(uint number) {
	if(number >=size) {
		SecureMessage("wrongindex", number);
		return NULL;
	}
	else{
		return pattern_set[number];
	}
}

/**
 * Returns a specific pattern of the best pattern set
 * @return best_pattern[number], specific string
 */
std::string Pattern::GetBestPattern(uint number) {
	if(number >=size) {
		SecureMessage("wrongindex", number);
		return NULL;
	}
	else{
		return best_pattern[number];
	}
}

/**
 * Returns the current variance
 * @return returns variance
 */
double Pattern::GetVariance() {
	return variance;
}

/**
 * Returns the current best variance
 * @return returns best variance
 */
double Pattern::GetBestVariance() {
	return best_variance;
}

/**
 * Returns the current variance, normalized
 * @return returns norm_variance
 */
double Pattern::GetNormVariance() {
	return variance/Gauss();
}

/**
 * Returns the current best variance, normalized
 * @return returns  best norm_variance
 */
double Pattern::GetBestNormVariance() {
	return best_variance/Gauss();
}

/**
 * Returns the weight of each Pattern, the match positions
 * @return returns weight
 */
int Pattern::GetWeight() {
	return weight;
}

/**
 * Returns the amount of patters; number of Patterns
 * @return returns size
 */
int Pattern::GetSize() {
	return size;
}

/**
 * Returns the length of each pattern
 *
 * @return returns length
 */
int Pattern::GetLength() {
	return length;
}

/**
 * Returns the position of the worst Pattern, estimated by the maximum variancepart
 * 	for each Pattern pair
 *
 * @return returns position worst matrix by max_value
 */
int Pattern::GetWorstPatMaxVal() {
	return WorstPattern_max_val();
}

/**
 * Returns the position of the worst Pattern, estimated by the summation for each Pattern Pi
 * 	and the corresponding variance parts.
 *
 * @return returns position worst matrix by max_pat
 */
int Pattern::GetWorstPatMaxPat() {
	return WorstPattern_max_pat();
}

/**
 * Creates random a set of pattern. For convention a pattern has to start and end with  '1'
 * 	Reason 10010 ~ 1001
 *
 * @return Returns a randomly created pattern set
 */
std::vector<std::string> Pattern::CreateRandomPattern() {
	std::vector<std::string> pattern;
	std::string prototype(length,'0'), tmp;
	uint position, counter;
	// generator.seed(0);

	prototype[0] = '1';
	prototype[length-1] = '1';

	std::uniform_int_distribution<int> distribution(1,length-2);	/*better random generator then time, uses likelihood for evenly random distribution*/

	for(int i = 0; i < size; i++) {
		counter = 2;
		tmp = prototype;					/*Prototype which saves the start and end '1'*/
		while(counter < weight) {
			position = distribution(generator);
			if(tmp[position] == '0') {			/*have fun and fill with '1' randomly :D */
				counter ++;
				tmp[position] = '1';
			}
		}
		pattern.push_back(tmp);
	}
	return pattern;
}

/**
 * Copy the pattern set
 *
 * @param old_pattern
 * 		The pattern which has to be copied
 *
 * @return A new created copy of the pattern set
 */
std::vector<std::string> Pattern::PatternCopy(std::vector<std::string>old_pattern) {
	std::vector<std::string> new_pattern;
	for(uint i = 0; i < old_pattern.size(); i++) {
		new_pattern.push_back(old_pattern[i]);
	}
	return new_pattern;
}

/**
 * Changes two different positions ('1' and '0') in a specific pattern
 * Start and end are excluded
 *
 * @param number
 * 		The pattern which has to be modified
 */
void Pattern::ChangePatternRandom(uint number) {
	uint pos1, pos2;
	bool flag = true;
	char c;

	std::uniform_int_distribution<int> distribution(1,length-2);

	while(flag && improve) {
		pos1 = distribution(generator);
		pos2 = distribution(generator);
		if(pattern_set[number][pos1] != pattern_set[number][pos2]) {		/*Changes two positions, switching '1' and 1' or '0' and '0'...*/
			flag = false;							/*... does not make really sense ;) */
			c = pattern_set[number][pos1];
			pattern_set[number][pos1] = pattern_set[number][pos2];
			pattern_set[number][pos2] = c;
		}
		if(!UniqPattern(number)) {
			flag = true;
		}
	}

	return;
}

/**
 * Scans if there is another pattern in the same format
 *
 * @return returns boolean if there is another same pattern
 */
bool Pattern::UniqPattern(uint number) {
	bool uniq = true;
	for(uint i = 0; i < size; i++) {
		if( number != i) {
			if(pattern_set[i] == pattern_set[number]) {
				uniq = false;
			}
		}
	}
	return uniq;
}

/*---Variance-----------------------------------------------------------------*/
/**
 * Wether there is an alignment or not it decides on its own which variance
 * 	calculation is now neccessary
 *
 * @return returns current variance
 */
double Pattern::Variance() {
	return CalcVariance();
}

/**
 * Calculates the variance for a pattern set
 * In this case l_hom = l1 = l2 = l (commandline parameter S)
 *
 * @return Calculates and returns current variance
 */
double Pattern::CalcVariance() {
	double homologue, background, var_hom, var_bac;
	int shift;

	homologue = 0.0;
	background = 0.0;
	for(int i = 0; i < size; i++) {							/*i and j represents Pi and Pj of the set of pattern*/
		for(int j = i; j < size; j++) {
			var_hom = 0.0;
			var_bac = 0.0;
			for(int s = -1*length+1; s < length; s++) {			/*As in the formula, the shift goes from max shift left to max shift right*/
				shift = ShiftPos(i, j, s);				/*At least one position has to overlap*/
				var_hom += (pow(p, shift) - pow(p, 2*weight));		/*summation of the homologue first part*/
				var_bac += (pow(q, shift) - pow(q, 2*weight));		/*summation of the background second part*/
			}
			var_sum[i][j]=(length - length + 1)*var_hom + (length - length + 1)*(length - length)*var_bac;		/*For each pair Pi and Pj this is the direct share of the complete variance...*/
			var_sum[j][i]=var_sum[i][j];									/*...which we can use to estimate a "worst" pattern*/
											/*... and save it in an size x size matrix*/
			homologue += var_hom;
			background += var_bac;

		}
	}
	this->variance = (length - length + 1)*homologue + (length - length + 1)*(length - length)*background;
	if(this->variance < this->best_variance) {
		this->best_pattern = PatternCopy(pattern_set);
		this->best_variance = variance;
	}
	return variance;

}

/**
 * Shifts to pattern to see the number of maximum match positions
 *
 * @param p1
 *		Position of the first used pattern of the pattern set
 *
 * @param p2
 * 		Position of the second used pattern of the pattern set
 *	NOTE: possible is p1 = p2
 *
 * @param s
 * 		The shift of the second pattern, s < 0 := shift left pattern 2, s > := shift right pattern 2
 * @return Calculates and returns current variance
 */
int Pattern::ShiftPos(int p1, int p2, int s) {
	int counter;

	std::string pat1, pat2;

	counter = 0;

	if(s < 0) {
		s = 0 - s;
		pat1 = pattern_set[p2];				/*If s < 0 for pat2 it is in the point of view for pat1 like s > 0*/
		pat2 = pattern_set[p1];				/*Therefore changing the pattern and take the absolute value from s is easier*/
	}
	else{
		pat1 = pattern_set[p1];
		pat2 = pattern_set[p2];
	}
								/*Afterwards the second pattern (may changed) is only shifted left, its easier*/

	for(int i = 0; i < s; i++) {				/*The pattern section from only pattern 1*/
		if(pat1[i] == '1') {
			counter++;
		}
	}
	for(int i = 0; i < (length-s); i++) {			/*The common section of the pattern*/
		if(pat1[i+s] == '1' || pat2[i] == '1') {
			counter++;
		}
	}
	for(int i = length-s; i < length; i++) {			/*The pattern section from only pattern 2*/
		if(pat2[i] == '1') {
			counter++;
		}
	}

	return counter;
}

/**
 * Returns the position of the worst pattern, estimated by the maximum variancepart
 * 	for each pattern pair
 *
 * @return returns position worst matrix by max_value
 */
int Pattern::WorstPattern_max_val() {
	double max_var_val, i_sum, j_sum;
	int i_max; int j_max;

	max_var_val = 0.0;
	i_sum = 0.0;
	j_sum = 0.0;
	i_max = 0;
	j_max = 0;


	for(int i = 0; i < size; i++) {
		for(int j = 0; j < size; j++) {
			if(var_sum[i][j] > max_var_val) {
				max_var_val = var_sum[i][j];		/*looking for the maximum value and so for the maximum pattern pair*/
				i_max = i;
				j_max = j;
			}
		}
	}


	for(int i = 0; i < size; i++) {
		i_sum += var_sum[i][j_max];				/*If Pi has the higher summation, the part variance of Pi with each corresponding pattern is higher, so more worse*/
	}

	for(int j = 0; j < size; j++) {
		j_sum += var_sum[i_max][j];				/*If Pj has the higher summation, the part variance of PJ with each corresponding pattern is higher, so more worse*/
	}

	if(i_sum > j_sum) {
		return i_max;
	}
	else{
		return j_max;
	}
}

/**
 * Returns the position of the worst pattern, estimated by the maximum summation
 * 	for each pattern with each corresponding variance pattern part
 *
 * @return returns position worst matrix by max_pattern
 */
int Pattern::WorstPattern_max_pat() {
	double max_var_pat, i_sum;
	int i_max;

	max_var_pat = 0.0;
	i_max = 0;

	for(int i = 0; i < size; i++) {
		i_sum = 0.0;
		for(int j = 0; j < size; j++) {				/*If the summation for Pi has the higest value, by adding the corresponding variance part of each other pattern,...*/
			i_sum+=var_sum[i][j];				/*...Pi has the highest part at the variance summation. In this case if Pi is changed mayby the variance parts to...*/
		}							/*...each other pattern is lowered, therefor the variance, which means the variance is better*/
		if(max_var_pat < i_sum) {
			max_var_pat = i_sum;
			i_max = i;
		}
	}
	return i_max;
}

/**
 * Standardmethod for improvment, using the WorstPattern_max_pat optimization
 *
 * @param limit
 *	Number of allowed improvement steps/position changes for all pattern
 */
void Pattern::Improve(int limit) {
	DoImprove(limit, false, true, false);
}

/**
 * Method for improvment, using the loop optimization, just start with
 *	first pattern, try to improve, go on to the next till end, then
 * 	start again with first Pattern.
 *
 * @param limit
 *	Number of allowed improvement steps/position changes for all pattern
 */
void Pattern::ImproveLoop(int limit) {
	DoImprove(limit, false, false, true);
}

/**
 * Method for improvment, using the WorstPattern_max_val optimization
 *
 * @param limit
 *	Number of allowed improveement steps/position changes for all pattern
 */
void Pattern::ImproveMaxValue(int limit) {
	DoImprove(limit, true, false, false);
}

/**
 * Method for improvement, using both the WorstPattern_max_pat and
 *	WorstPattern_max_val optimization. Might be slow.
 *
 * @param limit
 *	Number of allowed improvement steps/position changes for all pattern
 */
void Pattern::ImproveMaxValuePattern(int limit) {
	DoImprove(limit, true, true, false);
}

/**
 * Activating secure mode: rebuild best_patternset, if improved pattern is
 * 	not better than best_patternset, for each improvement step.
 */
void Pattern::ImproveSecure() {
	this->secure = true;
}

/**
 * The improvement method, using the submitted booleans for correct used improvement
 * 	option.
 * @param limit
 *	Number of allowed improvement steps/position changes for all pattern
 *
 * @param max_val
 * 	Possible for the WorstPattern_max_val mode
 *
 * @param max_pat
 *	Possible for the WorstPattern_max_pat mode
 *
 * @param loop
 *	Possible for the loop mode
 */
void Pattern::DoImprove(int limit, bool max_val, bool max_pat, bool loop) {
	std::string patsave1_new, patsave1_old, patsave2_new, patsave2_old, least_pattern;
	double tmp_variance, tmp_best_variance, var1, var2;
	int worst_max_val, worst_max_pat, least_pos, steps, pat_modulo, counter_best_pat, better_pattern;
	bool flag_better;

	tmp_best_variance = GetBestVariance();
	tmp_variance = GetVariance();
	better_pattern = 0;
	counter_best_pat = 0;
	pat_modulo = counter_best_pat % size;
	steps = 0;

	if(!silent && quiet) {
		//std::cout << "First variance: \t" << GetBestVariance() << std::endl;
		//std::cout << "First norm_variance: \t" << GetBestNormVariance() << std::endl << std::endl;
	}

	if(improve) {
		for(int i = 1; i <= limit; i++) {
			if(!loop) {								/*!!!--> have a look at: 	int pattern::WorstPattern_max_pat()	and	int pattern::WorstPattern_max_val()	*/
				flag_better = false;
				if(max_val && max_pat) {						/*Possible to optimize by using both options*/
					worst_max_val = GetWorstPatMaxVal();			/*The maybe worst pattern*/
					worst_max_pat = GetWorstPatMaxPat();			/*The maybe worst pattern*/
					patsave1_old = pattern_set[worst_max_val];
					patsave2_old = pattern_set[worst_max_pat];		/*Saving them to undo this step and have the previous pattern for correct calculation*/

					ChangePatternRandom(worst_max_val);
					patsave1_new = pattern_set[worst_max_val];		/*Get new changed pattern at worst position max_value...*/
					var1 = Variance();

					pattern_set[worst_max_val] = patsave1_old;		/*... and undo it*/

					ChangePatternRandom(worst_max_pat);
					patsave2_new = pattern_set[worst_max_pat];		/*Get new changed pattern at worst position max_value...*/
					var2 = Variance();

					pattern_set[worst_max_pat] = patsave2_old;			/*... and undo it*/

					if(var1 < var2) {					/*lower variance is the better variance. Set this pattern and calculate again*/
						pattern_set[worst_max_val] = patsave1_new;
						tmp_variance = Variance();
						least_pattern = patsave1_old;
						least_pos = worst_max_val;
					}
					else{
						pattern_set[worst_max_pat] = patsave2_new;
						tmp_variance = Variance();
						least_pattern = patsave2_old;
						least_pos = worst_max_pat;
					}
				}
				else if(max_val && !max_pat) {					/*Option just to use the worst pattern by using maximum variance part value*/
					worst_max_val = GetWorstPatMaxVal();
					least_pattern = pattern_set[worst_max_val];
					least_pos = worst_max_val;
					ChangePatternRandom(worst_max_val);
					tmp_variance = Variance();
				}
				else{								/*Option just to use the worst pattern by using maximum pattern with corresponding variance parts*/
					worst_max_pat = GetWorstPatMaxPat();
					least_pattern = pattern_set[worst_max_pat];
					least_pos = worst_max_pat;
					ChangePatternRandom(worst_max_pat);
					tmp_variance = Variance();
				}
				if(tmp_best_variance > tmp_variance) {
					better_pattern++;
					if(quiet && !silent) {
					//	std::cout << "\r### BETTER PATTERN " << better_pattern <<  " ###";
					//	std::cout.flush();
					}
					else if(!quiet && !silent) {
						//std::cout << "### BETTER PATTERN " << better_pattern <<  " ###" << std::endl;
					}
					tmp_best_variance = GetBestVariance();
					flag_better = true;
				}
				if(secure && flag_better) {					/*Undo pattern change if it was not bettern than the best pattern up to now*/
					pattern_set[least_pos] = least_pattern;
					tmp_variance = Variance();
				}
			}
			else{									/*Random part, just go directly through each pattern and try to lower the variance*/
				flag_better = false;
				patsave1_old = pattern_set[pat_modulo];				/*If it is a better variance, we want to keept it and therefore go on*/
				ChangePatternRandom(pat_modulo);
				tmp_variance = Variance();

				if(steps++ > (limit/size)) {
					counter_best_pat++;
					pat_modulo = counter_best_pat % size;
					steps = 0;
				}
				if(tmp_best_variance > tmp_variance) {
					better_pattern++;
					if(quiet && !silent) {
						//std::cout << "\r### BETTER PATTERN " << better_pattern <<  " ###";
						//std::cout.flush();
					}
					else if(!quiet && !silent) {
						//std::cout << "### BETTER PATTERN " << better_pattern <<  " ###" << std::endl;
					}
					tmp_best_variance = GetBestVariance();
					counter_best_pat++;
					pat_modulo = counter_best_pat % size;			/*If it is a better variance, we want to keept it and therefore go on*/
					steps = 0;						/*Reset steps cause of patternnumber changing*/
					flag_better = true;
				}
				if(!quiet && !silent) {
					std::cout << "Step " << i << " / " << limit  << std::endl << "Patternset: \t";
					Print();
					std::cout << "variance: " << GetVariance() << std::endl;
					std::cout << "norm_variance: " << GetNormVariance() << std::endl << std::endl;
				}
				if(secure && flag_better) {
					pattern_set[pat_modulo] = least_pattern;
					tmp_variance = Variance();
				}
			}
		}
		if(!silent && quiet) {
			std::cout << "\n\nBest pattern:\t";
			std::cout << "Best variance:\t" << GetBestVariance() << std::endl;
			std::cout << "Best norm_variance: " << GetBestNormVariance() << std::endl;
		}
	}
	else{
		SecureMessage("noimprove", -1);
	}
}

/*---stuff-------------------------------------------------------------------*/
/**
 * Method to print the current pattern
 */
void Pattern::Print() {
	for(int i = 0; i < size-1; i++) {
		std::cout << pattern_set[i] << ",";
	}
	std::cout << pattern_set[size-1] << std::endl;
}

/**
 * Method calculate the number of all pattern combinations
 * Actually the gauss summation (n*(n+1)/2)
 *
 * @return the number of maximum pattern combinations
 */
double Pattern::Gauss() {
	return 0.5*size*size + 0.5*size;
}

/**
 * Determines the maximum number of patterns with weight and length
 *
 *@param weight
 *	Complete pattern weight-2, start and end have to be match and do not change
 *
 *@param length
 *	Complete pattern lengt-2, start and end have to be match and do not change
 *
 *@return max number of possible pattern
 */
double Pattern::MaxNumberPattern(int p_weight, int p_length) {
	double tmpa = Faculty(p_length);
	double tmpb = Faculty(p_length-p_weight);
	double tmpc = Faculty(p_weight);
	double tmp = tmpa / (tmpb*tmpc);
	if (int(tmp) < 0) {
		return std::numeric_limits<int>::max();
	}
	return tmp;
}

/**
 * Calculates the faculty
 *
 *@param value
 *	Calculating the faculty for value
 *
 *@return faculty
 */
double Pattern::Faculty(int value) {
	double tmp;
	tmp = 1;
	for(int i = 1; i <= value; i++) {
		tmp *=i;
	}
	return tmp;
}

void Pattern::Quiet() {
	this->quiet = true;
}

void Pattern::Silent() {
	this->quiet = true;
	this->silent = true;
}

/**
 * A Method to collect all errormessages. Just easier for programmer to change
 *  	the text or extend.
 *
 * @param errmsg
 *	Due to a few possible errormessages, this is the option, which has to be printed.
 *
 * @param pos
 *	The position of the incorrect patterns.
 *
 */
void Pattern::SecureMessage(std::string errmsg, int pos) {
	if(errmsg == "wrongindex") {
		std::cerr << "ERROR! Pattern " << pos << " does not exist... do nothing\n" << std::endl;
		return;
	}
	if(errmsg == "nkl") {
		std::cerr << "ERROR! Wrong values for weight, pattern number or pattern length!" << std::endl;
		std::cerr << "Return to default values\n" << std::endl;
		return;
	}
	if(errmsg == "weight_pat") {
		std::cerr << "ERROR! Weight of a pattern cannot be above the pattern length!" << std::endl;
		std::cerr << "Return to submittet or default values\n" << std::endl;
		return;
	}
	if(errmsg == "max_number_pattern") {
		std::cerr << "Using your pattern conditions, we can create all possible pattern directly!" << std::endl;
		std::cerr << "Updating your number of patterns to n = " << pos << "\n" << std::endl;
		return;
	}
	if(errmsg == "noimprove") {
		std::cerr << "Using your pattern conditions it is not sensible to improve your pattern, sorry!" << std::endl;
		std::cerr << "Deactivating improve mode\n" << std::endl;
		return;
	}
}
