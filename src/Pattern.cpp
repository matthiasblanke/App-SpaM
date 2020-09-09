/**
 * This programm calculates the variance of a set of pattern with the same length and weight.
 * It is possible to improve your patternset, estimate values for p, q and S, l_hom, li, lj
 *   from a multiple alignment file in fasta format, and also read patterns from a file.
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

/*===Main-Part===============================================================*/
/*---Constructor-& Init------------------------------------------------------*/
/**
 * Default constructor, sets the default vaulues, pattern will be generated automatically.
 */
Pattern::Pattern() {
	this->size = 10;
	this->length = 14;
	this->weight = 8;
	Pattern(NULL, NULL, size, length, weight, 10000, 10000, 10000, 0.9, 0.25);
}

/**
 * File constructor, sets values only from files; resets automatically if there are problems.
 */
Pattern::Pattern(char* pattern_file, char* align_file) {
	this->pattern_file = pattern_file;
	this->align_file = align_file;
	Pattern(pattern_file, align_file, 10, 14, 8, 10000, 10000, 10000, 0.9, 0.25);
}

/**
 * Long constructor, sets the values; resets automatically, if there are problems.
 *
 * @param pattern_file
 *		File, that may contains submitted pattern
 *
 * @param align_file
 * 		File, that may contains an alignment file to estimate p, q, l_hom, li and lj
 *
 * @param size
 * 		The amount of patterns; pattern number.
 *
 * @param length
 *		The pattern length for each pattern of the pattern set.
 *
 * @param weigth
 * 		The weight (match positions; '1') for each pattern of the pattern set.
 *
 * @param l_hom
 *		In theory, the amount of homologous positions of two sequences in an multiple alignment.
 *
 * @param l1
 * 		In theory, the first(represents each sequence i) sequence length of two observed sequences.
 *
 * @param l2
 * 		In theory, the second(represents each sequence j) sequence length of two observed sequences.
 *
 * @param p
 * 		The match probability ( = #matches / #l_hom)
 *
 * @param q
 * 		The background probability for each nucleotide A,C,G,T
 */
Pattern::Pattern(char* pattern_file, char* align_file, int size, int length, int weight, int l_hom, int l1, int l2, double p, double q) {
	this->pattern_file = pattern_file;
	this->align_file = align_file;
	this->size = size;
	this->length = length;
	this->weight = weight;
	this->l_hom = l_hom;
	this->l1 = l1;
	this->l2 = l2;
	this->p = p;
	this->q = q;
	this->variance = 0;
	this->best_variance = 0;
	this->quiet = false;
	this->silent = false;
	this->secure = false;
	ReinitPattern();
}

Pattern::Pattern(char* pattern_file, char* align_file, int size, int length, int weight, int l_hom, int l1, int l2, double p, double q, int seed) {
	this->pattern_file = pattern_file;
	this->align_file = align_file;
	this->size = size;
	this->length = length;
	this->weight = weight;
	this->l_hom = l_hom;
	this->l1 = l1;
	this->l2 = l2;
	this->p = p;
	this->q = q;
	this->variance = 0;
	this->best_variance = 0;
	this->quiet = false;
	this->silent = false;
	this->secure = false;
	generator.seed(seed);
	ReinitPattern();
}

/**
 * Short constructor, sets some default vaulues, just pattern dimension is set.
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
Pattern::Pattern(int size, int length, int weight) {
	Pattern(NULL,NULL,size,length,weight, 10000, 10000, 10000, 0.9, 0.25);
}

/**
 * Default destructor, deletes all vectors and matrices in the object.
 */
Pattern::~Pattern() {
	for (int i = 0; i < seq_matrix.size(); i++) {
		l_hom_val[i].clear();
		p_values[i].clear();
		seq_matrix[i].clear();
	}
	l_hom_val.clear();
	p_values.clear();
	seq_matrix.clear();
	seq_leng.clear();

	for (int i = 0; i < q_values.size(); i++) {
		q_values[i].clear();
	}
	q_values.clear();

	for (int i = 0; i < var_sum.size(); i++) {
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
	std::ifstream patternfile, check_align;
	std::vector<std::string> pattern_tmp;
	std::string tmp;
	char tokens[4] = {'.',' ',',',';'};			/*These tokens are allowed to seperate patterns*/
	bool start, flag_align;
	size_t f_size;

	start = false;
	flag_align = false;

	pattern_set.clear();
	best_pattern.clear();


	if(pattern_file != NULL) {
		patternfile.open(pattern_file);
		patternfile.seekg(0,std::ios::end);
    		f_size = patternfile.tellg();
		if (!patternfile) {
			SecureMessage("file", -1);
		} 
		else if (f_size == 0) {				/*[FILE].eof() does not recognize empty files...*/
			SecureMessage("empty", -1);
		}
		else {
			patternfile.close();			/*..therefore it has also to be closed and opened -.-** */
			patternfile.open(pattern_file);
			//std::cout << "Reading pattern from submitted patternfile ...\n" << std::endl;
			patternfile >> tmp;
			while(!patternfile.eof()) {
				if(!ValidatePatternsFormat(tmp)) {
					patternfile >> tmp;	/*Ignoring each incorrect pattern. It is easier to calculate with all the rests*/
				}
				else {
					pattern_tmp = SplitString(tmp, tokens);
					for(unsigned int i = 0; i < pattern_tmp.size(); i++) {
						pattern_set.push_back(pattern_tmp[i]);			/*Each pattern needs to be saved in its own std::string for comparison*/
						best_pattern.push_back(pattern_tmp[i]);
					}
					patternfile >> tmp;
				}
			}
			if (size > 0) {				/*For an empty set we do not have to validate*/
				start = ValidatePatternConditions();
			}
			if (!start) {				/*Some conditions, like changed weight is too much much amount of work*/
				SecureMessage("pattern", -1);	/*Therefore we just use our default values*/
				this->size = 10;
				this->weight = 8;
				this->length = 14;
				pattern_set.clear();		/*Do not forget to reset, or the pattern set will not replaces, just increased*/
				best_pattern.clear();
			}
			else {
				this->size = pattern_set.size();
				this->length = pattern_set[0].length();
				this->weight = PatternWeight(pattern_set[0]);
				//std::cout << "\n... Done!\n" << std::endl;
			}
		}
		patternfile.close();
	}

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
	if (size >= MaxNumberPattern(weight-2, length-2)) {	/*Match positions at the start and end do not alterate, therefore -2*/
		size = MaxNumberPattern(weight-2, length-2);
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
		for(int i = 0; i < size; i++) {
			ChangePatternRandom(i);			/*Creating just uniq Pattern!*/
		}
		this->best_pattern = PatternCopy(pattern_set);
		//std::cout << "\n ... Done!\n" << std::endl;
	}


	if(align_file != NULL) {					/*Mayby alignfile set, but if it doesnt work, we want to use the normal... */
		check_align.open(align_file);
		check_align.seekg(0,std::ios::end);		/*... variance calculation. There will be no sequences matrix*/
    		f_size = check_align.tellg();
		if(!check_align) {
			SecureMessage("align", -1);
			flag_align = false;
		}
		else if(f_size == 0) {
			SecureMessage("empty", -1);
			flag_align = false;
		}
		else{
			check_align.close();
			check_align.open(align_file);
			check_align >> tmp;
			if(tmp[0] != '>') {
				SecureMessage("fasta", -1);
				flag_align = false;
			}
			else{
				flag_align = true;		/*shows the used start variance calculation */
			}
		}
		check_align.close();
	}


	InitMatrix();

	if(!flag_align) {
		this->variance = CalcVariance();
		this->best_variance = variance;
		this->best_pattern = PatternCopy(pattern_set);
	}
	else{
		ReadAlign();
		InitPValues();					/*Neccessary just with an alginment file...*/
		InitQValues();					/*... otherwise too much RAM and junk*/
		this->variance = CalcVarianceAlign();
		this->best_variance = variance;
		this->best_pattern = PatternCopy(pattern_set);

	}

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
std::string Pattern::GetPattern(int number) {
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
std::string Pattern::GetBestPattern(int number) {
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
 * Returns the length of the homologous sequence pair
 *
 * @return returns homologous positions
 */
int Pattern::GetLHom() {
	return l_hom;
}

/**
 * Returns the length of the first observed sequence
 *
 * @return returns length sequence 1
 */
int Pattern::GetL1() {
	return l1;
}

/**
 * Returns the length of the second observed sequence
 *
 * @return returns length sequence 2
 */
int Pattern::GetL2() {
	return l2;
}

/**
 * Returns the match probability
 *
 * @return returns p value
 */
double Pattern::GetP() {
	return p;
}

/**
 * Returns the background probabillity, summation over all nucleotids
 *
 * @return returns q value
 */
double Pattern::GetQ() {
	return q;
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


/*---PatternCreateFunc-------------------------------------------------------*/
/**
 * Splits a string, read by a pattern file. Mayby each pattern does not get
 * 	a new line, it has to be parsed, when a pattern starts and ends
 *
 * @param pattern_split
 * 		The string containing a few pattern
 *
 * @param tokens
 *		The allowed tokens which can be used to seperate patterns in a line
 */
std::vector<std::string> Pattern::SplitString(std::string pattern_split, char* tokens) {
	std::vector<std::string> patternset;
	std::string tmp = "";
	bool flag_token = false;

	for(unsigned int i = 0; i < pattern_split.length(); i++) {
		for(unsigned int j = 0; j < strlen(tokens); j++) {
			if(pattern_split[i]==tokens[j]) {
				flag_token = true;
			}
		}
		if(flag_token) {							/*token found, which means in one line more patterns  --> start new pattern, save last pattern*/
			patternset.push_back(tmp);
			tmp = "";
		}
		else{
			tmp = tmp + pattern_split[i];					/*concatenating patternparts*/
		}
		flag_token = false;
	}
	patternset.push_back(tmp);
	return patternset;
}

/**
 * Validates a pattern, if it contains only pattern symbols and seperation tokens
 *
 * @param pattern_form
 * 		The pattern which has to be investigate for symbols and tokens
 *
 * @return Returns true if this one pattern is in right format, false else
 */
bool Pattern::ValidatePatternsFormat(std::string pattern_form) {
	bool flag = true;

	for(unsigned int i = 0; i < pattern_form.length(); i++) {			/*allowed tokens in patternformat, also separating tokens*/
		if(pattern_form[i] != '1' && pattern_form[i] != '0' && pattern_form[i] != ' ' && pattern_form[i] != ',' && pattern_form[i] != '.' && pattern_form[i] != ';') {
			flag=false;
			SecureMessage("format", -1);
			break;
		}
	}
	if(pattern_form[0] != '1' || pattern_form[pattern_form.length()-1] != '1') {
		flag=false;
		SecureMessage("startend", -1);
	}

	return flag;
}

/**
 * Validates a pattern set, if all patterns have the same length and weight
 *
 * @return Returns true if this one pattern is in right format, false else
 */
bool Pattern::ValidatePatternConditions() {
	int com_length, com_weight, com_size, leng;
	bool condition;

	com_size = pattern_set.size();
	com_length = pattern_set[0].length();
	com_weight = PatternWeight(pattern_set[0]);					/*as fixpoint saving first patternweight*/
	condition = true;

	for(int i = 1; i < com_size; i++) {
		leng = pattern_set[i].length();
		if(leng != com_length) {					/*for the variance each pattern should have the same weight*/
			SecureMessage("size", i);
			condition = false;
		}
		if(PatternWeight(pattern_set[i]) != com_weight) {
			SecureMessage("weight", i);
			condition = false;
		}
	}
	return condition;
}

/**
 * Estimates the weight of a pattern
 *
 * @param pattern_str
 * 		The pattern which has to be investigate for the weight
 *
 * @return Returns true if this one pattern is in right format, false else
 */
int Pattern::PatternWeight(std::string pattern_wght) {
	int str_weight;
	int str_length;

	str_weight = 0;
	str_length = pattern_wght.length();

	for(int i = 0; i < str_length; i++) {
		if(pattern_wght[i] == '1') {
			str_weight++;
		}
	}
	return str_weight;
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
	int position, counter;
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
	for(unsigned int i = 0; i < old_pattern.size(); i++) {
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
void Pattern::ChangePatternRandom(int number) {
	int pos1, pos2;
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
bool Pattern::UniqPattern(int number) {
	bool uniq = true;
	for(int i = 0; i < size; i++) {
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
	if(seq_matrix.size() < 2) {							/*At least 2 sequences are neccessary for an alignment!*/
		return CalcVariance();
	}
	else{
		return CalcVarianceAlign();
	}
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
			var_sum[i][j]=(l_hom - length + 1)*var_hom + (l1 - length + 1)*(l2 - length)*var_bac;		/*For each pair Pi and Pj this is the direct share of the complete variance...*/
			var_sum[j][i]=var_sum[i][j];									/*...which we can use to estimate a "worst" pattern*/
											/*... and save it in an size x size matrix*/
			homologue += var_hom;
			background += var_bac;

		}
	}
	this->variance = (l_hom - length + 1)*homologue + (l1 - length + 1)*(l2 - length)*background;
	if(this->variance < this->best_variance) {
		this->best_pattern = PatternCopy(pattern_set);
		this->best_variance = variance;
	}
	return variance;

}


/**
 * Calculates the variance for a pattern set with alignment file
 * In this case l_hom, l1, l2, p and q are estimated
 *
 * @return Calculates and returns current variance
 */
double Pattern::CalcVarianceAlign() {
	double var_value;
	int s_size, counter;

	s_size = seq_matrix.size();
	counter = 0;
	var_value = 0.0;

	for(int i = 0; i < s_size-1; i++) {
		for(int j = i+1; j < s_size; j++) {
			q = 0.0;
			for(int k = 0; k < 4; k++) {
				q += q_values[i][k]*q_values[j][k];
			}
			l_hom = l_hom_val[i][j];		/*Previously estimated ....*/
			p = p_values[i][j];			/*... so that by improving ...*/
			l1 = seq_leng[i];			/*it is already calculated and does not take ane time for another caculation*/
			l2 = seq_leng[j];
			var_value += CalcVariance();
			counter++;
		}
	}
	variance = var_value / counter++;
	if(variance < best_variance) {
		best_pattern = PatternCopy(pattern_set);
		best_variance = variance;
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
			PrintBest();
			std::cout << "Best variance:\t" << GetBestVariance() << std::endl;
			std::cout << "Best norm_variance: " << GetBestNormVariance() << std::endl;
		}
	}
	else{
		SecureMessage("noimprove", -1);
	}
}


/*---Alignment---------------------------------------------------------------*/
/**
 * Reads and creates sequence matrix from a submitted multiple alignment file
 */
void Pattern::ReadAlign() {
	std::ifstream alignfile;
	std::vector<std::string> tmp_seq;
	std::string tmp;

	alignfile.open(align_file);
	if(!alignfile.eof()) {
		alignfile >> tmp;
		if(tmp[0] == '>') {
			alignfile >> tmp;						/*First line is a Fasta-Header!*/
		}
		else{
			return;
		}
		while(!alignfile.eof()) {						/*seq_matrix contains std::vector, in which each line is saved as std::string*/
			if(tmp[0] == '>') {						/*up to the next fasta-header it is the same sequence*/
				seq_matrix.push_back(tmp_seq);
				tmp_seq.clear();					/*It is easier to handle the sequence, instead of concatenating everything*/
			}
			else{
				tmp_seq.push_back(tmp);
			}
			alignfile >> tmp;
		}
		seq_matrix.push_back(tmp_seq);
	}
	alignfile.close();
	return;
}

/**
 * Estimates li/lj for the variance calculation with an alignment file by increment over
 *	each position without a gap '-'
 *
 * @param seq
 *		For this sequence the previous (not alignment) length is calculated
 *
 * @return returns length of the sequence without gap positions
 */
int Pattern::LengthSeq(std::vector<std::string> seq) {
	std::string read;
	int s_size, s_length;

	s_size = seq.size();
	s_length = 0;

	for(int i = 0; i < s_size; i++) {
		read = seq[i];
		for(unsigned int j = 0; j < read.length(); j++) {
			if(read[j] != '-') {
				s_length++;					/*Everything instead of a gap belongs to the original sequences length*/
			}
		}
	}
	return s_length;
}
/**
 * Estimates all values for q for the variance calculation with an alignment file for each
 * 	sequence by scanning for each nucleotide and dividing by all nucleotides
 */
void Pattern::InitQValues() {
	int s_size;

	s_size = seq_matrix.size();

	for(int i = 0; i < s_size; i++) {
		q_values.push_back(BackgroundProb(seq_matrix[i]));
	}
}

/**
 * Estimates all values for p for the variance calculation with an alignment file for each
 * 	sequence pair by estimating and saving each homologous sequence pair length and estemating
 * 	pairwise each number of match positions
 *
 * Estimates and saves each original sequence length
 */
void Pattern::InitPValues() {
	std::vector<double> tmp_p;
	std::vector<int> tmp_hom;
	double match_tmp;
	int s_size, l_hom_tmp, s_leng;

	s_size = seq_matrix.size();

	for(int i = 0; i < s_size-1; i++) {
		for(int j = 0; j < s_size; j++) {
			tmp_p.push_back(0.0);
			tmp_hom.push_back(0);
		}						/*Was neccessary due to a bug, where using push_back directly changes all variables if used...*/
	}
	for(int i = 0; i < s_size; i++) {
		l_hom_val.push_back(tmp_hom);			/*... in the for loop below.*/
		p_values.push_back(tmp_p);
	}


	for(int i = 0; i < s_size-1; i++) {
		for(int j = i+1; j < s_size; j++) {
			l_hom_tmp = CountHom(i, j);
			l_hom_val[i][j] = l_hom_tmp;
			match_tmp = CountMatch(i, j);
			p_values[i][j] = match_tmp / ((double)l_hom_tmp);	/*probability of matches for the sequence pair Si and Sj*/
		}
	}

	for(int i = 0; i < s_size; i++) {
		s_leng= LengthSeq(seq_matrix[i]);
		this->seq_leng.push_back(s_leng);
	}

	return;
}

/**
 * Estimates one specific background probability q for a sequence
 *	q[0] = A , q[1] = C, q[2] = G, q[3] = T
 *
 * @param sequence
 *		For this sequence all background probabilities for each nucleotide is calculated
 *
 * @return returns a vector containing for each nucleatide the probability
 */
std::vector<double> Pattern::BackgroundProb(std::vector<std::string> sequence) {
	std::vector<double> q(4,0.0);
	std::string read;
	int s_length;

	s_length = 0;

	for(unsigned int i = 0; i < sequence.size(); i++) {
		read = sequence[i];						/*For each sequencepart we count the number of every nucleotide*/
		for(unsigned int j = 0; j < read.length(); j++) {
			if(read[j] == 'A' || read[j] == 'a') {
				q[0]++;
				s_length++;
			}
			else if(read[j] == 'C' || read[j] == 'c') {
				q[1]++;
				s_length++;
			}
			else if(read[j] == 'G' || read[j] == 'g') {
				q[2]++;
				s_length++;
			}
			else if(read[j] == 'T' || read[j] == 't') {
				q[3]++;
				s_length++;
			}
		}
	}
	for(int i = 0; i < 4; i++) {
		q[i] = q[i] / (double) s_length; 					/*Sequencelength (NOT alignment) is used to get the relative values*/
	}
	return q;
}


/**
 * Estimates for two specific sequences the common match positions
 *
 * @param pos1
 *		Number of the position 1 in the pattern set
 *
 * @param pos2
 *		Number of the position 2 in the pattern set
 *
 * @return returns the number of common matches
 */
double Pattern::CountMatch(int pos1, int pos2) {
	std::string read1, read2;
	double count;
	int s_size;

	count = 0;
	s_size = seq_matrix[pos1].size();

	for(int i = 0; i < s_size; i++) {
		read1 = seq_matrix[pos1][i];
		read2 = seq_matrix[pos2][i];
		for(unsigned int j = 0; j < read1.length(); j++) {
			if(read1[j] != '-' && read1[j] == read2[j]) {
				count++;
			}
		}
	}
	return count;
}

/**
 * Estimates for two specific sequences the homologous positions
 * 	where is not a gap, but different nucleotides are allowed
 *
 * @param pos1
 *		Number of the position 1 in the pattern set
 *
 * @param pos2
 *		Number of the position 2 in the pattern set
 *
 * @return returns the number of homologous positions
 */
int Pattern::CountHom(int pos1, int pos2) {
	std::string read1, read2;
	int count;
	int s_size;

	count = 0;
	s_size = seq_matrix[pos1].size();

	for(int i = 0; i < s_size; i++) {					/*Homologous when there is no gap in the alignment*/
		read1 = seq_matrix[pos1][i];
		read2 = seq_matrix[pos2][i];
		for(unsigned int j = 0; j < read1.length(); j++) {
			if(read1[j] != '-' && read2[j] != '-') {
				count++;
			}
		}
	}
	return count;
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
 * Method to print the current best pattern
 */
void Pattern::PrintBest() {
	for(int i = 0; i < size-1; i++) {
		//std::cout << best_pattern[i] << ",";
	}
	//std::cout << best_pattern[size-1] << std::endl;
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
	if(errmsg == "file") {
		std::cerr << "ERROR! Pattern file \'" << pattern_file << "\' could not be found!" << std::endl;
		std::cerr << "Return to submitted or default values\n" << std::endl;
		return;
	}
	if(errmsg == "empty") {
		std::cerr << "ERROR! File \'" << pattern_file << "\' is an empty file!" << std::endl;
		std::cerr << "Return to submitted or default values\n" << std::endl;
		return;
	}
	if(errmsg == "pattern") {
		std::cerr << "ERROR! Patternconditions from pattern file were not correct (different weight or length)!" << std::endl;
		std::cerr << "Return to submitted or default values\n" << std::endl;
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
	if(errmsg == "format") {
		std::cerr << "FORMAT-ERROR: Pattern containes illegal characters!" << std::endl;
		std::cerr << "Allowed characters: '0','1' for pattern; ','|'.'|';'|' ' to seperate patterns" << std::endl;
		std::cerr << "Go on to next pattern.\n" << std::endl;
		return;
	}
	if(errmsg == "max_number_pattern") {
		std::cerr << "Using your pattern conditions, we can create all possible pattern directly!" << std::endl;
		std::cerr << "Updating your number of patterns to n = " << pos << "\n" << std::endl;
		return;
	}
	if(errmsg == "align") {
		std::cerr << "ERROR! Alignment file \'" << align_file << "\' could not be found!" << std::endl;
		std::cerr << "Return to submitted or default values\n" << std::endl;
		return;
	}
	if(errmsg == "fasta") {
		std::cerr << "ERROR! Alignment file \'" << align_file << "\' seems not be in fasta format!" << std::endl;
		std::cerr << "Return to submitted or default values\n" << std::endl;
		return;
	}
	if(errmsg == "size") {
		std::cerr << "By comparing with the first pattern, the " << pos+1 << ". pattern has a different length!\n" << std::endl;
		return;
	}
	if(errmsg == "weight") {
		std::cerr << "By comparing with the first pattern, the " << pos+1 << ". pattern has a different weight!\n" << std::endl;
		return;
	}
	if(errmsg == "startend") {
		std::cerr << "FORMAT-ERROR: Pattern has to start and end with a match position '1' !\n" << std::endl;
		return;
	}
	if(errmsg == "noimprove") {
		std::cerr << "Using your pattern conditions it is not sensible to improve your pattern, sorry!" << std::endl;
		std::cerr << "Deactivating improve mode\n" << std::endl;
		return;
	}
}
