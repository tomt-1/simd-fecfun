//define structure for simulation parameters

struct sim_struct {
	unsigned outer_count;		//outer loop count (results reported evey outer loop iteration)
	unsigned loop_count;		//inner loop count
	unsigned num_codeword;		//number of codewords per inner loop (allocates block of memory)
	unsigned max_iter;			//maximum number of LDPC decode iterations per codeword
	float max_row_metric;		//max(abs(row_metric)).  Or delta from 1 if sum-prod.
	int puncture_8023ca_flag;	//1 if puncturing the last 512 parity bits of 802.3ca matrix
	int use_given_rand_state;	//1 if using fixed random state.  (Otherwise uses rdrand)
	float LLR_scaling_factor;	//LLR scaling factor (for fixed point metrics. Use 1.0 for floats)
	float max_abs_LLR;			//maximum allowed a priori LLR
	float beta_offset;			//value for beta offsets (when using min-sum and compiled with USE_BETA_OFFSET)
	float SNR_start;			//Starting SNR value for simulations
	float SNR_end;				//Ending SNR value
	float SNR_step;				//Step size.  Don't let total number of steps exceed MAX_SNR_ITER
	uint64_t rand_state[4*4*4];	//random state for 4 SIMD PRNGs, each having 4 SIMD words of state
	char enc_raw_file[128];
	char dec_raw_file[128];
};
