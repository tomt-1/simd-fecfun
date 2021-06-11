//define SIMD_CLASS outside

struct Hdec_struct {
	uint32_t z_value;		//Z value for the quasi-cyclic H matrix
	uint32_t mrow;			//Num of major rows in H matrix (rows in compressed Hc matrix)
	uint32_t mcol;			//Num of major columns in H matrix/cols in Hc matrix
	uint32_t total_Hc_ones;	//Total 1's in compressed H matrix
	uint32_t max_row_weight;	//Max 1's in a given column of full H matrix
	uint32_t *row_weights;	//weights (number of active columns) for each major row
	uint32_t *row_active_cols;	//column number of the active columns for each major row
	uint32_t *row_offset;	//offsets (from cyclic shift) for these columns in each row
	uint32_t *f2b_copy_cnt;	//count of column metrics in which to do a f2b copy
	uint32_t *f2b_cols;		//column numbers involved in doing f2b copies
	uint32_t *f2b_offsets;	//offsets within each column 
	uint32_t *b2f_copy_cnt;	//same as above, but for b2f (back-to-front) copies
	uint32_t *b2f_cols;
	uint32_t *b2f_offsets;
	SIMD_CLASS *decoded_metrics;
	SIMD_CLASS *local_metrics;
	SIMD_CLASS *group_min;
	SIMD_CLASS *fin_min;
	SIMD_CLASS *row_xor;
	uint32_t *iterations;	//decoder iterations.  Could make this 16 or 8b if desired
};

/*
//structure for items calculated from z_value (when Z != 256)
struct z_items {
	unsigned zmod8;
    unsigned zmod8_rev;
    unsigned last_byte_offset;
    unsigned last_bitpos;
    uint8_t endmask;
    uint8_t begmask;
};
*/
