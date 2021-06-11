//Store Global summed values in a manner aligned with normal major rows/columns
//  Need to access these in un-aligned manner when doing decode_row
//Store Local values aligned such that each major row starts with the correct offset (since Local values
// are specific per-row, different rows can have different offsets and still be aligned accesses)


//USE_SAT_ADD_SUB moved to simd headers
//if integer and need saturating adds/subtracts, set this to one,
//if float or setting maxes to insure no saturation, set to zero

void decode_row( unsigned num_active_cols, unsigned *row_active_ptr, 
	int *row_offset_ptr, 
	float max_row_metric,	//could initialize group_min outside function
	float beta_offset,		//only used if USE_BETA_OFFSET is 1 during compile
	SIMD_CLASS *row_xor,
	SIMD_CLASS *group_min,
	SIMD_CLASS *fin_min,
	SIMD_CLASS *local_metrics,
	SIMD_CLASS *decoded_metrics,
	void (*min_fn_ptr)(SIMD_CLASS *,SIMD_CLASS *,SIMD_CLASS,int),
	//a given row only ever needs one min_fn, even if we use different ones on different rows
	unsigned SIMD_PER_Z,
	int Z_VALUE,
	int Z_remainder
) {

	SIMD_CLASS Max_row_simd;
	if ( max_row_metric < 1e-6 ) { //a bit of a hack.  Assume sum-product if max_row_metric is small
		Max_row_simd = 1.0 - (CTYPE)max_row_metric;
	} else {
		Max_row_simd = max_row_metric;
	}
	//initialize row XOR to zero
	for (unsigned i=0; i < SIMD_PER_Z; ++i) {
		row_xor[i] ^= row_xor[i];
	}

	//find decoded_metrics - local_val for this row.  Prepare for min-finding
	for (unsigned j=0; j < num_active_cols; ++j) {
		int col_num = row_active_ptr[j];
		SIMD_CLASS *local_blk_start = local_metrics+j*SIMD_PER_Z;
		int local_idx = row_offset_ptr[j];

		CTYPE *glb_ptr = (CTYPE *)(decoded_metrics + (SIMD_PER_Z+2)*col_num);

		for (unsigned i=0; i < SIMD_PER_Z; ++i) {
			SIMD_CLASS Global_val, loc_val, loc_residual, loc_abs;
			Global_val.load( glb_ptr + local_idx );
			loc_val = local_blk_start[i];
			#if USE_SAT_ADD_SUB
				loc_residual = sub_saturated(Global_val, loc_val);
			#else
				loc_residual = Global_val - loc_val;
			#endif
			local_blk_start[i] = loc_residual; 	//save the adjusted residual value 
			loc_abs = abs(loc_residual);
			//it is likely that option to *not* have a beta will be removed after
			//measuring the cost of doing the subtraction/max versus not having it at all
			#if USE_BETA_OFFSET
				loc_abs = loc_abs - beta_offset;	//TBD: is implicit float->int conversion handled
				loc_abs = max(loc_abs, 0.0f);
			#endif
			//TBD: handle grouping, likely by pulling 1st iteration on it's own
			group_min[j*SIMD_PER_Z + i] = loc_abs;
			row_xor[i] ^= loc_residual;
			local_idx += SIMD_PAR;
			if ( local_idx >= Z_VALUE ) local_idx -= Z_VALUE;
		}
	}

	//find minimum val applied to each group/element
	min_fn_ptr( group_min, fin_min, Max_row_simd, SIMD_PER_Z );

	//Calc new Global_val and local_val
	for (unsigned j=0; j < num_active_cols; ++j) {
		int col_num = row_active_ptr[j];
		SIMD_CLASS *local_blk_start = local_metrics+j*SIMD_PER_Z;
		int local_idx = row_offset_ptr[j];
		CTYPE *glb_ptr = (CTYPE *)(decoded_metrics + (SIMD_PER_Z+2)*col_num);

		//Option: rewrite as i=0 to SIMD_PER_Z-2.  Then final i=SIMD_PER_Z-1 with
		//store_partial on last tmp_glb store for correct number of elements
		for (unsigned i=0; i < SIMD_PER_Z; ++i) {
			SIMD_CLASS new_loc_signs, new_loc_val,tmp_glb;
			new_loc_signs = row_xor[i] ^ local_blk_start[i];
			new_loc_val = fin_min[j*SIMD_PER_Z + i];
			new_loc_val = sign_combine( new_loc_val, new_loc_signs );
			#if USE_SAT_ADD_SUB
				tmp_glb = add_saturated( local_blk_start[i], new_loc_val );
			#else
				tmp_glb = local_blk_start[i] + new_loc_val;	//at this point, local_blk_start[i] has the Global-old_local val
			#endif
			local_blk_start[i] = new_loc_val;
			if ( i != (SIMD_PER_Z-1) ) {
				tmp_glb.store( glb_ptr+local_idx );
			} else {
				tmp_glb.store_partial( Z_remainder, glb_ptr+local_idx );
			}
			local_idx += SIMD_PAR;
			if ( local_idx >= Z_VALUE ) local_idx -= Z_VALUE;
		}
	}
}

//Calculate partial syndrome for one major row.  Collapse into single SIMD word
//sign bit is zero on all elements means entire major row was zero
//assumes that global_metric has copied first SIMD_WORD to end where needed
SIMD_CLASS	calc_syndrome_row( unsigned num_active_cols, unsigned *row_active_ptr,
		int *row_offset, SIMD_CLASS *decoded_metrics, SIMD_CLASS *row_xor,
		unsigned SIMD_PER_Z, int Z_VALUE, SIMD_CLASS Z_remainder_mask ) {
	//zero out the row_xor.  
    for (unsigned i=0; i < SIMD_PER_Z; ++i) {
		row_xor[i] ^= row_xor[i]; //set to zero
    }
	for (unsigned j=0; j < num_active_cols; ++j) {
		int blk_num = row_active_ptr[j];
		int local_idx = row_offset[j];
		CTYPE *glb_ptr = (CTYPE *)(decoded_metrics + (SIMD_PER_Z+2)*blk_num);
		for (unsigned i=0; i < SIMD_PER_Z; ++i ) {
			SIMD_CLASS glb_val;
			glb_val.load( glb_ptr + local_idx );
			row_xor[i] ^= glb_val;
			local_idx += SIMD_PAR;
			if ( local_idx >= Z_VALUE ) local_idx -= Z_VALUE;
		}
	}

	//possible to have a matrix with SIMD_PER_Z equal to 1.
	SIMD_CLASS ret_val;
	ret_val ^= ret_val; //set to zero
	for (unsigned i=0; i < SIMD_PER_Z-1; ++i) {
		ret_val |= row_xor[i];
	}
	SIMD_CLASS last_or_value = row_xor[SIMD_PER_Z-1] & Z_remainder_mask;
	ret_val |= last_or_value;
	return ret_val;
}
