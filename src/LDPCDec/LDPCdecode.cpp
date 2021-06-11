#include <stdio.h>
// scalars
// majrow, majcol, z_value, max_row_weight
// total Hc ones (ones in compressed H matrix)
// max_iterations
// max_row_metric
// group_cnt
// (sum of f2b copy counts?) - can be obtained by summing while reading
// (sum of b2f copy coutns?) - can be obtained by summing while reading

// arrays
// row_weights [majrow] - row weight for each majrow (currently "act_cols")
// row_active [total Hc ones] - rows active in each majrow (currently "blk_tbl")
// row_offset [total Hc ones]

// f2b_copy_cnt [mrow]
// f2b_cols [sum f2b copy cnt]
// f2b_offsets [sum f2b copy cnt]

// b2f_copy_cnt [mrow]
// b2f_cols [sum b2f copy cnt]
// b2f_offsets [sum b2f copy cnt]

// grp_ary [max_row_weight]

// SIMD pointers/arrays
// apriori_metrics[ MCOL * (SIMD_PER_Z) * WIDTH/8 ] - read only - copy these to decoded_metrics
// decoded_metrics[ MCOL * (SIMD_PER_Z+2) * WIDTH/8 ]
//       was "global_metrics" - this is main return value
//       caller knows that there is an ignored SIMD word before/after each major block
// local_metrics [ total_H_ones * SIMD_PER*Z * WIDTH/8 ]
// group_min[GROUP_CNT][SIMD_PER_Z]  -- restructure as one dimension (GROUP_CNT is usually MAX_ROW_WEIGHT)
// min_grp[GROUP_CNT][SIMD_PER_Z] -- restructure as one dimension
// row_xor[SIMD_PER_Z]
// combine group_min, min_grp, and row_xor into single work_area interface?

// min_function_ptr (if ONLY_ONE_MIN_FUNCTION)
// pointer to list of min_functions (if not ONLY_ONE_MIN_FUNCTION)

//pull all f2b and b2f data into required structure:
// sum f2b copy cnt, sum b2f copy cnt, f2b_copy_cnt[mrow], b2f_copy_cnt[mrow],
// etc

//TBD: also want input for check syndrome at iteration 0 (not high priority -- useless for 802.3ba)
//also, might want option to check syndrome at fractional iteration

//putting ??? on things that decoder might not need
//int return can be actual iterations
//
//compile flag for: ALWAYS_CHECK_SYNDROME (ignore syndrome_check_array or remove from interface)
// one_min_function - just use one min function regardless of row weight - changes interface?
// caller_copies_metrics - apriori_metrics are already copied to decoded_metrics by caller
// no_grouping_min_vals

// for now, remove grouping options.  will need multiple grouping patterns and array for which ones
// to use on which iteration.  Expect to use grouping in early iterations, then no grouping

#include "vectorclass.h"
#include "vectorclass_extensions.h"

#include "simd_defs.h"

#include "dec_struct_LDPC.h"

const int caller_copies_metrics = 1;
const int zero_local_metrics = 1;
const int always_check_syndrome = 1;

//including subfunctions via "include", just in case compiling as separate object limits optimization
#include "decode_functions.cpp"

typedef void (*min_fn_ptr_type)(SIMD_CLASS *,SIMD_CLASS *,SIMD_CLASS,int);

//collapse matrix params into struct: majrow, majcol, z_value, total_Hc_ones, *row_weights
//     *row_active, *row_offset, f2b_ pointers, b2f_ pointers
//could also collapse decoded_metrics_all, local_metrics, group_min, fin_min, row_xor pointers
//no real need for explicit SIMD struct -- kind of want SIMD_CLASS-related item as arg for possible overloading

int LDPCdecode ( const unsigned majrow, const unsigned majcol, const unsigned z_value,
	const unsigned total_Hc_ones, const unsigned max_iter,
	const float max_row_metric, const float beta_offset,
	unsigned *row_weights, unsigned *row_active, unsigned *row_offset,
	//could pass *f2b_b2f_data with below structure to reduce interface count
	unsigned *f2b_copy_cnt, unsigned *f2b_cols, unsigned *f2b_offsets,
	unsigned *b2f_copy_cnt, unsigned *b2f_cols, unsigned *b2f_offsets,
	SIMD_CLASS *apriori_metrics,
	SIMD_CLASS *decoded_metrics_all, SIMD_CLASS *local_metrics,
	SIMD_CLASS *group_min, SIMD_CLASS *fin_min, SIMD_CLASS *row_xor,
	min_fn_ptr_type *min_funcs
	//could add syndrome check rules (per iteration)
	//could add grouping counts and grouping lists per iteration
) {

	unsigned SIMD_PER_Z = (z_value+SIMD_PAR-1) / SIMD_PAR; 	//ceiling function
	int Z_remainder = z_value - (SIMD_PER_Z-1)*SIMD_PAR;	//utilized elements in last SIMD word of each Z col/row
	SIMD_INT AllOnes = -1;
	SIMD_CLASS Z_remainder_mask = 0;
	Z_remainder_mask.load_partial(Z_remainder,(CTYPE *)&AllOnes );	//bit-mask of utilized elements in last SIMD word

	SIMD_CLASS *decoded_metrics = decoded_metrics_all + 1; 	//1st SIMD word used for partial copies
	if ( !caller_copies_metrics ) {
		SIMD_CLASS *apriori_ptr = (SIMD_CLASS *)apriori_metrics;
		for (unsigned dcol=0; dcol < (majcol-majrow); ++dcol) {
			for (unsigned i=0; i < SIMD_PER_Z; ++i) {
				*decoded_metrics++ = *apriori_ptr++;
			}
			decoded_metrics += 2;
		}
		decoded_metrics = decoded_metrics_all + 1;	//reset to start
	}

	//zero the local_metrics
	//have option to not zero these -- in case doing multiple 1-iteration calls
	//     and changing options, or grafting in/out interleaved codewords
	if ( zero_local_metrics ) {
		for (unsigned i=0; i < total_Hc_ones * SIMD_PER_Z; ++i) {
			local_metrics[i] ^= local_metrics[i];
		}
	}

	unsigned done = 0;
	unsigned iter = 0;
	while ( !done ) { //do an iteration if not done
		unsigned *row_active_ptr = (unsigned *)row_active;
		int *row_offset_ptr = (int *)row_offset;
		SIMD_CLASS *local_metrics_ptr = local_metrics;
		CTYPE *dec_metric_elem_ptr = (CTYPE *)decoded_metrics;

		unsigned f2b_idx = 0;
		unsigned b2f_idx = 0;
		for (unsigned mrow=0; mrow < majrow; ++mrow) {
			//do f2b and b2f copies for this row of the decoded metrics
			for (unsigned i=0; i < f2b_copy_cnt[mrow]; ++i) {
				unsigned act_col = f2b_cols[f2b_idx];
				unsigned offset  = f2b_offsets[f2b_idx];
				CTYPE *col_base_ptr = dec_metric_elem_ptr + (SIMD_PER_Z+2)*SIMD_PAR*act_col;
				SIMD_CLASS val_to_copy;
				val_to_copy.load(  col_base_ptr + offset );
				val_to_copy.store( col_base_ptr + offset + z_value );
				++f2b_idx;
			}
			for (unsigned i=0; i < b2f_copy_cnt[mrow]; ++i) {
				unsigned act_col = b2f_cols[b2f_idx];
				unsigned offset = b2f_offsets[b2f_idx];
				CTYPE *col_base_ptr = dec_metric_elem_ptr + (SIMD_PER_Z+2)*SIMD_PAR*act_col;
				SIMD_CLASS val_to_copy;
				val_to_copy.load(  col_base_ptr - SIMD_PAR + offset + z_value );
				val_to_copy.store( col_base_ptr - SIMD_PAR + offset );
				++b2f_idx;
			}

			//decode row
			decode_row( row_weights[mrow], row_active_ptr, row_offset_ptr, max_row_metric, beta_offset,
				row_xor, group_min, fin_min, local_metrics_ptr, decoded_metrics, min_funcs[row_weights[mrow]], 
				SIMD_PER_Z, z_value, Z_remainder );
			row_active_ptr += row_weights[mrow];
			row_offset_ptr += row_weights[mrow];
			local_metrics_ptr += row_weights[mrow] * SIMD_PER_Z;
		}
		//do b2f copies to make all columns have data at offset 0.  (b2f_copy_cnt has an extra entry for this)
		for (unsigned i=0; i < b2f_copy_cnt[majrow]; ++i) {
			unsigned act_col = b2f_cols[b2f_idx];
			unsigned offset  = b2f_offsets[b2f_idx];
			CTYPE *col_base_ptr = dec_metric_elem_ptr + (SIMD_PER_Z+2)*SIMD_PAR*act_col;
			SIMD_CLASS val_to_copy;
			val_to_copy.load(  col_base_ptr - SIMD_PAR + offset + z_value );
			val_to_copy.store( col_base_ptr - SIMD_PAR + offset );
			++b2f_idx;
		}

		//TBD: skip check if last iteration?? seems reasonable
		if ( always_check_syndrome ) { // || syndrome_check[iter] ) {
			//copy offset 0 (at front) to end
			for (unsigned i=0; i < majcol; ++i) {
				SIMD_CLASS *col_simd_ptr = decoded_metrics + (SIMD_PER_Z+2)*i;
				CTYPE *col_elem_ptr = (CTYPE *)col_simd_ptr;
				SIMD_CLASS val_to_copy = col_simd_ptr[0];
				val_to_copy.store( col_elem_ptr + z_value );
			}
			//check syndrome
			row_active_ptr = (unsigned *)row_active;
			row_offset_ptr = (int *)row_offset;
			SIMD_CLASS syn_xor = 0;
			for (unsigned mrow=0; mrow < majrow; ++mrow ) {
				SIMD_CLASS ret_val;
				ret_val = calc_syndrome_row( row_weights[mrow], row_active_ptr, row_offset_ptr, decoded_metrics, row_xor, SIMD_PER_Z, z_value, Z_remainder_mask );
				syn_xor |= ret_val;
				row_active_ptr += row_weights[mrow];
				row_offset_ptr += row_weights[mrow];
			}
			// SIMD_INT is integer version of SIMD_CLASS.  Functionally needed for float types
			SIMD_INT t1 = reinterpret_i( syn_xor );
			MTYPE synmask = to_bits( t1 < 0 );
			done = (synmask == 0);
		}
		++iter;
		done = done || (iter == max_iter);
	} //end "while (!done)"

	return iter;
}

unsigned LDPCdecode_set( unsigned num_codeword, unsigned max_iter, float max_row_metric, float beta_offset, struct Hdec_struct *Hdec, min_fn_ptr_type *min_funcs ) {
	unsigned mrow = Hdec->mrow;
	unsigned mcol = Hdec->mcol;
	unsigned z_value = Hdec->z_value;
	unsigned total_Hc_ones = Hdec->total_Hc_ones;
	unsigned *row_weights = Hdec->row_weights;
	unsigned *row_active = Hdec->row_active_cols;
	unsigned *row_offset = Hdec->row_offset;
	unsigned *f2b_copy_cnt = Hdec->f2b_copy_cnt;
	unsigned *f2b_cols = Hdec->f2b_cols;
	unsigned *f2b_offsets = Hdec->f2b_offsets;
    unsigned *b2f_copy_cnt = Hdec->b2f_copy_cnt;
    unsigned *b2f_cols = Hdec->b2f_cols;
    unsigned *b2f_offsets = Hdec->b2f_offsets;
	SIMD_CLASS *apriori_metrics = NULL;
	SIMD_CLASS *decoded_metrics_all = Hdec->decoded_metrics;
	SIMD_CLASS *local_metrics = Hdec->local_metrics;
	SIMD_CLASS *group_min = Hdec->group_min;
	SIMD_CLASS *fin_min = Hdec->fin_min;
	SIMD_CLASS *row_xor = Hdec->row_xor;

	unsigned simd_per_z = (z_value+SIMD_PAR-1) / SIMD_PAR;
	unsigned metric_increment = Hdec->mcol * (simd_per_z + 2);
	unsigned tot_iter = 0;
	for (unsigned cw=0; cw < num_codeword; ++cw) {
		//ignoring the iteration count return value for now
		unsigned iter_cnt;
		iter_cnt = LDPCdecode( mrow, mcol, z_value, total_Hc_ones,
			max_iter, max_row_metric, beta_offset, row_weights, row_active, row_offset,
			f2b_copy_cnt, f2b_cols, f2b_offsets,
			b2f_copy_cnt, b2f_cols, b2f_offsets,
			apriori_metrics, decoded_metrics_all, local_metrics,
			group_min, fin_min, row_xor, min_funcs );
		tot_iter += iter_cnt;
        //tot_iter = (iter_cnt > tot_iter) ? iter_cnt : tot_iter; //just getting max iter for now
		//DEBUG
		if ( (0) && (iter_cnt == max_iter) ) {
			printf("cw=%d iter=%d\n",cw,iter_cnt);
		}
		decoded_metrics_all += metric_increment;
	}
	return tot_iter;
}
