#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "vectorclass.h"
#include "simd_defs.h"

#include "dec_struct_LDPC.h"

//load H matrix items for decode

void decode_matrix_allocate (const char *filename, struct Hdec_struct *H, unsigned num_cw) {

	FILE *fp;
	int fread_return;

	//Consider adding error checking on opening file and all fread return values

	//read the Decode Matrix data (Decode H matrix dimensions must match Encode)
	fp = fopen(filename,"r");

    uint32_t H_scalars[7];
    fread_return = fread(H_scalars,sizeof(int32_t),7,fp);
    H->mrow = H_scalars[0];
    H->mcol = H_scalars[1];
    H->z_value = H_scalars[2];
    H->total_Hc_ones = H_scalars[3];
    H->max_row_weight = H_scalars[4];
    uint32_t len_f2b_items = H_scalars[5];
    uint32_t len_b2f_items = H_scalars[6];

    uint32_t total_read_count = H->mrow + H->total_Hc_ones*2 + H->mrow + len_f2b_items*2 +
                        H->mrow+1 + len_b2f_items*2;

    uint32_t *H_arrays = (uint32_t *)malloc(sizeof(uint32_t)*total_read_count);
    fread_return = fread(H_arrays,sizeof(uint32_t),total_read_count,fp);
	if (fread_return != (int)total_read_count) {
		printf("Error in reading decode Matrix data\n");
		exit(1);
	}
    fclose(fp);
    H->row_weights     = H_arrays;
    H->row_active_cols = H_arrays+H->mrow;
    H->row_offset      = H_arrays+H->mrow+H->total_Hc_ones;
    H->f2b_copy_cnt    = H->row_offset+H->total_Hc_ones;
    H->f2b_cols        = H->f2b_copy_cnt+H->mrow;
    H->f2b_offsets     = H->f2b_cols + len_f2b_items;
    H->b2f_copy_cnt    = H->f2b_offsets + len_f2b_items;
    H->b2f_cols        = H->b2f_copy_cnt + H->mrow+1;
    H->b2f_offsets     = H->b2f_cols + len_b2f_items;

	//Allocate memory for items
	unsigned SIMD_PER_Z = (H->z_value+SIMD_PAR-1) / SIMD_PAR;
	//SIMD_CLASS *apriori_metrics = (SIMD_CLASS *)_mm_malloc( num_cw*mcol*SIMD_PER_Z*WIDTH/8,WIDTH/8 );
    H->decoded_metrics = (SIMD_CLASS *)_mm_malloc( num_cw*H->mcol*(SIMD_PER_Z+2)*WIDTH/8,WIDTH/8 );
    H->local_metrics = (SIMD_CLASS *)_mm_malloc( H->total_Hc_ones * SIMD_PER_Z * WIDTH/8,WIDTH/8 );
    H->group_min = (SIMD_CLASS *)_mm_malloc( H->max_row_weight * SIMD_PER_Z * WIDTH/8,WIDTH/8 );
    H->fin_min = (SIMD_CLASS *)_mm_malloc( H->max_row_weight * SIMD_PER_Z * WIDTH/8,WIDTH/8 );
    H->row_xor = (SIMD_CLASS *)_mm_malloc( SIMD_PER_Z * WIDTH/8, WIDTH/8 );
	H->iterations = (uint32_t *)malloc(sizeof(uint32_t)*num_cw);
}
