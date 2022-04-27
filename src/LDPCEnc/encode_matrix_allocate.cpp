#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "vectorclass.h"

#include "enc_defs.h"
#include "enc_struct_LDPC.h"

//load H matrix items for encode

//function for filling puncture and filler masks (punctured bits will have 0 LLR, filler bits will have max LLR)
int enc_fill_mask( uint32_t cnt, uint32_t *ptr_list, ENC_CLASS *ptr_mask, uint32_t Z, uint32_t mcol ) {
    ENC_CLASS *mask = ptr_mask;
	ENC_CLASS val0 = 0;
    for (uint32_t i=0; i < mcol; ++i) {
        *mask = val0;  //requires Z <= ENC_CLASS bit width
        ++mask;
    }
	int sum_fill_punc = 0;
    for (uint32_t i_cnt=0; i_cnt < cnt; i_cnt += 2) {
        int start_idx = ptr_list[i_cnt];
        int end_idx = ptr_list[i_cnt+1];
		sum_fill_punc += (end_idx-start_idx)+1;
        for (int idx=start_idx; idx <= end_idx; ++idx) {
            int col = idx / Z;
            int offset = idx - col*Z;
            int q_offset = offset/64;
            int bit_offset = offset - q_offset*64;
            uint64_t qval = (ptr_mask+col)->extract(q_offset);
            qval ^= (1ULL << bit_offset);
            (ptr_mask+col)->insert(q_offset,qval);
        }
    }
	return sum_fill_punc;
}

void encode_matrix_allocate (const char *filename, struct Henc_struct *H, unsigned num_cw) {

	FILE *fp;
	int fread_return;

	//Consider adding error checking on opening file and all fread return values

	//read the Encode Matrix data
	fp = fopen(filename,"r");
	if ( fp == NULL ) {
		fprintf(stderr,"Cannot open file %s for reading\n",filename);
		exit(1);
	}
    fread_return  = fread(&(H->z_value),sizeof(int32_t),1,fp);
	fread_return += fread(&(H->puncture_cnt),sizeof(int32_t),1,fp);
	fread_return += fread(&(H->filler_cnt),sizeof(int32_t),1,fp);
	fread_return += fread(&(H->num_par_cw),sizeof(int32_t),1,fp);
    fread_return += fread(&(H->mrow),sizeof(int32_t),1,fp);
    fread_return += fread(&(H->mcol),sizeof(int32_t),1,fp);
    if ( fread_return != 6 ) {
        fprintf(stderr,"Error in fread of initial encode integers\n");
        exit(1);
    }
    H->par_cnts = (uint32_t *)malloc(sizeof(int32_t)*H->mrow);
    H->par_colnum = (uint32_t *)malloc(sizeof(int32_t)*H->mrow);
    fread_return = fread(H->par_cnts,sizeof(int32_t),H->mrow,fp);
    fread_return = fread(H->par_colnum,sizeof(int32_t),H->mrow,fp);
    unsigned sum_cnt = 0;
    for (unsigned i=0; i < H->mrow; ++i) {
        sum_cnt += H->par_cnts[i];
    }
    H->par_offset = (uint32_t *)malloc(sizeof(int32_t)*sum_cnt);
    fread_return = fread(H->par_offset,sizeof(int32_t),sum_cnt,fp);

	int punc_fill_cnt;
	if ( H->puncture_cnt > 0 ) {
		H->punc_list = (uint32_t *)malloc(sizeof(int32_t) * H->puncture_cnt);
		fread_return = fread(H->punc_list, sizeof(int32_t),H->puncture_cnt,fp);
	}
	H->puncture_mask = (ENC_CLASS *)_mm_malloc( H->mcol*sizeof(ENC_CLASS),sizeof(ENC_CLASS) );
	punc_fill_cnt = enc_fill_mask( H->puncture_cnt, H->punc_list, H->puncture_mask, H->z_value, H->mcol );
	if ( H->filler_cnt > 0 ) {
		H->fill_list = (uint32_t *)malloc(sizeof(int32_t) * H->filler_cnt);
		fread_return = fread(H->fill_list, sizeof(int32_t),H->filler_cnt,fp);
	}
	H->filler_mask = (ENC_CLASS *)_mm_malloc( H->mcol*sizeof(ENC_CLASS),sizeof(ENC_CLASS) );
	punc_fill_cnt = enc_fill_mask( H->filler_cnt, H->fill_list, H->filler_mask, H->z_value, H->mcol );

    fclose(fp);
	printf("Matrix Params: Z=%d  Hcol=%d  Hrow=%d num_par_cw=%d fill_cnt=%d\n",H->z_value,H->mcol,H->mrow,H->num_par_cw,punc_fill_cnt);

	//allocate memory
	H->cw_bits   = (ENC_CLASS *)_mm_malloc( num_cw*H->mcol*sizeof(ENC_CLASS),sizeof(ENC_CLASS) );
    H->work_area = (ENC_CLASS *)_mm_malloc( (H->mcol*8*2)*sizeof(ENC_CLASS), sizeof(ENC_CLASS) );
}
