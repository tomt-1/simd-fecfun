#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "vectorclass.h"

#ifndef ENC_CLASS
#define ENC_CLASS Vec4uq
#endif

#include "enc_struct_LDPC.h"

//load H matrix items for encode

void encode_matrix_allocate (const char *filename, struct Henc_struct *H, unsigned num_cw) {

	FILE *fp;
	int fread_return;

	//Consider adding error checking on opening file and all fread return values

	//read the Encode Matrix data
	fp = fopen(filename,"r");
    fread_return  = fread(&(H->z_value),sizeof(int32_t),1,fp);
	fread_return += fread(&(H->num_par_cw),sizeof(int32_t),1,fp);
    fread_return += fread(&(H->mrow),sizeof(int32_t),1,fp);
    fread_return += fread(&(H->mcol),sizeof(int32_t),1,fp);
    if ( fread_return != 4 ) {
        printf("Error in fread mrow or mcol\n");
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
    fclose(fp);
	printf("Matrix Params: Z=%d  Hcol=%d  Hrow=%d num_par_cw=%d\n",H->z_value,H->mcol,H->mrow,H->num_par_cw);

	//allocate memory
	H->cw_bits   = (ENC_CLASS *)_mm_malloc( num_cw*H->mcol*sizeof(ENC_CLASS),sizeof(ENC_CLASS) );
    H->work_area = (ENC_CLASS *)_mm_malloc( (H->mcol*8*2)*sizeof(ENC_CLASS), sizeof(ENC_CLASS) );
}
