#include "vectorclass.h"

//use simd_defs to define the Encoder class
//enc_defs can use 256b or 512b vectors

#include "enc_defs.h"
#include "enc_struct_LDPC.h"

Vec4uq rand_xoshiro_N(int rand_sel);

void gen_rand_databits( unsigned num_codeword, struct Henc_struct *Henc ) {
	//assumption is that there is one SIMD word per major col (ie Z <= ENC_BW)
	//Need to fill just the data pieces with random bits, but zero out others
	ENC_CLASS z_bitmask;
	//code for z_bitmask duplicated in LDPCencode
	uint64_t *zmask_ptr = (uint64_t *)&z_bitmask;
	int z_val_copy = Henc->z_value;
	for (int i=0; i<(ENC_BW/64); ++i) {
		int shift_val = (z_val_copy >= 64) ? 0 : (64 - z_val_copy);
		zmask_ptr[i]  = (z_val_copy ==  0) ? 0 : 0xffff'ffff'ffff'ffff >> shift_val;
		z_val_copy    = (z_val_copy >= 64) ? z_val_copy - 64 : 0;
	}
	ENC_CLASS *databuff = Henc->cw_bits;
	for (unsigned cw=0; cw < num_codeword; ++cw) {
		ENC_CLASS *filler_mask = Henc->filler_mask;
		for (unsigned dcol=0; dcol < (Henc->mcol-Henc->mrow); ++dcol) {
#if ENC_BW == 256
			ENC_CLASS rand_dbits = rand_xoshiro_N(1);
#else
			Vec4uq rand_dbits_256a = rand_xoshiro_N(1);
			Vec4uq rand_dbits_256b = rand_xoshiro_N(1);
			ENC_CLASS rand_dbits(rand_dbits_256a,rand_dbits_256b);
#endif
			rand_dbits &= z_bitmask;
			rand_dbits &= ~(*filler_mask++);
			*databuff++ = rand_dbits;
		}
		databuff += Henc->mrow;
	}
}
