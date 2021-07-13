#include "vectorclass.h"

//use simd_defs to define the Encoder class
//vec4uq is only one supported ATM.  Expect to add vec8uq

//#include "simd_defs.h"
//#define ENC_CLASS SIMD_CLASS
#define ENC_CLASS Vec4uq
#include "enc_struct_LDPC.h"

ENC_CLASS rand_xoshiro_N(int rand_sel);

void gen_rand_databits( unsigned num_codeword, struct Henc_struct *Henc ) {
	//assumption is that there is one SIMD word per major col (ie Z <= 256)
	//Need to fill just the data pieces with random bits, but zero out others
	ENC_CLASS z_bitmask;
	//code for z_bitmask duplicated in LDPCencode
	uint64_t *zmask_ptr = (uint64_t *)&z_bitmask;
	int z_val_copy = Henc->z_value;
	for (int i=0; i<4; ++i) {  //assumes 256b SIMD words 4*64=256
		int shift_val = (z_val_copy >= 64) ? 0 : (64 - z_val_copy);
		zmask_ptr[i]  = (z_val_copy ==  0) ? 0 : 0xffff'ffff'ffff'ffff >> shift_val;
		z_val_copy    = (z_val_copy >= 64) ? z_val_copy - 64 : 0;
	}
	ENC_CLASS *databuff = Henc->cw_bits;
	for (unsigned cw=0; cw < num_codeword; ++cw) {
		for (unsigned dcol=0; dcol < (Henc->mcol-Henc->mrow); ++dcol) {
			ENC_CLASS rand_dbits = rand_xoshiro_N(1);
			rand_dbits &= z_bitmask;
			*databuff++ = rand_dbits;
		}
		databuff += Henc->mrow;
	}
}
