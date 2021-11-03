#include <stdio.h>
#include <string.h>
#include "vectorclass.h"

#include "simd_defs.h"

#include "../../src/LDPCDec/dec_struct_LDPC.h"

#include "enc_defs.h"
#include "../../src/LDPCEnc/enc_struct_LDPC.h"


#if ( SIMD_PAR <= 8 )
	#define DB_TYPE uint8_t
#elif ( SIMD_PAR == 16 )
	#define DB_TYPE uint16_t
#elif ( SIMD_PAR == 32 )
	#define DB_TYPE uint32_t
#else
	//assume that decoder must be 8-bit metrics and 512-bit SIMD word
	#define DB_TYPE uint64_t
#endif

//maximum num_par_cw allowed
#define MAX_PAR_CW 16

//1s every Nth bit, for N=1 to 16.  Where N is number of codewords in parallel in H matrix
static const uint64_t par_cw_mask[MAX_PAR_CW+1] = {
	0x0, 
	0xffff'ffff'ffff'ffff, 0x5555'5555'5555'5555,
	0x9249'2492'4924'9249, 0x1111'1111'1111'1111,
	0x1084'2108'4210'8421, 0x1041'0410'4104'1041,
	0x8102'0408'1020'4081, 0x0101'0101'0101'0101,
	0x8040'2010'0804'0201, 0x1004'0100'4010'0401,
	0x0080'1002'0040'0801, 0x1001'0010'0100'1001,
	0x0010'0080'0400'2001, 0x0100'0400'1000'4001,
	0x1000'2000'4000'8001, 0x0001'0001'0001'0001
	};

void LDPC_count_errors( unsigned num_codeword, struct Henc_struct *Henc, struct Hdec_struct *Hdec,
	unsigned *tot_cw_err, unsigned *tot_bit_err ) { //for now, ignoring iter_cnts: unsigned *iter_cnts ) {

	unsigned simd_per_z = (Hdec->z_value+SIMD_PAR-1) / SIMD_PAR;
	SIMD_CLASS *dec_ptr = Hdec->decoded_metrics + 1;
	ENC_CLASS *enc_data = Henc->cw_bits;
	ENC_CLASS z_bitmask;	//identical to z_bitmask calculated in LDPCencode - perhaps add to Henc
	DB_TYPE *dec_local_ptr;
	ENC_CLASS dec_local[2]; //adding an extra word in case any write can go past 1st (might remove)
	unsigned cw_with_error = 0;
	unsigned local_tot_bit_errors = 0;
	unsigned num_par_cw = Henc->num_par_cw;

	//assign z_bitmask
	//due to pointer aliasing, memcopy z_mask_ptr results to z_bitmask.
	uint64_t zmask_ptr[ENC_BW/64];
	//uint64_t *zmask_ptr = (uint64_t *)&z_bitmask;
	int z_val_copy = Henc->z_value;
	for (int i=0; i<(ENC_BW/64); ++i) {	
		int shift_val = (z_val_copy >= 64) ? 0 : (64 - z_val_copy);
		zmask_ptr[i]  = (z_val_copy ==  0) ? 0 : 0xffff'ffff'ffff'ffff >> shift_val;
		z_val_copy    = (z_val_copy >= 64) ? z_val_copy - 64 : 0;
	}
	memcpy(&z_bitmask,zmask_ptr,(ENC_BW/64)*sizeof(uint64_t));
	const uint64_t par_mask = par_cw_mask[num_par_cw];

	//get error counts over all codwords
	for (unsigned cw=0; cw < num_codeword; ++cw) {
		unsigned cw_bit_errors = 0;
		unsigned cw_err_cnt_by_slot[MAX_PAR_CW] = {0};
		
		//compress signs of decoded codeword data (not parity) to bits
		DB_TYPE dec_local_store[(ENC_BW/8)*2]; //room for 2 256b or 512b words if 8b metrics.  More than enough for others
		for (unsigned col=0; col < (Hdec->mcol - Hdec->mrow); ++col) {
			dec_local_ptr = dec_local_store;
			//dec_local_ptr = (DB_TYPE *)dec_local;
			unsigned cw_col_errors = 0;
			// SIMD_PAR == 4 if double prec and 256b => only 4 decoded bits per vector
			#if ( SIMD_PAR == 4 )
				DB_TYPE prev_dec_bits = 0;
			#endif
			for (unsigned i=0; i < simd_per_z; ++i) {
				SIMD_CLASS dec_word = *dec_ptr++;
				//if SIMD_CLASS is a float-type, assume no -0s
				DB_TYPE dec_bits = to_bits( dec_word < 0 );
				#if ( SIMD_PAR == 4 )
					if ( i & 1 ) {
						dec_bits <<= 4;
						dec_bits |= prev_dec_bits;
					}
					if ( (i & 1) || (i == (simd_per_z-1)) ) {
						*dec_local_ptr++ = dec_bits;
					}
					prev_dec_bits = dec_bits;
				#else
					*dec_local_ptr++ = dec_bits;
				#endif
			}
			dec_ptr += 2;	//skip the extra SIMD words in each major col
			//compare to original codeword data (assumes Z <= 256)
			//avoid pointer aliasing
			memcpy(dec_local,dec_local_store,sizeof(ENC_CLASS));
			ENC_CLASS diff_tmp1 = *enc_data ^ dec_local[0];
			ENC_CLASS diff = diff_tmp1 & z_bitmask;
			++enc_data;
			uint64_t diff_64ptr[ENC_BW/64];
			//use memcpy to avoid illegal pointer aliasing.  Should optimize away
			memcpy(diff_64ptr,&diff,(ENC_BW/64)*sizeof(uint64_t));
			for (unsigned i=0; i<(ENC_BW/64); ++i) {
				uint64_t mask_shift = par_mask;
				cw_col_errors += _popcnt64( diff_64ptr[i] );
				for (unsigned par=0; par < num_par_cw; ++par) {
					uint64_t diff_par = diff_64ptr[i] & mask_shift;
					unsigned err_cnt = _popcnt64( diff_par );
					unsigned cw_slot = (i*64+par) % num_par_cw;
					cw_err_cnt_by_slot[cw_slot] += err_cnt;
					mask_shift = (mask_shift << 1);
				}
			}
			cw_bit_errors += cw_col_errors;
		}
		//advance past parity bits 
		enc_data += Hdec->mrow;
		dec_ptr  += Hdec->mrow * (simd_per_z+2);
		if ( cw_bit_errors != 0 ) {
			//DEBUG
			//printf("  cw error at cw number %d, bit_errors=%d\n",cw,cw_bit_errors);
			for (unsigned par=0; par < num_par_cw; ++par) {
				cw_with_error += (cw_err_cnt_by_slot[par] > 0);
			}
			//++cw_with_error;
			local_tot_bit_errors += cw_bit_errors;
		}
	}
	*tot_cw_err = cw_with_error;
	*tot_bit_err = local_tot_bit_errors;
}
