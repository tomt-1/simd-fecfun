#include "vectorclass.h"
#include "vectorclass_extensions.h"

//ENC_CLASS uses 256b vectors.  Might add support for 512b (Vec8uq) in the future
#include "enc_defs.h"

#include "enc_struct_LDPC.h"

//Restriction: Z must be less than or equal to 256
//should be straightforward to add support for 512b vectors
//support for Z values > vector size is harder -- substantial re-work

//Various values derived from Z (the overall Z, including any parallel multiplier)
//for completeness: parity matrix H is divided into ZxZ subarrays.
struct z_param_st {
	ENC_CLASS z_bitmask;
    unsigned zmod8;
    unsigned zmod8_rev;
	int byte_offset;
	int byte_offset_rev;
};

//This function creates 2Z-sized copies of input data at all 8 bit offsets
//This allows directly using (unaligned) SIMD words in XOR computations 
//for AVX512VBMI, it is likely faster to use permutexvar_epi8. Otherwise, the 128b lanes lead
//   to "lookup" function requiring several instructions
void copy_rotations( ENC_CLASS inp_data, struct z_param_st z, ENC_CLASS *prework, ENC_CLASS *dest );

//Prework offsets hold codeword data in SIMD words at follows:
//01: cw data 
//03: cw shift left 1
//05: cw shift left 2, ... 
//15: cw shift left 7
//16: cw data
//18: cw shift right 1 ...
//30: cw shift right 7
//other offsets are zero.  Function copy_rotations fills in the above shifted values
static ENC_CLASS Prework_area[32] = {0};

void LDPCencode( const int num_codewords, struct Henc_struct H )
{
	//num_codewords is number of codewords to encode
	//H is the Encoder structure, including codeword data bits
	//note that the codeword parity bits for each codeword are placed after the data bits,
	//to form a complete codeword
	//also, bits are grouped in one major column per 256b word
	//
	//work_area has 8bit-rotations of 2x codewords for each column
	//  = majcol*8*2 codewords

	//compute the Z-based items, including mask
	struct z_param_st z;
	z.zmod8 = H.z_value & 7;
	z.zmod8_rev = (8-z.zmod8) & 7;
	z.byte_offset = (H.z_value >> 3);
	z.byte_offset_rev = (ENC_BW-H.z_value) >> 3;

	uint64_t *zmask_ptr = (uint64_t *)&z.z_bitmask;
	int z_val_copy = H.z_value;
	for (int i=0; i < (ENC_BW/64); ++i) {  //Num of 64b words: 4 if 256b ENC_BW, 8 if 256b ENC_BW
		int shift_val = (z_val_copy >= 64) ? 0 : (64 - z_val_copy);
		zmask_ptr[i]  = (z_val_copy ==  0) ? 0 : 0xffff'ffff'ffff'ffff >> shift_val;
		z_val_copy    = (z_val_copy >= 64) ? z_val_copy - 64 : 0;
	}

	//define pointers, number of data (major) columns
	ENC_CLASS *rotation_data = H.work_area;
	ENC_CLASS *prework  = Prework_area;
	const int tot_data_col = H.mcol - H.mrow;

	ENC_CLASS *cw_data = H.cw_bits;
	ENC_CLASS *cw_parity = H.cw_bits + tot_data_col;
	for (int cw=0; cw < num_codewords; ++cw) {
		//for each cw, copy shifts of cw data to work area
		//(NOTE: probably faster execution if parity columns were built one column at a time)
		for (int col=0; col < tot_data_col; ++col) {
			//get all 8 bit-offsets of back-to-back (2Z) copies of each mcol
			copy_rotations( cw_data[col], z, prework, rotation_data+16*col );
		}
		//calculate all parity major columns (and rotations for parity columns)
		int8_t *wa_byte = (int8_t *)rotation_data;
		int idx = 0;
		for (unsigned par_idx=0; par_idx < H.mrow; ++par_idx) {
			ENC_CLASS xor_accum = 0;
			for ( unsigned i=0; i < H.par_cnts[par_idx]; ++i) {
				int offset = H.par_offset[idx];
				ENC_CLASS data;
				data.load( wa_byte + offset );
				xor_accum ^= data;
				++idx;
			}
			xor_accum &= z.z_bitmask;
			int par_col = H.par_colnum[par_idx];
			cw_parity[par_col] = xor_accum;
			copy_rotations( xor_accum, z, prework, rotation_data+16*(tot_data_col+par_col) );
		}
		cw_data += H.mcol;
		cw_parity += H.mcol;
	}
}

void copy_rotations( ENC_CLASS inp_data, struct z_param_st z, ENC_CLASS *prework, ENC_CLASS *dest ) {
	//calculate the left shifts and right shifts of the current column of cw data
	ENC_CLASS curr_cw_ldata = inp_data;
	ENC_CLASS curr_cw_rdata = inp_data;
	prework[1] = curr_cw_ldata;
	prework[16+0] = curr_cw_rdata;
	for (int rot=1; rot<8; ++rot) {
		curr_cw_ldata = shift_left1(curr_cw_ldata);
		curr_cw_ldata &= z.z_bitmask;
		curr_cw_rdata = shift_right1(curr_cw_rdata);
		prework[rot*2+1] = curr_cw_ldata;
		prework[16+rot*2] = curr_cw_rdata;
	}
	//calculate the combined 2Z copies at 8 single-bit offsets
	//first SIMD word is cw_sr_0..7 OR with (byte-offset)cw_sr_N
	int sl_idx = z.zmod8;
	int sl_offset = z.byte_offset;
	int sr_idx = z.zmod8_rev;
	int sr_offset = z.byte_offset_rev;
	uint8_t *sl_byteptr = (uint8_t *)(prework+1);    //shift-left byte pointer
	uint8_t *sr_byteptr = (uint8_t *)(prework+16);   //shift-right byte pointer

	for (int i=0; i<8; ++i) {
		ENC_CLASS cw1,cw2;
		cw1.load( sl_byteptr+sl_idx*2*(ENC_BW/8)-sl_offset );
		cw1 |= prework[16+i*2];
		dest[i*2] = cw1;
		--sl_idx;
		sl_offset -= ((sl_idx & 7) == 7);
		sl_idx = sl_idx & 7;

		//note that for ENC_BW==256 and Z<=128, this is 0.  (Should also be 0 for 512 and Z<=256)
		cw2.load( sr_byteptr+sr_idx*2*(ENC_BW/8)+sr_offset );
		dest[i*2+1] = cw2;
		++sr_idx;
		sr_offset += (sr_idx == 8);
		sr_idx = sr_idx & 7;
	}
}
