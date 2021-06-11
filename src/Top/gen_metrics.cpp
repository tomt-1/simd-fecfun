#include <math.h>
#include "vectorclass.h"

#include "simd_defs.h"

#define SIMD_F Vec8f

//Because the Encode and WGN buffers are traversed in a manner that
//doesn't depend on the underlying SIMD types, simd_defs.h is 
//used here for the Decoder Metrics definitions

//Since this isn't a major cycle-time contributor, encoded bits are
//accessed a byte at a time, and WGN floats are accessed 8 at a time
//Decoder Metrics could be one of several SIMD classes

//return uncoded_err_cnt (might want to make it 64b)
unsigned gen_metrics(unsigned num_codeword, float SNR, unsigned char *enc_cw_byteptr, int z_value, int mcol, 
	SIMD_F *wgn_fltptr, SIMD_CLASS *simd_metric_ptr, float LLR_scaling_factor, float max_abs_LLR, int puncture_8023ca_flag) {

	//expand encoded bits to 8 float values of +1.0 (enc bit was 0) or -1.0 (enc bit was 1)
	//scale WGN and add to this.  Get a count of sign changes (these are uncoded bit errors due to noise)
	const Vec8i bit_select_mask(0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80);

	//SNRdB = -20*log10(sigma) for real-valued noise, subtract 10*log10(2) for complex-valued noise
	//sigma = 10^( SNRdb/-20 ) for pure real-valued noise, or 10^( (SNR+10*log10(2))/-20 ) for complex.  (BPSK)
	float Three_dB = 10*log10f(2.0);
	float sigma_scalar = powf( 10.0, ((SNR+Three_dB)/-20.0) );
	float sigma_factor_scalar = 2*LLR_scaling_factor/(sigma_scalar*sigma_scalar);
	SIMD_F sigma = sigma_scalar;
	SIMD_F sigma_factor = sigma_factor_scalar;

	int last_byte_iter = (z_value+7)/8;
	int last_count = ((z_value-1) & 7) + 1;
	int enc_col_offset = (256/8) - last_byte_iter;	//encode is always 256b SIMD.  Offset to next major col
	int dec_simd_per_z = (z_value+SIMD_PAR-1) / SIMD_PAR;  //ceiling function
	int dec_col_offset = SIMD_PAR * dec_simd_per_z - 8*last_byte_iter;  //par*simd_per_z = total elements per major column stored
	dec_col_offset += 2*SIMD_PAR; //and skip 2 additional SIMD words due to metric layout (empty SIMD word before and after each mcol)

	CTYPE *metric_ptr = (CTYPE *)simd_metric_ptr;	//TBD:make sure there is no need to add 1 (input starts at actual metric, not blank)

	Vec8ui err_sum = 0;

	for (unsigned cw=0; cw < num_codeword; ++cw) {
		for (int col=0; col < mcol; ++col) {
			int no_puncture_flag = ( puncture_8023ca_flag && (col >= (69-2)) ) - 1;
			Vec8i nopunc_int = no_puncture_flag;
			Vec8f nopunc_mask = reinterpret_f(nopunc_int);
			for (int byte_cnt=0; byte_cnt < last_byte_iter; ++byte_cnt) {
				uint8_t byte_val = *enc_cw_byteptr++;
				Vec8i PosOne =  1;
				Vec8i NegOne = -1;
				Vec8i bit_expand = byte_val;
				bit_expand &= bit_select_mask;
				bit_expand = select( bit_expand == 0, PosOne, NegOne );
				Vec8f LLR_orig = to_float( bit_expand );
				Vec8f Rx_sig = mul_add( *wgn_fltptr++, sigma, LLR_orig );
				//if last iteration, zero the unused elements of Rx_sig
				if ( byte_cnt == (last_byte_iter-1) ) {
					Rx_sig.cutoff( last_count );
					//no need to cutoff LLR_orig (i.e. byte_val is zero beyond z_value)
				}
				//NOTE: err_cnt based on floating pt addition of noise.  If converted to integer metrics
				//some individual "0" metrics will be different sign.  Same total input BER though.
				//check on bit errors
				Vec8ui chk_orig = reinterpret_i( LLR_orig );
				Vec8ui chk_Rx   = reinterpret_i( Rx_sig );
				Vec8ui err_cnt  = (chk_orig ^ chk_Rx) >> 31;	//signs differ => mismatch.  Shift makes integer value 1 if mismatch
				err_sum += err_cnt;
				//Consider adding sum for this codeword.  If no bit errors after encoding, skip decoding
				//(would only apply to codeword with punctured bits, of course.  And very high SNR)

				Vec8f LLR = Rx_sig * sigma_factor;	//this is the true LLR if LLR_scaling_factor is 1.0
				LLR = min( max_abs_LLR, LLR );
				LLR = max(-max_abs_LLR, LLR );
				LLR = LLR & nopunc_mask;

				//convert to 8 of whatever the Decode Metric type is, and store it
				#if ( SIMD_ELEM_BW == -64 )
					Vec8d Metric = to_double( LLR );
					Metric.store( metric_ptr );	//note that this can go past end of normal metric buffer
				#elif ( SIMD_ELEM_BW == -32 )
					Vec8f Metric = LLR;
					Metric.store( metric_ptr );
				#else
					Vec8i Metrics_int = roundi(LLR);
					#if ( SIMD_ELEM_BW == 16 )
						Vec8s Metric = compress( Metrics_int );
						Metric.store( metric_ptr );
					#elif ( SIMD_ELEM_BW == 8 )
						Vec8s Metric_s = compress( Metrics_int );
						Vec16c Metric_c = (Vec16c)Metric_s;
						Metric_c = permute16<0,2,4,6,8,10,12,14, -1,-1,-1,-1, -1,-1,-1,-1>(Metric_c);
						Metric_c.store_partial( 8, metric_ptr );
					#else //must be 32b integer
						Metrics_int.store( metric_ptr );
					#endif
				#endif
				metric_ptr += 8;
			} //end of byte_cnt in Z loop
			//add to metric_ptr, enc_cw_byteptr, and wgn_buff to get to next major column
			enc_cw_byteptr += enc_col_offset;
			metric_ptr += dec_col_offset;
		} //end of col loop
		//shouldn't need any adjustment to pointer here -- should already be at next codeword (?)
	} //end of codeword loop
	unsigned err_total = horizontal_add( err_sum );
	return err_total;
}
