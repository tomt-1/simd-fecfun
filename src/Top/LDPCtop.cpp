#include <stdio.h>

#include "vectorclass.h"

#define DEBUG 0
#define ENC_CLASS Vec8uq

#include "../../src/LDPCEnc/enc_struct_LDPC.h"
#include "simd_defs.h"
#include "../../src/LDPCDec/dec_struct_LDPC.h"
#include "sim_struct_LDPC.h"

#define MAX_SNR_ITER 20

unsigned low_tsc0,high_tsc0,low_tsc1,high_tsc1;
#include "timing_macro.h"

#define SIMD_F Vec8f
#define SIMD_F_PAR 8

typedef void (*min_fn_ptr_type)(SIMD_CLASS *,SIMD_CLASS *,SIMD_CLASS,int);
extern min_fn_ptr_type *all_minprod_func_ary;
min_fn_ptr_type *min_funcs = all_minprod_func_ary;

//function prototypes
void read_sim_params( struct sim_struct * );
void init_randstate(uint64_t *inp);
void encode_matrix_allocate( const char *, struct Henc_struct*, unsigned );
void decode_matrix_allocate(const char *filename, struct Hdec_struct *H, unsigned num_cw);
void normdist2( SIMD_F *result, int simd_count );
void gen_rand_databits( unsigned num_codeword, struct Henc_struct *Henc );
void LDPCencode( const int num_codewords, struct Henc_struct H );
unsigned gen_metrics(unsigned num_codeword, float SNR, unsigned char *enc_cw_byteptr, int z_value, int mcol,
    SIMD_F *wgn_fltptr, SIMD_CLASS *simd_metric_ptr, float LLR_scaling_factor, float max_abs_LLR,
	ENC_CLASS *puncture_mask, ENC_CLASS *filler_mask);
unsigned LDPCdecode_set( unsigned num_codeword, unsigned max_iter, float max_row_metric, 
	float beta_offset, struct Hdec_struct *Hdec, min_fn_ptr_type *min_funcs );
void LDPC_count_errors( unsigned num_codeword, struct Henc_struct *Henc, struct Hdec_struct *Hdec,
    unsigned *tot_cw_err, unsigned *tot_bit_err );

int main() {
	//get top-level parameter information 
		//number of different Encode/Decode matrices?  For now just one
		//Encode/Decode raw matrix info filenames
		//SNR min/max/step. number of outer_loop iterations.  number of codewords (inner loop)
		//  max iter, max_local_metric, scaling/sat values for a priori metrics
		//init seed values  

	
	struct Henc_struct Henc;
	struct Hdec_struct Hdec;
	struct sim_struct Sim;		//simulation parameters for the current run

	read_sim_params( &Sim );
	unsigned num_codeword = Sim.num_codeword;

	encode_matrix_allocate( Sim.enc_raw_file, &Henc, num_codeword );
	decode_matrix_allocate( Sim.dec_raw_file, &Hdec, num_codeword );
	//allocate memory for WGN values
	int simd_per_floatz = (Henc.z_value + SIMD_F_PAR -1) / SIMD_F_PAR;
	int wgn_count = num_codeword * Henc.mcol * simd_per_floatz;
	if ( wgn_count & 1 ) {
		++wgn_count;	//make wgn_count even
	}
	SIMD_F *wgn_buff  = (SIMD_F *)_mm_malloc( wgn_count * sizeof(SIMD_F), sizeof(SIMD_F) );
	float *wgn_fltptr = (float *)wgn_buff;

	unsigned tot_bit_err_ary[MAX_SNR_ITER] = {0};
	unsigned tot_cw_err_ary[MAX_SNR_ITER] = {0};
	unsigned tot_cw_ary[MAX_SNR_ITER] = {0};
	unsigned tot_inp_err_ary[MAX_SNR_ITER] = {0};
	unsigned tot_iter[MAX_SNR_ITER] = {0};
	uint64_t tot_dec_ticks[MAX_SNR_ITER] = {0};
	//DEBUG
	SIMD_CLASS start_simd_metric[(64+2)*69+2];
	unsigned SNRiter_start = 0;
	unsigned SNRiter_end = ( Sim.SNR_end - Sim.SNR_start + 1e-6 ) / Sim.SNR_step;	//1-e6 to avoid fp rounding errors

	unsigned char *enc_cw_byteptr = (unsigned char *)Henc.cw_bits;
	SIMD_CLASS *simd_metric_ptr = Hdec.decoded_metrics+1;

	//initialize the random state
	init_randstate( Sim.rand_state );

	//main loops
	for (unsigned outloop=0; outloop < Sim.outer_count; ++outloop) {
		unsigned uncoded_err_cnt = 0;
		for (unsigned inloop=0; inloop < Sim.loop_count; ++inloop) {
			//get WGN values
			normdist2( wgn_buff, wgn_count );
			gen_rand_databits( num_codeword, &Henc );
			LDPCencode( num_codeword, Henc );
			for (unsigned SNRiter = SNRiter_start; SNRiter <= SNRiter_end; ++SNRiter) {
				unsigned tot_cw_err;
				unsigned tot_bit_err;
				float SNR = SNRiter * Sim.SNR_step + Sim.SNR_start;
				uncoded_err_cnt = gen_metrics( num_codeword, SNR, enc_cw_byteptr, Henc.z_value, Henc.mcol,
					wgn_buff, simd_metric_ptr, Sim.LLR_scaling_factor, Sim.max_abs_LLR, Henc.puncture_mask, Henc.filler_mask);
				//DEBUG: copy all starting metrics in first word to "start_simd_metric"
				if ( DEBUG ) {
					for (unsigned i=0; i < (64+2)*Henc.mcol; ++i) {
						start_simd_metric[i] = simd_metric_ptr[i];
					}
				}
				tot_inp_err_ary[SNRiter] += uncoded_err_cnt;
				//decode codewords (and save statistics - maybe option of brief or detailed)
				RDTSC_RAW(low_tsc0,high_tsc0)
				unsigned iteration_cnt;
				iteration_cnt = LDPCdecode_set( num_codeword, Sim.max_iter, Sim.max_row_metric, Sim.beta_offset, &Hdec, min_funcs );
				RDTSC_RAW(low_tsc1,high_tsc1)
				tot_iter[SNRiter] += iteration_cnt;
				uint64_t tsc_ticks = (((uint64_t)high_tsc1  << 32) | low_tsc1 ) -
                        (((uint64_t)high_tsc0 << 32) | low_tsc0 );
				tot_dec_ticks[SNRiter] += tsc_ticks;

				//evaluate decoded word (BER etc)
				LDPC_count_errors( num_codeword, &Henc, &Hdec, &tot_cw_err, &tot_bit_err );

				tot_bit_err_ary[SNRiter] += tot_bit_err;
				tot_cw_err_ary[SNRiter] += tot_cw_err;
				tot_cw_ary[SNRiter] += num_codeword;
			}
		}
		//report statistics (clear counts)
		for (unsigned SNRiter = SNRiter_start; SNRiter <= SNRiter_end; ++SNRiter) {
			float SNR = SNRiter * Sim.SNR_step + Sim.SNR_start;
			printf("loop %d: SNR=%f  inp_bit_err=%d tot_iter=%d tot_bit_err=%d tot_cw_err=%d tot_dec_ticks=%ld\n",
				outloop, SNR, tot_inp_err_ary[SNRiter], tot_iter[SNRiter], tot_bit_err_ary[SNRiter], 
				tot_cw_err_ary[SNRiter], tot_dec_ticks[SNRiter]);
			fflush(stdout);
			tot_inp_err_ary[SNRiter] = 0;
			tot_iter[SNRiter] = 0;
			tot_bit_err_ary[SNRiter] = 0;
			tot_cw_err_ary[SNRiter] = 0;
			tot_dec_ticks[SNRiter] = 0;
		}
		//adjust SNR_start based on parameters (number of total cw seen with error exceed a threshold?)
		//    might also do one-time adjustment of SNR_end based on keeping N SNR buckets that had 0 errors
	}
	if (DEBUG) {
	//print data for first cw of last loop
	//print cw bits (in groups of bytes)
	unsigned char *e_ptr = (unsigned char *)Henc.cw_bits;
	for (unsigned col=0; col < Henc.mcol; ++col) {
		for (unsigned byte=0; byte < 32; ++byte) {
			if ( byte < (Henc.z_value+7)/8 ) {
				printf("%02x ",*e_ptr);
			}
			++e_ptr;
		}
		printf("\n");
	}
	printf("end of cw bits\n");
	//print metrics
	float *m_ptr = wgn_fltptr;
	for (unsigned col=0; col < Henc.mcol; ++col) {
		for (unsigned i=0; i < (8*((Henc.z_value+7)>>3)); ++i) {
			printf("%6.4f ",*m_ptr++);
		}
		printf("\n");
	}
	printf("end of WGN data\n");
	//double *d1_ptr = (double *)start_simd_metric;
	int32_t *d1_ptr = (int32_t *)start_simd_metric;
	for (unsigned col=0; col < Henc.mcol; ++col) {
		for (unsigned i=0; i < Henc.z_value; ++i) {
			//printf("%4.12f ",*d1_ptr++);
			printf("%d ",*d1_ptr++);
		}
		printf("\n");
		d1_ptr += 2*8+6;	//2 SIMD words of 8 ints per word, plus 6 if z_value is 42
		//d1_ptr += 2*4; //2 SIMD words, 4 doubles per SIMD word
	}
	printf("end of start metrics\n");
	//double *d_ptr = (double *)simd_metric_ptr;
	int32_t *d_ptr= (int32_t *)simd_metric_ptr;
	for (unsigned col=0; col < Henc.mcol; ++col) {
		for (unsigned i=0; i < Henc.z_value; ++i) {
			//printf("%7.5f ",*d_ptr++);
			printf("%d ",*d_ptr++);
		}
		printf("\n");
		d_ptr += 2*8+6;	//Z=42, 32b elem
		//d_ptr += 2*4; //2 SIMD words, 4 doubles per SIMD word
	}
	printf("end of dec metrics\n");
	}
}
