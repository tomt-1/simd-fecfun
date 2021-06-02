#include <stdio.h>

#include "vectorclass.h"

//ENC_CLASS is Vec4uq (support for Vec8uq will be investigated if any need for Z>256)
#define ENC_CLASS Vec4uq
#define ENC_SIZE 256

#include "enc_struct_LDPC.h"

//unsigned rdtsc_start_high, rdtsc_start_low, rdtsc_stop_high, rdtsc_stop_low;
//uint64_t rdtsc_counts[10];
//unsigned rdtsc_idx = 0;
//#include "timing_macro.h"

//function prototypes
void encode_matrix_allocate( const char *, struct Henc_struct*, unsigned );
void LDPCencode( const int num_codewords, struct Henc_struct H );

int main() {
	const unsigned num_codeword = 10; //MUST MATCH VALUE OF num_cw_per IN LDPCenc_test.m!!!
	struct Henc_struct H;
	char raw_file[] = "LDPCencode.raw.dat";

	encode_matrix_allocate( raw_file, &H, num_codeword );

	uint8_t *cw_data = (uint8_t *)H.cw_bits;
	uint8_t *cw_parity = (uint8_t *)(H.cw_bits + (H.mcol-H.mrow));

	for (unsigned i=0; i<num_codeword; ++i) {
		//read bits from stdin
		char linedata[H.z_value];
		char dummy;
		for (unsigned col=0; col<(H.mcol-H.mrow); ++col) {
			int fread_return = fread(linedata,sizeof(char),H.z_value,stdin);
			fread(&dummy,sizeof(char),1,stdin); //skip LF (probably need to skip 2 if Windoze)
			char *lineptr = linedata;
			for (unsigned bytecnt=0; bytecnt < ENC_SIZE/8; ++bytecnt) {
				cw_data[bytecnt] = 0;
				for (int j=0; j<8; ++j) {
					if ( (bytecnt*8+j) < H.z_value ) {
						cw_data[bytecnt] |= ( (*lineptr++ - '0') << j );
					}
				}
			}
			cw_data += ENC_SIZE/8;
		}
		cw_data += (ENC_SIZE/8)*H.mrow;	//increment cw_data to next codeword (skip parity)
	}

	LDPCencode( num_codeword, H );

	for (unsigned i=0; i<num_codeword; ++i) {
		//write encoded bits to stdout
		for (unsigned col=0; col < H.mrow; ++col) {
			for (unsigned bytecnt=0; bytecnt < ENC_SIZE/8; ++bytecnt) {
				for (int j=0; j<8; ++j) {
					if ( (bytecnt*8+j) < H.z_value ) {
						uint8_t temp = (cw_parity[bytecnt] >> j) & 1;
						printf("%d",temp);
					}
				}
			}
			printf("\n");
			cw_parity += ENC_SIZE/8;
		}
		cw_parity += (ENC_SIZE/8)*(H.mcol-H.mrow);
	}
}
