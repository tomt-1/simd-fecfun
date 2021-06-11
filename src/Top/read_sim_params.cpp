#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <immintrin.h>

#include "sim_struct_LDPC.h"

void read_sim_params( struct sim_struct *Sim ) {
	//read simulation parameters from stdin (pipe from file)
	//using fread to directly read items
	int fread_return;
	fread_return  = fread(&(Sim->outer_count),sizeof(unsigned),1,stdin);
	fread_return += fread(&(Sim->loop_count),sizeof(unsigned),1,stdin);
	fread_return += fread(&(Sim->num_codeword),sizeof(unsigned),1,stdin);
	fread_return += fread(&(Sim->max_iter),sizeof(unsigned),1,stdin);
	fread_return += fread(&(Sim->puncture_8023ca_flag),sizeof(int),1,stdin);
	fread_return += fread(&(Sim->use_given_rand_state),sizeof(int),1,stdin);
	fread_return += fread(&(Sim->max_row_metric),sizeof(float),1,stdin);
	fread_return += fread(&(Sim->LLR_scaling_factor),sizeof(float),1,stdin);
	fread_return += fread(&(Sim->max_abs_LLR),sizeof(float),1,stdin);
	fread_return += fread(&(Sim->beta_offset),sizeof(float),1,stdin);
	fread_return += fread(&(Sim->SNR_start),sizeof(float),1,stdin);
	fread_return += fread(&(Sim->SNR_end),sizeof(float),1,stdin);
	fread_return += fread(&(Sim->SNR_step),sizeof(float),1,stdin);
	fread_return += fread(&(Sim->rand_state),sizeof(uint64_t),4*4*4,stdin);
	fread_return += fread(&(Sim->enc_raw_file),sizeof(char),128,stdin);
	fread_return += fread(&(Sim->dec_raw_file),sizeof(char),128,stdin);
	if ( fread_return != 4+2+7+4*4*4+128+128 ) {
		printf("Error in reading sim parameters from stdin!\n");
		exit(1);
	}
	printf("CW Count: %d codewords in inner loop * num_codeword\n",Sim->loop_count * Sim->num_codeword);
	if (Sim->use_given_rand_state != 1) {
		//generate the initial rand state (and print them for potential repeated use)
		for (int i=0; i < 4*4*4; ++i) { //4 quads per SIMD word, 4*4 rand state words
			int retval = _rdrand64_step( (long long unsigned int*)(Sim->rand_state + i) );
			printf("rand seed: %016lx\n",Sim->rand_state[i]);
			if ( retval != 1 ) {
				printf("ERROR: rdrand retval=%d\n",retval);
				exit(1);
			}
		}
	}
}
