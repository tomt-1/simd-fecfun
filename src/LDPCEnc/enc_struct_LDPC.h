//Hencoder structure
//major rows/cols is total number of rows/cols divided by Z
//                   also number of rows/cols in compressed Hc matrix
//There are P major columns of parity bits (so P*Z total parity bits)
//par_cnts points to array of P integers with number of offsets to compute the current parity column
//par_colnum points to array of P integers in range of 0 to P-1 with the parity column number being computed
//par_offset points to array of offset in memory to use for next XOR operation (in computing parity)

//note that the offsets assume 2Z-sized copies of input data at all 8 bit offsets are available
//see LDPCencode.cpp for exact creation of the array that par_offset refers to

struct Henc_struct {
	uint32_t z_value;		//Z value for the pseudo-cyclic H matrix
	uint32_t num_par_cw;    //number of interleaved sub-codewords in the defined H matrix
	uint32_t mrow;			//Num of major rows in H matrix (rows in compressed Hc matrix)
	uint32_t mcol;			//Num of major columns in H matrix/cols in Hc matrix
	uint32_t *par_cnts;		//pointer to count-of-offsets for each parity column
	uint32_t *par_colnum;	//pointer to parity column number being computed
	uint32_t *par_offset;	//pointer to offsets used in computing parity
	ENC_CLASS *cw_bits;		//pointer to memory for all codeword bits (data+parity)
	ENC_CLASS *work_area;	//pointer to work area used by LDPC encoder
};
