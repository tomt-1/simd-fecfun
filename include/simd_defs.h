//ELEM_NUM: 8=>8b int, 16=>16b int, 32=>32b int, -32=>float, -64=>double
//define SIMD_PAR and SIMD_ELEM_BW, all others are defined from these

#ifndef SIMD_PAR
	#define SIMD_PAR -64
#endif
#ifndef SIMD_ELEM_BW
	#define SIMD_ELEM_BW 4
#endif

#ifndef USE_BETA_OFFSET
	#define USE_BETA_OFFSET 1
#endif

#if SIMD_ELEM_BW < 0
	#define FLOAT_NOT_INT 1
	#define USE_SAT_ADD_SUB 0
#else
	#define FLOAT_NOT_INT 0
	#define USE_SAT_ADD_SUB 1
#endif

#if (SIMD_PAR*SIMD_ELEM_BW == 256) || (SIMD_PAR*SIMD_ELEM_BW == -256)
	#define WIDTH 256
	#define MTYPE uint32_t
#else
	#define WIDTH 512
	#define MTYPE uint64_t
#endif

#define THREE_CAT1(x,y,z) x ## y ## z
#define THREE_CAT(x,y,z) THREE_CAT1(x,y,z)

#if SIMD_ELEM_BW == 8
	#define SIMD_CLASS THREE_CAT(Vec,SIMD_PAR,c)
	#define CTYPE int8_t
	#define SIMD_INT SIMD_CLASS
#elif SIMD_ELEM_BW == 16
	#define SIMD_CLASS THREE_CAT(Vec,SIMD_PAR,s)
	#define CTYPE int16_t
	#define SIMD_INT SIMD_CLASS
#elif SIMD_ELEM_BW == 32
	#define SIMD_CLASS THREE_CAT(Vec,SIMD_PAR,i)
	#define CTYPE int32_t
	#define SIMD_INT SIMD_CLASS
#elif SIMD_ELEM_BW == -32
	#define SIMD_CLASS THREE_CAT(Vec,SIMD_PAR,f)
	#define CTYPE float
	#define SIMD_INT THREE_CAT(Vec,SIMD_PAR,i)
#elif SIMD_ELEM_BW == -64
	#define SIMD_CLASS THREE_CAT(Vec,SIMD_PAR,d)
	#define CTYPE double
	#define SIMD_INT THREE_CAT(Vec,SIMD_PAR,q)
#endif

//testcpp: sc=SIMD_CLASS sp=SIMD_PAR flag=FLOAT_NOT_INT w=WIDTH int=SIMD_INT ctype=CTYPE
