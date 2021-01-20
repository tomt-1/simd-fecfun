/*
Functions for generating random bit sequences and WGN

Using 256b SIMD classes as there isn't expected to be significant performance gains by explicit use of 
	512b classes where AVX-512 is available

PRNG is either xoshiro256++ or sfc64
*/

#include "vectorclass.h"

//allow up to N parallel SIMD-words of generators
//(currently, 1 for codeword generation and 1 for WGN generation)
#define PAR_GEN 2

//State elements.  Must be initialized such that each set-of-4 are not all zero
//Note that original usage was to copy a 64b initial value 4 times to fill 256b state
//Here, all 256b are initialized

static Vec4uq s[4*PAR_GEN];

void init_randstate(uint64_t *inp) {
	for (int i=0; i < 4*PAR_GEN; ++i) {
		s[i].load( inp+4*i );
	}
}

//Generate random SIMD word using Nth set of state elements 
//caller must guarantee 0 <= rand_sel < PAR_GEN
//Vec4uq => running 4 64b generators in parallel

//4 parallel PRNGs using sfc64
inline Vec4uq rand_sfc64_N(int rand_sel) {
	Vec4uq *sl = s+rand_sel*4;
	Vec4uq temp = sl[0] + sl[1] + sl[3];
	sl[3] = sl[3] + 1;
    sl[0] = sl[1] ^ (sl[1] >> 11);
    sl[1] = sl[2] + (sl[2] << 3);
    sl[2] = ((sl[2] << 24) | (sl[2] >> (64-24))) + temp;
    return temp;
}

//4 parallel PRNGs using xoshiro256++
inline Vec4uq rand_xoshiro_N(int rand_sel) {
	Vec4uq *sl = s+rand_sel*4;
	Vec4uq a = sl[0] + sl[3];
	Vec4uq b = (a << 23);
	Vec4uq c = (a >> (64-23));
	Vec4uq result = (b | c) + sl[0];

	Vec4uq t1 = (sl[1] << 17);
	sl[2] = sl[2] ^ sl[0];
	sl[3] = sl[3] ^ sl[1];
    sl[1] = sl[1] ^ sl[2];
    sl[0] = sl[0] ^ sl[3];
    sl[2] = sl[2] ^ t1;
	Vec4uq e = (sl[3] << 45);
	Vec4uq f = (sl[3] >> (64-45));
	sl[3] = e | f;
	return result;
}

//note that the output of PRNG functions is inverted to minimize "zero-land" effects
//(This puts near-zeros in the middle of Gaussian distribution instead of the edges)

//Fast approximate SIMD log function
//requires that x is positive and not subnormal ( > 1.2e-38 roughly )
//Based on approach describe by David Goldberg in Fast Approximate Logarithms
//computes -2*log(x), since this is used in Box-Muller transform
inline Vec8f fastm2log( const Vec8f x ) {
	const Vec8f coeff_a( -0.6296735f );
	const Vec8f coeff_b(  1.466967f );
	const Vec8f m2log2 ( -2.0f * 0.6931472f );	//minus 2 * log(2)

	Vec8i x_as_int = reinterpret_i( x );
	Vec8i exp = x_as_int >> 23;	//shift will remove significand

	Vec8i more1p5 = (x_as_int & 0x0040'0000) ^ 0x0040'0000;
	Vec8i div2_flag = (more1p5 >> 22);
	Vec8i expbits = (more1p5 << 1) | 0x3f00'0000;
	Vec8i frac_as_int = ( x_as_int & 0x007f'ffff ) | expbits;
	const Vec8f one_flt(1.0f);
	Vec8f frac = reinterpret_f( frac_as_int ) - one_flt;

	exp = exp - ( 126 | div2_flag );	//sub 126 or 127
	Vec8f exp_f = to_float( exp );

	//the approximate log2(x) result = exp_f + frac*(coeff_a*frac + coeff_b);
	Vec8f result = mul_add( mul_add(coeff_a,frac,coeff_b), frac, exp_f );
	result = result * m2log2;	//multiplying by -2*log(2)
	return result;
}

//x is normalized: x=-1 is -pi radians, x=+1 is +pi radians.
//there is no checking for x outside of [-1,+1] range
//fast approximation -- accurate to 10 bits
// (might go even faster by approximating cos as 1-4*x^2)
inline void sincos_norm(const Vec8f x, Vec8f *sin, Vec8f *cos) {
	const Vec8f coeff_a =  0.99940307f ;
	const Vec8f coeff_b = -4.8911858f  ;
	const Vec8f coeff_c =  3.5838442f  ;

	const Vec8f msbbit = -0.f;

	Vec8f x_sign = x & msbbit;
	Vec8f abs_x = abs(x);
	Vec8f abs_x_minus_p5 = abs_x - 0.5f;
	Vec8i tempa = reinterpret_i(msbbit);
	Vec8i tempb = reinterpret_i(abs_x_minus_p5);
	Vec8i mag_gt_p5_i = reinterpret_i(andnot(tempa,tempb));
	Vec8f mag_gt_p5 = reinterpret_f(mag_gt_p5_i);
	//using ~ operator is very picky about arg types

	//since floating point -0.0 equals +0.0, must use int for sign-bit comparison
	//if_sub(cond,a,b) is equivalent to: result = (cond) ? a-b : a
	Vec8i msbbit_i = reinterpret_i(msbbit);
	Vec8f xcos = if_sub( mag_gt_p5_i == msbbit_i, abs_x, 1.0f );
	xcos = xcos * xcos;
	Vec8f xsin = abs_x_minus_p5 * abs_x_minus_p5;
	//lcos = coeff_a + xcos * (coeff_b + coeff_c * xcos);
	//lsin = coeff_a + xsin * (coeff_b + coeff_c * xsin);
	Vec8f lcos = mul_add( mul_add(xcos,coeff_c,coeff_b), xcos, coeff_a );
	Vec8f lsin = mul_add( mul_add(xsin,coeff_c,coeff_b), xsin, coeff_a );
	lcos = lcos | mag_gt_p5;
	*cos = lcos;
	lsin = lsin | x_sign;
	*sin = lsin;
	//NOTE: calculating sin as sqrt(1-cos^2) loses accuracy (and probably slower)
	return;
}

//Generate normal distribution of simd_count*8 number of values
//assumes seed for random number generator has been intialized
//simd_count should be even! (else will return one short)
//puts values into "result" -- caller must guarantee enough space

void normdist2( Vec8f *result, int simd_count ) {
	for (int i=0; i < simd_count/2; ++i) {
		//could replace rand_sfc64_N(0) with rand_randr() or rand_xoshiro_N(0) if desired.  Will be slower
		Vec4uq rand_val_as64 = rand_sfc64_N(0);
		Vec8i rand_val = (Vec8i)rand_val_as64;

		rand_val = rand_val & 0x7fff'ffff;	//zero the sign bit
		rand_val = rand_val ^ 0x7fff'ffff;	//invert the bits to avoid zero-land impact
		Vec8f rand_flt = to_float( rand_val );
		//instead of zeroing the sign bit, we could add 1.0f to values below zero.  Cannot just add 0.5f due to lost significance
		rand_flt = rand_flt * (0x1.0p-31f); //random value between 0 and 1
		//rand_flt is in range from 0 to 1-2^-31 in steps of 2^-31
		//if rand_val is 0000'0000 to 0000'0007 => very small number, use 64b to extend beyond ~6.6 sigma
		rand_val = rand_val & 0x7fff'fff8;

		if ( to_bits(rand_val == 0) != 0 ) {
			Vec4uq rand_val_as64 = rand_sfc64_N(0); //rand_randr() and rand_sfc64_N(0) are also options
			Vec8i rand_val = (Vec8i)rand_val_as64;
			rand_val = rand_val & 0x7fff'ffff;  //zero the sign bit
			rand_val = rand_val ^ 0x7fff'ffff;	//invert the bits to avoid zero-land impact
			Vec8f rand_flt2 = to_float( rand_val );
			rand_flt2 = mul_add(rand_flt2 , (0x1.0p-62f), (0x1.0p-63f));
			//rand_flt2 is 2^-63 to 2^-31-(2^-63) (2^-63 offset from zero)
			//Note that even the smallest value (2^-63) is not subnormal in single precision float
			rand_flt = rand_flt + rand_flt2;
		}

		Vec8f gn_term1 = sqrt( fastm2log(rand_flt) ); //fastm2log is -2.0*log(rand_flt)
		Vec4uq rand_cs_as64 = rand_sfc64_N(0); //again, rand_randr() and rand_sfc64_N(0) are options
		Vec8i rand_cs = (Vec8i)rand_cs_as64;
		Vec8f cvtm1top1 = to_float( rand_cs );
		//get random value from -1 to +1
		cvtm1top1 = cvtm1top1 * (0x1.0p-31f);
		Vec8f sinval,cosval;
		sincos_norm( cvtm1top1, &sinval, &cosval );
		*result++ = sinval * gn_term1;
		*result++ = cosval * gn_term1;
	}
}
