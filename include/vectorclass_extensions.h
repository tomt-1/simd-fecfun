
//included file for additional inline functions

//rotate left (from LSB to MSB), treating 256b vector as LSB at lowest address
//aka little-endian, just like other IA items
inline Vec4uq rot256_left1(const Vec4uq a) {
	Vec4uq result,tmp1,tmp2;
	tmp1 = (a << 1);
	tmp2 = (a >> 63);
	tmp2 = permute4<3,0,1,2>(tmp2);
	result = tmp1 | tmp2;
	return result;
}

//rotate right (from MSB to LSB), treating 256b vector as little-endian
inline Vec4uq rot256_right1(const Vec4uq a) {
	Vec4uq result,tmp1,tmp2;
	tmp1 = (a >> 1);
	tmp2 = (a << 63);
	tmp2 = permute4<1,2,3,0>(tmp2);
	result = tmp1 | tmp2;
	return result;
}

//shift left 1 (same as rot left, but lsb is zero)
inline Vec4uq shift256_left1(const Vec4uq a) {
	Vec4uq result,tmp1,tmp2;
	tmp1 = (a << 1);
	tmp2 = (a >> 63);
	tmp2 = permute4<-1,0,1,2>(tmp2);
	result = tmp1 | tmp2;
	return result;
}

//shift right 1 (same as rotate right, but msb is zero)
inline Vec4uq shift256_right1(const Vec4uq a) {
	Vec4uq result,tmp1,tmp2;
	tmp1 = (a >> 1);
	tmp2 = (a << 63);
	tmp2 = permute4<1,2,3,-1>(tmp2);
	result = tmp1 | tmp2;
	return result;
}

//extending sign_combine for 8/16/32b integer cases
//NOTE: no element of vector "a" should hold a "most negative value"
//      e.g.: if using 8-bit elements, a value of -128 should be avoided
//use 256b _sign_intrinsic where applicable, else "select( (b<0),-a,a )"
inline Vec32c sign_combine(const Vec32c a, const Vec32c b) {
	#if INSTRSET > 7 // AVX2 and later
	Vec32c b_valid = b | 1; //make sure b is not zero
	Vec32c result = _mm256_sign_epi8( a, b_valid );
	#else
	Vec32c result = select( (b < 0),-a,a );
	#endif
	return result;
}

inline Vec8i sign_combine(const Vec8i a, const Vec8i b) {
	#if INSTRSET > 7 // AVX2 and later
	Vec8i b_valid = b | 1;	//make sure b is not zero
	Vec8i result = _mm256_sign_epi32( a, b_valid );
	#else
	Vec8i result = select( (b < 0),-a,a );
	#endif
	return result;
}

inline Vec16s sign_combine(const Vec16s a, const Vec16s b) {
	#if INSTRSET > 7 // AVX2 and later
	Vec16s b_valid = b | 1;	//make sure b is not zero
	Vec16s result = _mm256_sign_epi16( a, b_valid );
	#else
	Vec16s result = select( (b < 0),-a,a );
	#endif
	return result;
}

inline Vec64c sign_combine(const Vec64c a, const Vec64c b) {
	Vec64c result = select( (b < 0),-a,a );
	return result;
}

inline Vec32s sign_combine(const Vec32s a, const Vec32s b) {
	Vec32s result = select( (b < 0),-a,a );
	return result;
}

inline Vec16i sign_combine(const Vec16i a, const Vec16i b) {
	Vec16i result = select( (b < 0),-a,a );
	return result;
}
