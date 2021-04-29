# simd-fecfun: Using SIMD to accelerate FEC-related functions

This is a set of functions for forward error correction (FEC).  The code uses Agner Fog's vector class library (https://github.com/vectorclass - TBD: add as git submodule)

This is definitely a work-in-progess, for a more polished FEC project, one might want to look at aff3ct (https://github.com/aff3ct/aff3ct).  The key difference here is to push the simulation and decoding rate higher.

Key Components:
1. [WGN generation](https://tomt-1.github.io/simd-fecfun/WGN/norm.dist.methods.html)
2. LDPC Encoding (TBD)
3. LDPC Decoding (TBD) 
