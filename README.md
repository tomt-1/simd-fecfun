# simd-fecfun: Using SIMD to accelerate FEC-related functions

This is a set of functions for forward error correction (FEC).  The code uses Agner Fog's vector class library (https://github.com/vectorclass - at a minimum get the "version2-master" files)

This is a work-in-progess, for a more polished FEC project, one might want to look at aff3ct (https://github.com/aff3ct/aff3ct).  The key difference here is to push the simulation and decoding rate higher.


## Key Components

1. [WGN generation](https://tomt-1.github.io/simd-fecfun/WGN/norm.dist.methods.html)
2. LDPC Encoding (TBD)
3. LDPC Decoding (TBD) 


## Prerequisites

Only tested on Ubuntu, but other Linux distributions should work.  Requires Intel or AMD processor with at
least SSE support, although AVX or AVX-512 is recommended.

Above-mentioned vector class library and environment variable pointing to it:
`VECINC=<path to vector class lib>/version2-master`

Either octave or Matlab on the path.  If using octave, set environment variable:
`USEOCT=1`

common programming tools:
g++/make/perl

## Build and Run

```bash
git clone https://github.com/tomt-1/simd-fecfun
mkdir simd-fecfun/build/tmp
cd simd-fecfun/build
source commands.build-all.src
./run_LDPCsim.pl num_proc=1
```
