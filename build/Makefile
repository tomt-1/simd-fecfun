#ARCHFLAG = -mavx512f -mavx512bw -mfma
#ARCHFLAG = -mavx -mavx2 -mfma
ARCHFLAG = -march=native

ifndef VECINC
$(error VECINC is not set - should point to version2-master directory of github.com/vectorclass clone)
endif

#override these defaults for different options
#for SIMD_ELEM: -32 is float, -64 is double.  postive values are integers
SIMD_PAR = 64
SIMD_ELEM = 8
MINSUM = 1
PROG_NAME = LDPCtop_Vec64c_min

ifeq ($(DEBUG),1)
	CXXFLAGS = -g -O0 -Wall -std=c++17 $(ARCHFLAG)
	CXXFLAGD = -g -O0 -Wall -std=c++17 $(ARCHFLAG)
else
	#-flto, -ffast-math: no benefit.  -funroll-loops on Decoder gives 3.6% speed improvement
	#CXXFLAGD is for the decoder 
	CXXFLAGS = -O3 -Wall -std=c++17 $(ARCHFLAG)
	CXXFLAGD = -O3 -funroll-loops -Wall -std=c++17 $(ARCHFLAG)
endif

ifeq ($(MINSUM),1)
	PRODMIN_OBJ = all_min_generated.o
	SVML_OPT = 
else
	PRODMIN_OBJ = all_prod_generated.o
	#if using SVL hacked version, uncomment below line and also need LD_LIBRARY_PATH=$(MKLROOT)
	#SVML_OPT = -L$(MKLROOT) -lsvml -lirc
	SVML_OPT =
endif

#Decode, Encode, and Top objects
DOBJ = decode_matrix_allocate.o $(PRODMIN_OBJ) LDPCdecode.o
EOBJ = encode_matrix_allocate.o LDPCencode.o
TOBJ = LDPC_count_errors.o LDPCtop.o

#.PHONY: all clean $(PROG_NAME) compile
.PHONY: all
all: $(PROG_NAME)

../../src/LDPCDec/all_min_generated.cpp  : ../../src/LDPCDec/min_find_wrap.pl ../../src/LDPCDec/min_find.pl
	cd ../../src/LDPCDec && ./min_find_wrap.pl min  > all_min_generated.cpp  && cd ../../build/tmp

../../src/LDPCDec/all_prod_generated.cpp : ../../src/LDPCDec/min_find_wrap.pl ../../src/LDPCDec/min_find.pl
	cd ../../src/LDPCDec && ./min_find_wrap.pl prod > all_prod_generated.cpp && cd ../../build/tmp

$(DOBJ) : %.o : ../../src/LDPCDec/%.cpp
	$(CXX) $(CXXFLAGD) -I$(VECINC) -I../../include -DSIMD_PAR=$(SIMD_PAR) -DSIMD_ELEM_BW=$(SIMD_ELEM) -c $< -o $@

prng_functions.o : ../../src/RNGfunc/prng_functions.cpp
	$(CXX) $(CXXFLAGS) -mrdrnd -I$(VECINC) -DSIMD_F=Vec8f -DSIMD_I=Vec8i -DSIMD_RAND=Vec4uq -c $< -o $@

$(EOBJ) : %.o : ../../src/LDPCEnc/%.cpp
	$(CXX) $(CXXFLAGS) -I$(VECINC) -I../../include -c $< -o $@

#always uses Vec4uq
gen_rand_data_bits.o : ../../src/Top/gen_rand_databits.cpp
	$(CXX) $(CXXFLAGS) -I$(VECINC) -I../../src/LDPCEnc -I../../include -DSIMD_PAR=4 -DSIMD_ELEM_BW=64 -c $< -o $@

read_sim_params.o : ../../src/Top/read_sim_params.cpp
	$(CXX) $(CXXFLAGS) -mrdrnd -c $< -o $@

gen_metrics.o : ../../src/Top/gen_metrics.cpp
	$(CXX) $(CXXFLAGS) -I$(VECINC) -I../../include -DSIMD_PAR=$(SIMD_PAR) -DSIMD_ELEM_BW=$(SIMD_ELEM) -c $< -o $@

$(TOBJ) : %.o : ../../src/Top/%.cpp
	$(CXX) $(CXXFLAGS) -mrdrnd -I$(VECINC) -I../../include -DSIMD_PAR=$(SIMD_PAR) -DSIMD_ELEM_BW=$(SIMD_ELEM) -c $< -o $@


$(PROG_NAME): $(DOBJ) prng_functions.o $(EOBJ) gen_rand_data_bits.o read_sim_params.o gen_metrics.o $(TOBJ)
	$(CXX) $^ $(SVML_OPT) -o ../$@
