%start is similar to encode script in build
addpath('../matlab/LDPC');
MatrixSet_cells = {  ...
    '11ad 1/2','11ad 5/8','11ad 3/4', '11ad 13/16', ...
    '16e 1/2','16e 2/3A', '16e 2/3B', '16e 3/4A', '16e 3/4B', '16e 5/6', ...
    '11n 1/2 z27', '11n 2/3 z27', '11n 3/4 z27', '11n 5/6 z27', ...
    '11n 1/2 z54', '11n 2/3 z54', '11n 3/4 z54', '11n 5/6 z54', ...
    '11n 1/2 z81', '11n 2/3 z81', '11n 3/4 z81', '11n 5/6 z81', ...
    'ITU G.hn 1/2 z14', '802.3ca' };

if (~exist('MatrixSet','var'))
    MatrixSet = '11ad 1/2';
end
if (~exist('num_par','var'))
    num_par = 1;
end
if ( ~exist('simd_size','var') )
    simd_size = 256;    %SIMD word size in bits.  256 or 512
end
if ( ~exist('elem_size','var') )
    elem_size = 32;     %SIMD element (or component) size in bits
end

if ( exist('OCTAVE_VERSION', 'builtin') )
	warning('off','Octave:data-file-in-path')
end

[Hc,Z_GetH,extra_5g] = GetHMatrix(MatrixSet);
[direct_cols subst_row subst_col inv_method inv_filename z_value] = encode_method(MatrixSet);
if ( MatrixSet(1:5) == "5GNR_" )
	par_mcol = extra_5g(1);
	punc_fill_lists = extra_5g(2:end);
	%adjust for num_par >1
	punc_fill_lists(3:2:end) = (punc_fill_lists(3:2:end)-1)*num_par+1;
	punc_fill_lists(4:2:end) = punc_fill_lists(4:2:end)*num_par;
	%and make indexing 0-based instead of 1-based
	punc_fill_lists(3:end) = punc_fill_lists(3:end)-1;
	subst_row = subst_row(1:par_mcol-2); %total parity mcol, minus the 2 in direct_col
	subst_col = subst_col(1:par_mcol-2);
	z_value = Z_GetH;  %Z is determined by GetHMatrix if 5GNR
else
	punc_fill_lists = [0 0];
end

z_value = z_value * num_par;
Hc = Hc*num_par;
Hc(Hc<0) = -1; %restore -1 values to -1

if (exist('z_over','var'))
    z_value = z_over; %override of actual Z
end

%this will create the encoder structure file LDPCencode.raw.dat
[H NumInfoBits invM2tM1] = gen_encode(Hc,z_value,num_par,inv_method,inv_filename,direct_cols,subst_row,subst_col,punc_fill_lists);

bypass = 0; %1=> bypass re-ordering within major rows
%derived params
elem_per_simd   = (simd_size / elem_size);  %elements per SIMD word
                                            %aka SIMD parallelism
simd_per_z      = ceil( z_value / elem_per_simd );  %SIMD words in Z

Hc = decode_method(MatrixSet,Hc,z_value,elem_per_simd,bypass);
%this will create the decoder structure file LDPCdecode.raw.dat
total_f2b_b2f = gen_decode(Hc,elem_per_simd, z_value );
