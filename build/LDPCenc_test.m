%run from 'build' directory.  But, script will cd down to 'tmp'
cd('tmp');
addpath('../../matlab/LDPC');
MatrixSet_cells = {  ...
	'11ad 1/2','11ad 5/8','11ad 3/4', '11ad 13/16', ...
	'16e 1/2','16e 2/3A', '16e 2/3B', '16e 3/4A', '16e 3/4B', '16e 5/6', ...
	'11n 1/2 z27', '11n 2/3 z27', '11n 3/4 z27', '11n 5/6 z27', ...
	'11n 1/2 z54', '11n 2/3 z54', '11n 3/4 z54', '11n 5/6 z54', ...
	'11n 1/2 z81', '11n 2/3 z81', '11n 3/4 z81', '11n 5/6 z81', ...
	'ITU G.hn 1/2 z14', '802.3ca' };
fid_datfile = 'enc_dat.txt';
fid_parfile = 'enc_par.txt';
fid_newfile = 'new_par.txt'; %file created by compiled encoder
rand('seed',sum(100*clock));
%for octave, one may want to suppress warning on where .mat file is found:
%warning('off','Octave:data-file-in-path');

%num_cw_per must match compiled code (line 17 of LPDCencode_example.cpp).  num_idx can be changed as desired
num_idx = 10;
num_cw_per = 10;

for i=1:num_idx
	MatrixIdx = randi(length(MatrixSet_cells), 1, 1);
	MatrixSet = MatrixSet_cells{MatrixIdx};
	Hc = GetHMatrix(MatrixSet);
	[direct_cols subst_row subst_col inv_method inv_filename z_value] = encode_method(MatrixSet);
	max_par = floor(256 / z_value);
	num_par = randi(max_par, 1, 1);
	z_value = z_value * num_par;
	Hc = Hc*num_par;
	Hc(Hc<0) = -1; %restore -1 values to -1, where -1 implies ZxZ zero submatrix

	[H NumInfoBits invM2tM1] = gen_encode(Hc,z_value,num_par,inv_method,inv_filename,direct_cols,subst_row,subst_col);

	fid_dat = fopen(fid_datfile,'w+');
	fid_par = fopen(fid_parfile,'w+');
	for cw=1:num_cw_per
		info = randi([0 1],NumInfoBits,1);
		p = mod(invM2tM1 * info,2);
		HxCW = mod( H * [info;p], 2 );
		if ( sum(abs(HxCW)) ~= 0 )
			display 'Error - computed parity does not give H times Codeword = zero-vector';
		end
		[mrow mcol] = size(Hc);
		for col=1:mcol-mrow
			col_data = info((col-1)*z_value+1:col*z_value);
			fprintf(fid_dat,'%d',col_data);
			fprintf(fid_dat,'\n');
		end
		for col=1:mrow
			col_data = p((col-1)*z_value+1:col*z_value);
			fprintf(fid_par,'%d',col_data);
			fprintf(fid_par,'\n');
		end
	end
	fclose(fid_dat);
	fclose(fid_par);

	%run in 'tmp', so executable is one up
	system(strcat(['../LDPCencode_example < ',fid_datfile,'| grep -v Matrix > ',fid_newfile]));
	fprintf('Testing MatrixSet %s, num_par=%d, Z=%d\n',MatrixSet,num_par,z_value);
	system(strcat(['diff -s --brief ',fid_newfile,' ',fid_parfile]));
end
