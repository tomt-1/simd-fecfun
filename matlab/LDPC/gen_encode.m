function [H NumInfoBits invM2tM1] = gen_encode(Hc,z_value,num_par,inv_method,inv_filename,direct_cols,subst_row,subst_col,punc_fill_list)
%create encode data for compiled function.  Return items so that
%example test vectors can be created and checked

%TODO: move simd_bw to an input of the function
simd_bw = 512;
byte_per_simd = simd_bw / 8;

Hlog = 1.0*(Hc >= 0);   %H logical.  1=> submatrix is not ZxZ zero array
[mrow mcol] = size(Hc);

[H MaxInRow MaxInCol] = ExpandH(Hc,z_value);

[tmprow tmpcol] = size(H);
NumInfoBits = tmpcol - tmprow;
M1 = H(:,1:NumInfoBits);
M2 = H(:,NumInfoBits+1:end);

if (inv_method == 1)
	load(inv_filename);  %name of matrix loaded must be 'invM2'
elseif (inv_method == 2)
    M2gf = gf(M2);
    invM2gf = inv(M2gf);  	%this can take a long time for large arrays
    invM2 = (invM2gf == 1); %back to real numbers
    save('recent_invM2','invM2');
else
	invM2 = mod(floor(inv(M2)+1e-12),2); %only works for simple M2
end
invM2tM1 = mod(invM2*M1,2);

%get the directly calculated columns of parity bits
row1_direct = invM2tM1(z_value*(direct_cols-1)+1,:); %assuming direct_cols is not 0-based
par_offset = [];
col_cnts = zeros(1,length(direct_cols));

for ii = 1:length(direct_cols)
	shift_idx = find(row1_direct(ii,:))-1;  %zero-based indexes of "1" values
	col_nums = floor( shift_idx / z_value );
	shift_idx = shift_idx - col_nums * z_value;

	offset = floor(shift_idx / 8);
	mod8 = bitand(shift_idx,7);

	wa_addr = col_nums*16*byte_per_simd + mod8*2*byte_per_simd + offset;
	%might want to sort the addrs to maximize cacheing
	par_offset = [par_offset wa_addr];
	col_cnts(ii) = length(wa_addr);
end

%substitution parity calcs

for r_idx=1:length(subst_row)
	tmp = find(Hlog(subst_row(r_idx),:));
	tmp = tmp(tmp ~= subst_col(r_idx)+mcol-mrow);
	col_cnts = [col_cnts length(tmp)];
	rotate_target = Hc(subst_row(r_idx),subst_col(r_idx)+mcol-mrow);
	row_offset = Hc(subst_row(r_idx),tmp);
	row_offset = mod(z_value-rotate_target+row_offset,z_value); %put substitution col to zero

	offset = floor(row_offset / 8);
	mod8 = bitand(row_offset,7);

	wa_addr = (tmp-1)*16*byte_per_simd + mod8*2*byte_per_simd + offset;
	par_offset = [par_offset wa_addr];
end

par_colnum = [direct_cols subst_col]-1; %make both 0-based for data file

%deal with puncture and fill locations
if ( length(punc_fill_list) > 0 )
	punc_cnt = punc_fill_list(1);
	fill_cnt = punc_fill_list(2);
	punc_list = punc_fill_list(3:2+punc_cnt);
	fill_list = punc_fill_list(3+punc_cnt:end);
else
	punc_cnt = 0;
	fill_cnt = 0;
	punc_list = [];
	fill_list = [];
end

%write out all information
fs = fopen('LDPCencode.raw.dat','w');
fwrite(fs,[z_value punc_cnt fill_cnt num_par mrow mcol col_cnts par_colnum par_offset punc_list fill_list],'int32');
fclose(fs);
