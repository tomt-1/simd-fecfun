function [total_f2b_b2f] = gen_decode( Hc, elem_per_simd, z_value )

Hlog = 1.0*(Hc >= 0);	%H logical.  1=> submatrix is not ZxZ zero array
[mrow mcol] = size(Hc);

z_remainder = mod(z_value,elem_per_simd);
z_perfect_flag = (mod(z_value,elem_per_simd) == 0);

row_active_cols = [];
row_offset = [];
for r=1:mrow
	cols_this_row = find(Hlog(r,:));
	row_active_cols = [row_active_cols cols_this_row];
	row_weights(r) = length(cols_this_row);
	row_offset = [row_offset Hc(r,cols_this_row)];
	%subtract 1 for 0-based indexing on row_active_cols
end
maxrow = max(row_weights);  %maximum row weight

total_Hc_ones = sum(sum(Hlog));

f2b_copy_cnt = [];
b2f_copy_cnt = [];
f2b_offsets = [];
b2f_offsets = [];
f2b_cols = [];
b2f_cols = [];
all_prev_offsets = zeros(1,mcol);
for r=1:mrow
	this_row = find(Hlog(r,:));
	offsets = Hc(r,this_row);
	prev_offsets = all_prev_offsets(this_row);
	if ( z_perfect_flag ) %same calc as below, since z is perfect.  consider simplifying code
		offsets = mod(offsets,elem_per_simd);
	else
		offsets(offsets >= z_remainder) = ...
			mod(offsets(offsets >= z_remainder)-z_value,elem_per_simd);
	end
	deltas = offsets - prev_offsets;
	%+ve => copy front to back, -ve => copy back to front, 0 => no copy
	%copy based on prev_offsets
	pos_delta_mask = (deltas > 0);
	neg_delta_mask = (deltas < 0);

	f2b_copy_cnt = [f2b_copy_cnt sum(pos_delta_mask)];
	b2f_copy_cnt = [b2f_copy_cnt sum(neg_delta_mask)];
	f2b_offsets = [f2b_offsets prev_offsets(pos_delta_mask)];
	b2f_offsets = [b2f_offsets prev_offsets(neg_delta_mask)];
	f2b_cols = [f2b_cols this_row(pos_delta_mask)];
	b2f_cols = [b2f_cols this_row(neg_delta_mask)];
	all_prev_offsets(this_row) = offsets;
end
%Final b2f copies for syndrome check or start of next iteration
deltas = -all_prev_offsets;
neg_delta_mask = (deltas < 0);
b2f_copy_cnt = [b2f_copy_cnt sum(neg_delta_mask)];
b2f_offsets = [b2f_offsets all_prev_offsets(neg_delta_mask)];
all_rows = (1:mcol);
b2f_cols = [b2f_cols all_rows(neg_delta_mask)];

fs = fopen('LDPCdecode.raw.dat','w');
fwrite(fs,[mrow mcol z_value total_Hc_ones maxrow length(f2b_cols) length(b2f_cols)],'int32');
fwrite(fs,[row_weights row_active_cols-1 row_offset f2b_copy_cnt f2b_cols-1 f2b_offsets ...
		   b2f_copy_cnt b2f_cols-1 b2f_offsets],'int32');
fclose(fs);
total_f2b_b2f = sum(f2b_copy_cnt) + sum(b2f_copy_cnt);
fprintf('total f2b+b2f elements: %d\n',total_f2b_b2f);

