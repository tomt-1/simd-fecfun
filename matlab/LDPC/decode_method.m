function Hc_reordered = decode_method(MatrixSet,Hc,z_value,elem_per_simd,bypass);

z_remainder = mod(z_value,elem_per_simd);
z_perfect_flag = (mod(z_value,elem_per_simd) == 0);

%Can add here for special row re-ordering of various matrix sets
if ( isequal(MatrixSet,'802.3ca') )
	%put row 2 first so that first mrow eval can correct about half
	%of the punctured bits
	tmp = Hc(1,:);
	Hc(1,:) = Hc(2,:);
	Hc(2,:) = tmp;
end

Hlog = (Hc >= 0);	%H logical.  1=> submatrix is not ZxZ zero array

[mrow mcol] = size(Hc);

%re-sort rows to minimize need for f2b and b2f copying
Hc_save = Hc; %temporary bypass
all_prev_offsets = zeros(1,mcol);
for row=1:mrow
	this_row = find(Hlog(row,:));
	prev_offsets = all_prev_offsets(this_row);
	max_no_copy = -1;
	max_no_copy_offsets = [];
	max_no_copy_Hrow = [];
	for i=0:z_value-1 %if z is perfect, the high limit could be changed to elem_per_simd-1
		new_Hrow = Hc(row,this_row) + i;
		new_Hrow = mod(new_Hrow,z_value);
		offsets = new_Hrow;
		offsets(offsets >= z_remainder) = mod(offsets(offsets >= z_remainder)-z_value,elem_per_simd);
		no_copy_count = sum( offsets == prev_offsets );
		if ( no_copy_count > max_no_copy )
			max_no_copy = no_copy_count;
			max_no_copy_offsets = offsets;
			max_no_copy_Hrow = new_Hrow;
		end
	end
	Hc(row,this_row) = max_no_copy_Hrow;
	all_prev_offsets(this_row) = max_no_copy_offsets;
end
Hc_reordered = Hc;
if (bypass) %bypass the re-sorting
	Hc_reordered = Hc_save;
end
