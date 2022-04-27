function [direct_cols subst_row subst_col inv_method inv_filename z_value] = encode_method(MatrixSet);

inv_file_exists = 0;
inv_filename = '';
inv_method = 0; %2 is gf-based inversion

switch MatrixSet
	case '11ad 1/2'
		direct_cols = [];
		subst_row = [1 2 3 4 5 6 7 8];
		subst_col = [1 2 3 4 5 6 7 8];
		z_value = 42;

	case '11ad 5/8'
		direct_cols = [];
		subst_row = [1 2 3 4 5 6];
		subst_col = [1 2 3 4 5 6];
		z_value = 42;

	case '11ad 3/4'
		direct_cols = [];
		subst_row = [1 2 3 4];
		subst_col = [1 2 3 4];
		z_value = 42;

	case '11ad 13/16'
		direct_cols = [];
		subst_row = [1 2 3];
		subst_col = [1 2 3];
		z_value = 42;

	case '11ay 1/2'
		direct_cols = [];
		subst_row = (1:16);
		subst_col = (1:16);
		z_value = 42;

	case '11ay 5/8'
		direct_cols = [];
		subst_row = (1:12);
		subst_col = (1:12);
		z_value = 42;

	case '11ay 3/4'
		direct_cols = [];
		subst_row = (1:8);
		subst_col = (1:8);
		z_value = 42;

	case '11ay 13/16'
		direct_cols = [];
		subst_row = (1:6);
		subst_col = (1:6);
		z_value = 42;

	case '11ay 7/8'
		direct_cols = [];
		subst_row = (1:4)
		subst_col = [2 1 3 4];
		z_value = 42;

	case '16e 1/2'
		direct_cols = [1];
		subst_row = [1 2 3 4 5 6 7 8  9 10 11];
		subst_col = [2 3 4 5 6 7 8 9 10 11 12];
		z_value = 24;  %this can be changed

	case 'ITU G.hn 1/2 z14'
		direct_cols = [1];
		subst_row = [1 2 3 4 5 6 7 8 9 10 11];
		subst_col = [2 3 4 5 6 7 8 9 10 11 12];
		z_value = 14;

	case '802.3ca'
		direct_cols = [8 11 12];
		subst_row = [2 6 12  1 4 5 7 11 8];
		subst_col = [4 2  3 10 7 1 9  6 5];
		inv_method = 1;  %1 => 'invM2' mat file exists
		inv_filename = 'invM2_size256';
		z_value = 256;

	case '16e 2/3A'
		direct_cols = [1];
		subst_row = [1 2 3 4 5 6 7];
		subst_col = [2 3 4 5 6 7 8];
		z_value = 24;

	case '16e 2/3B'
		direct_cols = [1];
		subst_row = [1 2 3 4 5 6 7];
		subst_col = [2 3 4 5 6 7 8];
		z_value = 24;

	case '16e 3/4A'
		direct_cols = [1];
		subst_row = [1 2 3 4 5];
		subst_col = [2 3 4 5 6];
		z_value = 24;

	case '16e 3/4B'
		direct_cols = [1];
		subst_row = [1 2 3 4 5];
		subst_col = [2 3 4 5 6];
		z_value = 24;

	case '16e 5/6'
		direct_cols = [1];
		subst_row = [1 2 3];
		subst_col = [2 3 4];
		z_value = 24;

	case '11n 1/2 z27'
		direct_cols = [1];
		subst_row = [1 2 3 4 5 6 7 8  9 10 11];
		subst_col = [2 3 4 5 6 7 8 9 10 11 12];
		z_value = 27;

	case '11n 2/3 z27'
		direct_cols = [1];
		subst_row = [1 2 3 4 5 6 7];
		subst_col = [2 3 4 5 6 7 8];
		z_value = 27;

	case '11n 3/4 z27'
		direct_cols = [1];
		subst_row = [1 2 3 4 5];
		subst_col = [2 3 4 5 6];
		z_value = 27;

	case '11n 5/6 z27'
		direct_cols = [1];
		subst_row = [1 2 3];
		subst_col = [2 3 4];
		z_value = 27;

    case '11n 1/2 z54'
        direct_cols = [1];
        subst_row = [1 2 3 4 5 6 7 8  9 10 11];
        subst_col = [2 3 4 5 6 7 8 9 10 11 12];
        z_value = 54;

    case '11n 2/3 z54'
        direct_cols = [1];
        subst_row = [1 2 3 4 5 6 7];
        subst_col = [2 3 4 5 6 7 8];
        z_value = 54;

    case '11n 3/4 z54'
        direct_cols = [1];
        subst_row = [1 2 3 4 5];
        subst_col = [2 3 4 5 6];
        z_value = 54;

    case '11n 5/6 z54'
        direct_cols = [1];
        subst_row = [1 2 3];
        subst_col = [2 3 4];
        z_value = 54;

    case '11n 1/2 z81'
        direct_cols = [1];
        subst_row = [1 2 3 4 5 6 7 8  9 10 11];
        subst_col = [2 3 4 5 6 7 8 9 10 11 12];
        z_value = 81;

    case '11n 2/3 z81'
        direct_cols = [1];
        subst_row = [1 2 3 4 5 6 7];
        subst_col = [2 3 4 5 6 7 8];
        z_value = 81;

    case '11n 3/4 z81'
        direct_cols = [1];
        subst_row = [1 2 3 4 5];
        subst_col = [2 3 4 5 6];
        z_value = 81;

    case '11n 5/6 z81'
        direct_cols = [1];
        subst_row = [1 2 3];
        subst_col = [2 3 4];
        z_value = 81;

	otherwise
		if ( MatrixSet(1:5) == "5GNR_" )
			direct_cols = [1 3];
			subst_row = [2 4 (5:99)]; %these will be truncated outside the function to correct size
			subst_col = [2 4 (5:99)]; 
			z_value = 0; %Z is set outside the function
		else
			display 'invalid MatrixSet argument - Error';
		end

end
