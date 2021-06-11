#!/usr/bin/perl -w

#expected to simply pipe stdout to all_min_functions.cpp when running this script
#always have at least two items in min_list

@min_list = (2 .. 53); #this is all supported min functions
#@min_list = (22 .. 23); #this is only 802.3ba matrix

$start_idx = 0;

print "#include \"vectorclass.h\"\n";
if ($ARGV[0] eq 'prod') {
	#Get SVML-based header if doing atanh/tanh products
	#print "#include \"vectormath_lib.h\"\n";
	print "#include \"vectormath_hyp.h\"\n";
}
print "#include \"vectorclass_extensions.h\"\n";
print "#include \"simd_defs.h\"\n\n";

if ($ARGV[0] eq 'prod') {
	print "//Original code defined MAXVALs for tanh argument to prevent infinite atanh.\n";
	print "//this is now handled by requiring the max_metric to be set appropriately\n";
	print "//(max metric <=18.02 for floats, <=38.12 for doubles)\n\n";
	#print "#if SIMD_ELEM_BW == -32\n";
	#print "\tconst SIMD_CLASS MAXVAL = 9.0f;    //prevent atanh from going +INF\n";
	#print "#elif SIMD_ELEM_BW == -64\n";
	#print "\tconst SIMD_CLASS MAXVAL = 19.0l;   //prevent atanh from going +INF\n";
	#print "#endif\n\n";
}

for $idx (@min_list[$start_idx .. $#min_list]) {
	if ( $ARGV[0] eq 'prod' ) {
		system("./min_find.pl $idx prod");
		$type='prod';
	} else {
		system("./min_find.pl $idx min");
		$type='min';
	}
}

$last_idx = $min_list[$#min_list];
print "typedef void (*min_fn_ptr_type)(SIMD_CLASS *,SIMD_CLASS *,SIMD_CLASS,int);\n\n";
print "#define FUNC_ARY THREE_CAT($ARGV[0],_fnary_,SIMD_CLASS)\n";
print "min_fn_ptr_type FUNC_ARY[$last_idx+1] = {NULL,NULL, ";

$min_idx = 0;
for $idx (2 .. $last_idx-1) {
	if ( $min_list[$min_idx] == $idx ) {
		print "&${type}_find_$idx, ";
		++$min_idx;
	} else {
		print "NULL, ";
	}
}
print "&${type}_find_$last_idx };\n";
#allow for us to change last_idx in only this script
print "#define MIN_FUNCS THREE_CAT($ARGV[0],_funcs_,SIMD_CLASS)\n";
#line below is can be used if building a library of all SIMD_CLASS options
#print "min_fn_ptr_type *MIN_FUNCS = FUNC_ARY;\n";
print "min_fn_ptr_type *all_minprod_func_ary = FUNC_ARY;\n";

