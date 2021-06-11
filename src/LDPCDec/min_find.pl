#!/usr/bin/perl -w

#sequence of min(a,b) functions to get the minimum for
#every element in a list excluding the "pivot" item
#as non-branching instructions
#  (so, if num_elem is 4 and these were labelled a,b,c,d
#   then result is mab = min(a,b), mcd=min(c,d)
#   mabc = min(mab,c) mabd=min(mab,d) macd=min(mcd,a)
#   and mbcd=min(mcd,b).  Requiring 6 min(a,b) evaluations)

#no guarantee this is optimal.  Feel free to replace
#step 1: progress down tree until only 2 elements
#        try to balance the number of items in args
#        (using sort to ensure last odd item isn't alone)
#
#step 2: combine items back up tree until all mins computed
#        if any "curr" row item unused, it gets copied to next
#        when "combine" row is at top, only keep finished items

#quick-and-dirty code to compute this.  Likely could be improved
#limit of 53 elements due to mantissa of double-precision floating pt.

#arguments are $num_elem and either 'min' or 'prod'
$num_elem=$ARGV[0];
$fn_type=$ARGV[1]; #either 'min' or 'prod' (min-sum or sum-product)

$hit = int(($num_elem-0.1)/4)+1; #hex digits in $num_elem bits

print "void ${fn_type}_find_$num_elem ( SIMD_CLASS *grp_min, SIMD_CLASS *fin_min, SIMD_CLASS Max_val, int SIMD_PER_Z ) {\n";

print "  for (int i=0; i < SIMD_PER_Z; ++i) {\n";
$res = '';
for $i (0 .. $num_elem-1) {
	push(@{$cnt1s[0]},1);
	push(@{$mask[0]},2**$i);
	$tmp = 2**$i;
	if ( $fn_type eq 'prod' ) {
		#printf "\tSIMD_CLASS m%0${hit}x = tanh( min(MAXVAL, grp_min[%d*SIMD_PER_Z + i] * 0.5) );\n",$tmp,$i
		printf "\tSIMD_CLASS m%0${hit}x = tanh( grp_min[%d*SIMD_PER_Z + i] * 0.5 );\n",$tmp,$i
	} else {
		printf "\tSIMD_CLASS m%0${hit}x = grp_min[%d*SIMD_PER_Z + i];\n",$tmp,$i;
	}
	$all1s = (2 ** $num_elem)-1;
	if ( $fn_type eq 'prod' ) {
		$res .= sprintf "\tfin_min[%d*SIMD_PER_Z + i] = 2.0*atanh( m%0${hit}x );\n",$i,$all1s ^ $tmp;
	} else {
		$res .= sprintf "\tfin_min[%d*SIMD_PER_Z + i] = m%0${hit}x;\n",$i,$all1s ^ $tmp;
	}
};
$res .= "  }\n}\n";

$lvl=1;
do {
	$num_elem_prev = scalar(@{$mask[$lvl-1]});
	$num_elem_d2 = ($num_elem_prev) >> 1;
	for $i (0 .. $num_elem_d2-1) {
		$mask = $mask[$lvl-1][2*$i] | $mask[$lvl-1][2*$i+1];
		$cnt = $cnt1s[$lvl-1][2*$i] + $cnt1s[$lvl-1][2*$i+1];
		push(@{$mask[$lvl]},$mask);
		push(@{$cnt1s[$lvl]},$cnt);
		$m1=$mask[$lvl-1][2*$i];
		$m2=$mask[$lvl-1][2*$i+1];
		if ( $num_elem != 2 ) {
			if ( $fn_type eq 'prod' ) {
				printf "\tSIMD_CLASS m%0${hit}x = m%0${hit}x * m%0${hit}x;\n",$mask,$m1,$m2;
			} else {
				printf "\tSIMD_CLASS m%0${hit}x = min( m%0${hit}x, m%0${hit}x );\n",$mask,$m1,$m2;
			}
		}
	}
	if ( $num_elem_prev & 1 ) {
		++$num_elem_d2;
		push(@{$mask[$lvl]},$mask[$lvl-1][$num_elem_prev-1]);
		push(@{$cnt1s[$lvl]},$cnt1s[$lvl-1][$num_elem_prev-1]);
	}
	#smooth the combined count-of-1s as much as possible
	#there might be a better way than sorting
	@idx = (0 .. $#{$mask[$lvl]});
	@sortidx = sort { $cnt1s[$lvl][$a] <=> $cnt1s[$lvl][$b] } @idx;
	@tmpm=();
	@tmpc=();
	for (0 .. $#sortidx) {
		$tmpm[$_] = $mask[$lvl][$sortidx[$_]];
		$tmpc[$_] = $cnt1s[$lvl][$sortidx[$_]];
	}
	@{$mask[$lvl]} = @tmpm;
	@{$cnt1s[$lvl]} = @tmpc;
	++$lvl;
} while ($num_elem_d2 > 2);

#now, comb last row with all elements of row 2 and higher
--$lvl;
$not_done = 1;
$comb_lvl = $lvl-1;

@currm = @{$mask[$lvl]};
@currc = @{$cnt1s[$lvl]};

#apply max allowed values here - the minimum number of arguments
for $i (0 .. 1) { #should *always* be 2 items (except for $num_elem == 2)
	if ( $num_elem == 2 ) {
		$tmp = $i+1;
	} else {
		$tmp = $currm[$i];
	}
	printf "\tm%0${hit}x = min( Max_val, m%0${hit}x );\n",$tmp,$tmp;
}

$end_cnt=0;
%allkey = ();
do {
	@nextm=();
	@nextc=();
	%key_thispass = ();
	$num_elem_prev = scalar(@{$mask[$comb_lvl]});
	$num_elem_here = scalar(@currm);
	for $i (0 .. $num_elem_here-1) {
		for $j (0 .. $num_elem_prev-1) {
			$t1 = $cnt1s[$comb_lvl][$j];
			$t2 = $currc[$i];
			$m1 = $mask[$comb_lvl][$j];
			$m2 = $currm[$i];
			$m = $m1 | $m2;
			next if $m1 & $m2;
			next if ($t1+$t2 == $num_elem);
			next if defined $allkey{$m};
			#next if ( ($comb_lvl == 0) && ($t1+$t2 != ($num_elem-1)) );
			++$end_cnt if ($t1+$t2 == $num_elem-1);
			#printf "\$m%0${hit}x = min( \$m%0${hit}x, \$m%0${hit}x );\n",$m,$m1,$m2;
			if ( $fn_type eq 'prod' ) {
				printf "\tSIMD_CLASS m%0${hit}x = m%0${hit}x * m%0${hit}x;\n",$m,$m1,$m2;
			} else {
				printf "\tSIMD_CLASS m%0${hit}x = min( m%0${hit}x, m%0${hit}x );\n",$m,$m1,$m2;
			}
			$allkey{$m} = 1;
			$key_thispass{$m2} = 1; #remove this from list
			if ( $t1+$t2 != $num_elem-1 ) {
				push(@nextm,$m1 | $m2);
				push(@nextc,$t1+$t2);
			}
		}
	}
	#need to copy any unused elements in curr level
	for $i (0 .. $#currm) {
		$del = 0;
		foreach $key (keys(%key_thispass)) {
			if ( $currm[$i] == $key ) {
				$del = 1;
			}
		}
		if ( $del == 0 ) {
			push(@nextm,$currm[$i]);
			push(@nextc,$currc[$i]);
		}
	}
	--$comb_lvl;
	#print "lvl=$lvl  pass1=$pass1 $end_cnt\n";
	@currm = @nextm;
	@currc = @nextc;
	$not_done = 1*( $end_cnt < $num_elem );
	#print "comb_lvl = $comb_lvl\n";
	$not_done = 0 if ( $comb_lvl < -1 );
} while ($not_done == 1);
print $res;
