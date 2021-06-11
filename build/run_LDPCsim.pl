#!/usr/bin/env -S perl -w

use IPC::Open2;
use IO::Select;

#default values - can be overridden by command line evals
$exec = 'LDPCtop_Vec16f_min'; #'LDPCtop_Vec64c_min';
$MatrixSet = '11ad 1/2';
$num_par = 1; #max is 9
$simd_size = 512; #SIMD word size for Decode words (256 or 512)
$elem_size = 32; #Element size in bits (32 for float, 64 for double, 8/16/32 for int)
$gen_dat_flag = 1; #generate encode and decode .dat files (else use existing ones)
#$accum_flag = 1; #TBD.  Want a flag to determine if all threads get same params except seed, and accumulate results
                   #or, if each thread has some different params.  And, maybe large number of unique params
#Realistically, $use_given_rand_state=1 implies $accum_flag is 0, no need for extra flag

($num_proc) = `grep processor /proc/cpuinfo | tail -1` =~ /: (\d+)/;
++$num_proc;  #default number of processes for decoding
$num_proc = 1; #override for now.  Remove this later

#@min_fer  = (1e5, 1000, 20);
#@time_fer = (30, 10*60, 24*3600); #last time is max wall clock - exit after report exceeding this

$outer_count = 10; #reports made on each outer iteration
$inner_count = 25; #number of times each block of "num_codeword" is decoded
$num_codeword = 500; #These are "super" codewords including parallelism in one chunk of memory
$max_iter = 10; #50;

#max_row_metric is max 1-tanh value if using sum-product
#6e-8 if float, 1e-16 if double (corresponds to ~18/2 or ~38/2 as tanh argument)
if ($exec =~ /_prod/) {
	$max_row_metric = ($exec =~ /f_prod/) ? 6e-8 : 1e-16;
	$max_abs_LLR = 1e9;
	$LLR_scaling_factor = 1.0;
	$beta_offset = 1.0;
} else {
	$max_row_metric = 31;
	$max_abs_LLR = 24;
	$LLR_scaling_factor = 2.0;
	$beta_offset = 1.0;
}

$puncture_8023ca_flag = 0;
$use_given_rand_state = 0;
$rand_seed = 100;
$SNR_start = 3.50; #-1.00;
$SNR_end   = 3.70; # 0.00;
$SNR_step  = 0.05; #0.25;
$enc_raw_file = "LDPCencode.raw.dat"; #"LDPCencode.11ad12.1x.raw.dat"; #11n.z27.raw.dat"; #16e12.z28.raw.dat";
$dec_raw_file = "LDPCdecode.raw.dat"; #"LDPCdecode.11ad12.512e8.raw.dat"; #n12.z27.256e32.raw.dat"; #16e12.256e32.raw.dat";

while ($_ = shift) {
	s/::/\"/g;
	eval("\$$_");	#blindly evaluating arg - use with caution
}

srand($rand_seed);
for $i (0..63) {
	$randval[$i] = 2**56 * rand;
}

$int_vals  = pack('L6',$outer_count,$inner_count,$num_codeword,$max_iter,$puncture_8023ca_flag,$use_given_rand_state);
$flt_vals = pack('f7',$max_row_metric,$LLR_scaling_factor,$max_abs_LLR,$beta_offset,$SNR_start,$SNR_end,$SNR_step);

$sim_vals = $int_vals . $flt_vals . pack('Q64',@randval) . pack('Z128',$enc_raw_file) . pack('Z128',$dec_raw_file);

$matlab_exe = ($ENV{'USEOCT'} == 1) ? 'octave --eval' : 'matlab -batch';
if ( $gen_dat_flag ) { #generate the Encoder and Decoder structure files
	$zstr = (defined $z_over) ? "z_over=$z_over" : '';
	$cmd = "$matlab_exe \"MatrixSet=\'$MatrixSet\';num_par=$num_par;simd_size=$simd_size;elem_size=$elem_size;$zstr;gen_structs\"";
	system($cmd);
}

$s = IO::Select->new();
for $proc (0 .. $num_proc-1) {
	open2($fho[$proc],$fhi[$proc],"./$exec") or die "open2() failed $!";
	$s->add($fho[$proc]);
	$fh = $fhi[$proc];
	print $fh "$sim_vals";
	$fh_map{$fho[$proc]} = $proc;
}
$child_count = $num_proc;
while ($child_count > 0) {
	@ready = $s->can_read(1);
	for $i (0 .. $#ready) {
		$_ = readline($ready[$i]);
		if (!defined $_) {
			$s->remove($ready[$i]);
			--$child_count;
		} else {
			$ch_num = $fh_map{$ready[$i]};
			#print "Child $ch_num: $_";
			#($SNR,$inp_err,$itercnt,$biterr,$cwerr,$ticks) =~ /SNR=(\S+)  inp_bit_err=(\d+) tot_iter=(\d+) tot_cw_err=(\d+) tot_dec_ticks=(\d+)/;
			&parse_line($_);
		}
	}
}

#grabbing TSC frequency this way is not universal method
@tsc_time = `egrep "model name" /proc/cpuinfo | head -1`;
($tick_freq) = $tsc_time[0] =~ m/(\S+)GHz/;

#cw_per_sec assumes that all processes run in parallel (some will end slightly before others, so not 100% accurate)
print "% value order: snr, in_ber, out_ber, cwer, tot_cw_with_error, ave_iterations, cw_decodes_per_sec, tot_cw\n";
print "sum_vals = [\n";
foreach $snr (sort { $a <=> $b } keys %in_be_tot) {
    $tot_in_ber = $in_be_tot{$snr} / ($tot_cw{$snr}*$bit_per_cw);
    $tot_out_ber = $out_be_tot{$snr} / ($tot_cw{$snr} * $data_bit_per_cw);
    $ave_iter = $iter_tot{$snr} * $num_par / $tot_cw{$snr};
    $tot_cw_er = $cw_err_tot{$snr} / $tot_cw{$snr};
    $cw_per_sec = $tot_cw{$snr} / ($tot_ticks{$snr} / ($num_proc * $tick_freq * 1e9));

    printf "    $snr $tot_in_ber $tot_out_ber $tot_cw_er $cw_err_tot{$snr} $ave_iter $cw_per_sec $tot_cw{$snr};\n";
}

sub parse_line {
	my ($inp) = @_;
	$_ = $inp;
	if ( /^CW Count/ ) {
		($cw_per_loop) = m/^CW Count: (\d+)/;
	}
	if ( /^Matrix Params/ ) {
		($Z,$Hcol,$Hrow,$num_par_cw) = m/Z=(\d+)  Hcol=(\d+)  Hrow=(\d+) num_par_cw=(\d+)/;
		$cw_per_loop = $cw_per_loop * $num_par_cw;
		$bit_per_cw = $Z * $Hcol;
		$data_bit_per_cw = $Z * ($Hcol - $Hrow);
	}
	if ( /^loop/ ) {
		@vals = m/SNR=(\S+)\s+inp_bit_err=(\S+)\s+tot_iter=(\S+)\s+tot_bit_err=(\S+)\s+tot_cw_err=(\S+)\s+tot_dec_ticks=(\S+)/;
		($SNR,$in_be,$iter,$out_be,$cw_err,$ticks) = @vals;
		$in_be_tot{$SNR} += $in_be;
		$out_be_tot{$SNR} += $out_be;
		$iter_tot{$SNR} += $iter;
		$cw_err_tot{$SNR} += $cw_err;
		$tot_ticks{$SNR} += $ticks;
		$tot_cw{$SNR} += $cw_per_loop;
	}
}
