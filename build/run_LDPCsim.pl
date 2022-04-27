#!/usr/bin/env -S perl -w

use IPC::Open2;
use IO::Select;
use Time::HiRes qw( gettimeofday tv_interval );

#default values - can be overridden by command line evals
$exec = 'LDPCtop_Vec64c_min';
$MatrixSet = '11ad 1/2';
$num_par = 6; #max is 9
$simd_size = 512; #SIMD word size for Decode words (256 or 512)
$elem_size = 8; #Element size in bits (32 for float, 64 for double, 8/16/32 for int)
$gen_dat_flag = 1; #generate encode and decode .dat files (else use existing ones)
#$accum_flag = 1; #TBD.  Want a flag to determine if all threads get same params except seed, and accumulate results
                   #or, if each thread has some different params.  And, maybe large number of unique params
#Realistically, $use_given_rand_state=1 implies $accum_flag is 0, no need for extra flag

($num_proc) = `grep processor /proc/cpuinfo | tail -1` =~ /: (\d+)/;
++$num_proc;  #default number of processes for decoding

#@min_fer  = (1e5, 1000, 20);
#@time_fer = (30, 10*60, 24*3600); #last time is max wall clock - exit after report exceeding this

$outer_count  = 10; #reports made on each outer iteration
$inner_count  = 100; #100000; #number of times each block of "num_codeword" is decoded
$num_codeword = 1; #These are "super" codewords including parallelism in one chunk of memory
$max_iter = 10; #50;

$max_row_metric = 31;
$max_abs_LLR = 24;
$LLR_scaling_factor = 2.0;
$beta_offset = 1.0;

$puncture_8023ca_flag = 0;
$use_given_rand_state = 0;
$rand_seed =   100; #only relevant if use_given_rand_state is 0
$SNR_start = -1.00;
$SNR_end   =  0.00;
$SNR_step  =  0.25;
$enc_raw_file = "LDPCencode.raw.dat"; #could rename this file after 1st time, put new name here, and set gen_dat_flag=0
$dec_raw_file = "LDPCdecode.raw.dat";

$save_stats  = 1;  #save stats to an Octave file
$mainid = 'm11ad12';  #some prefix for variables in Octave file
#grabbing TSC frequency this way is likely not robust.
@tsc_time = `dmesg | grep "tsc: Detected"`; #grabs first line, with "processor"
($tick_freq) = $tsc_time[0] =~ m/(\S+)\s*MHz/;

while ($_ = shift) {
	s/::/\"/g;
	eval("\$$_");	#blindly evaluating arg - use with caution
}

#default metrics for floating point (unless overridden)
#if manually setting, must set all of max_row_metric, max_abs_LLR, LLR_scaling_factor, and beta_offset
if ( ($exec =~ /_prod/) && ($max_row_metric == 31) ) { 
	#max_row_metric is max 1-tanh value if using sum-product
	#6e-8 if float, 1e-16 if double (corresponds to ~18/2 or ~38/2 as tanh argument)
	$max_row_metric = ($exec =~ /f_prod/) ? 6e-8 : 1e-16;
	$max_abs_LLR = 1e9;
	$LLR_scaling_factor = 1.0;
	$beta_offset = 0.0;
}
#compute simd_size and elem_size based on exec name
%elem_tab = ('c',8,'s',16,'w',32,'q',64,'f',32,'d',64);
($simd_par,$echar) = $exec =~ /Vec(\d+)([cswqfd])/;
$elem_size = $elem_tab{$echar};
$simd_size = $simd_par * $elem_size;
$simd_size = $simd_size_over if defined $simd_size_over;  #override simd_size
$elem_size = $elem_size_over if defined $elem_size_over; #override elem_size
print "Key Parameters: MatrixSet=\"$MatrixSet\" num_par=$num_par Decoder=$exec simd_size=$simd_size elem_size=$elem_size num_proc=$num_proc num_codeword=$num_codeword\n\n";
$target_loop_cnt = int( ($SNR_end-$SNR_start)/$SNR_step ) + 1;

srand($rand_seed);
for $i (0..63) {
	$randval[$i] = 2**56 * rand;
}

$int_vals  = pack('L6',$outer_count,$inner_count,$num_codeword,$max_iter,$puncture_8023ca_flag,$use_given_rand_state);
$flt_vals = pack('f7',$max_row_metric,$LLR_scaling_factor,$max_abs_LLR,$beta_offset,$SNR_start,$SNR_end,$SNR_step);

$sim_vals = $int_vals . $flt_vals . pack('Q64',@randval) . pack('Z128',$enc_raw_file) . pack('Z128',$dec_raw_file);
if (1) { #this is for running the decoder in debug mode with the saved simulation values
	open(SIM,">sim_vals.bin") || die;
	print SIM $sim_vals;
	close SIM
}

$matlab_exe = ($ENV{'USEOCT'} == 1) ? 'octave --eval' : 'matlab -batch';
if ( $gen_dat_flag ) { #generate the Encoder and Decoder structure files
	$zstr = (defined $z_over) ? "z_over=$z_over" : '';
	$cmd = "$matlab_exe \"MatrixSet=\'$MatrixSet\';num_par=$num_par;simd_size=$simd_size;elem_size=$elem_size;$zstr;gen_structs\"";
	print "$cmd\n";
	system($cmd);
}

$start_tstamp = [gettimeofday];
$parse_loop_cnt = 0;
$cursor_up = 0;
$cw_rate_update = 1;

$s = IO::Select->new();
for $proc (0 .. $num_proc-1) {
	open2($fho[$proc],$fhi[$proc],"./$exec") or die "open2() failed $!";
	$s->add($fho[$proc]);
	$fh = $fhi[$proc];
	print $fh "$sim_vals";
	#$fh_map{$fho[$proc]} = $proc;
}
$child_count = $num_proc;
while ($child_count > 0) {
	@ready = $s->can_read(1);
	for $i (0 .. $#ready) {
		$_ = readline($ready[$i]);
		if (!defined $_) {
			$s->remove($ready[$i]);
			if ( $child_count == $num_proc ) { #first child finished
				$etime1 = tv_interval( $start_tstamp, [gettimeofday] );
				$cw_rate_update = 0;
			}
			--$child_count;
		} else {
			#$ch_num = $fh_map{$ready[$i]};
			#print "Child $ch_num: $_";
			&parse_line($_);
		}
	}
	if ( $parse_loop_cnt >= $target_loop_cnt ) {
		$parse_loop_cnt =  0;
		&print_stats($cw_rate_update);
	}
}
$etime2 = tv_interval( $start_tstamp, [gettimeofday] );

&print_stats($cw_rate_update);
$percent_par = ( $num_proc == 1 ) ? 100 : 100*$etime1/$etime2;
$percent_dec = 100*$total_dec_time / $etime1;
printf("\nPercent time all proc=%6.1f%%    Approx. Decode percentage=%5.1f%%\n",$percent_par,$percent_dec);

if ($save_stats) { #append/create to fixed name.  Could change in future
	@cat = ('SNR','BER','iter','CWER','Dec_time','cwrate','inp_ber');
	open(OCT,">> ldpc_stat.m") || die;
	print OCT "accum_stat = [",join(' ',@all_stats),"];\n";
	for $i (0..6) {
		print OCT "${mainid}_$cat[$i] = accum_stat($i+1:7:end);\n";
	}
	print OCT "add_label=1;\n";
	print OCT "if exist('fig1')\nhold on;\nadd_label=0;\nlegend_str{end+1}='$mainid BER';\nlegend_str{end+1}='$mainid CWER';\nend\n";
	print OCT "fig1 = semilogy(${mainid}_SNR,${mainid}_BER,ls_tab{ls_idx},${mainid}_SNR,${mainid}_CWER,ls_tab{ls_idx+1});\n";
	print OCT "if (add_label)\nxlabel('SNR');\nylabel('BER/CWER');\ngrid on;\n";
	print OCT "title('Bit and Codeword Error Rate vs SNR');\nlegend_str={'$mainid BER','$mainid CWER'};\nend\n";
	print OCT "ls_idx = ls_idx+2; legend(legend_str);\n\n";
}

#parse line from decoder processes.  Uses global variables.
sub parse_line {
	my ($inp) = @_;
	$_ = $inp;
	if ( /^CW Count/ ) {
		($super_cw_per_loop) = m/^CW Count: (\d+)/;
	}
	if ( /^Matrix Params/ ) {
		($Z,$Hcol,$Hrow,$num_par_cw,$fill_cnt) = m/Z=(\d+)  Hcol=(\d+)  Hrow=(\d+) num_par_cw=(\d+) fill_cnt=(\d+)/;
		$cw_per_loop = $super_cw_per_loop * $num_par_cw;
		$bit_per_scw = $Z * $Hcol;
		$data_bit_per_scw = $Z * ($Hcol - $Hrow) - $fill_cnt;
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
		++$parse_loop_cnt;
	}
}

sub print_stats {
	my ($cw_rate_flag) = @_;
	if ( $cursor_up ) {
		$tmp = "\e[" . scalar(keys %in_be_tot) . "A";
		print $tmp;
	}
	$total_dec_time = 0;
	@all_stats = ();
	if ( $cw_rate_flag ) {
		@cw_rate_save = ();
	}
	foreach $snr (sort { $a <=> $b } keys %in_be_tot) {
		$tot_scw = $tot_cw{$snr}/$num_par; #total "super" codewords (includes num_par)
		$tot_out_ber = $out_be_tot{$snr} / ($tot_scw * $data_bit_per_scw);
		$tot_in_ber = $in_be_tot{$snr} / ($tot_scw * $bit_per_scw);
		$ave_iter = $iter_tot{$snr} / $tot_scw;
		$tot_cw_er = $cw_err_tot{$snr} / $tot_cw{$snr};
		$dec_time = $tot_ticks{$snr} / ($num_proc * $tick_freq * 1e6);
		if ( $cw_rate_flag ) {
			$cw_per_sec = $tot_cw{$snr} / $dec_time;
		} else {
			$cw_per_sec = shift(@cw_rate_save);
		}
		push(@cw_rate_save,$cw_per_sec);
		$total_dec_time += $dec_time;
		push(@all_stats,$snr,$tot_out_ber,$ave_iter,$tot_cw_er,$dec_time,$cw_per_sec,$tot_in_ber);
		printf("SNR=%6.2f   BER=%8.2e   CWER=%8.2e   cw dec/sec=%8.2e  tot_cw=%9d\n",$snr,$tot_out_ber,$tot_cw_er,$cw_per_sec,$tot_cw{$snr});
	}
	$cursor_up = 1;
}
