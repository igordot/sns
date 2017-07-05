#!/usr/bin/env perl

use strict;

my $HELP = <<HELP;

  Convert VCF file to tab-delimited table for merging with ANNOVAR output.

  usage: vcf-table.pl file.vcf sample_name > vcf_table.txt

HELP

if (!$ARGV[1]) {
	die $HELP;
}

main();

# main subroutine
sub main {
	my $vcf = $ARGV[0];
	my $sample = $ARGV[1];

	# get variant caller (sample can't be determined from some VCFs)
	# my ($caller, $sample) = get_info($vcf);
	my $caller = get_info($vcf);

	# print germline header
	if ($caller eq "HaplotypeCaller" || $caller eq "LoFreq") {
		print "#MUT\tSAMPLE\tCHR\tPOS\tQUAL\tDEPTH\tFREQ\n";
	}
	# print old MuTect somatic header
	if ($caller eq "MuTect114" || $caller eq "MuTect116") {
		print "#MUT\tSAMPLE N\tSAMPLE T\tCHR\tPOS\tt_lod_fstar ( log of ( likelihood event is real / likelihood event is seq error ) )\tN Depth\tN Freq\tT Depth\tT Freq\n";
	}
	# print somatic header (T and N are separate columns to make it clear which is which)
	if ($caller eq "MuTect2" || $caller eq "Strelka") {
		print "#MUT\tSAMPLE T\tSAMPLE N\tCHR\tPOS\tQUAL\tT DEPTH\tT FREQ\tN DEPTH\tN FREQ\n";
	}

	# remove header lines
	my @vcf_records = `cat $vcf | grep -v "^#" | LC_ALL=C sort -k1,1 -k2,2n`;
	while (my $line = shift(@vcf_records)) {
		chomp($line);
		my $clean_line;

		if ($caller eq "HaplotypeCaller") {
			$clean_line = format_haplotypecaller($line, $sample);
		}
		if ($caller eq "LoFreq") {
			$clean_line = format_lofreq($line, $sample);
		}
		if ($caller eq "MuTect114") {
			# $clean_line = format_mutect_114($line);
		}
		if ($caller eq "MuTect116") {
			# $clean_line = format_mutect_116($line);
		}
		if ($caller eq "MuTect2") {
			$clean_line = format_mutect2($line, $sample);
		}
		if ($caller eq "Strelka") {
			$clean_line = format_strelka($line, $sample);
		}


		# skip lines that returned empty
		if ($clean_line) {
			print $clean_line . "\n";
		}
	}

}

# get variant caller
sub get_info {
	my $vcf = $_[0];

	my $caller;

	open IN, $vcf or die $!;
	foreach my $line (<IN>) {
		chomp($line);

		# get caller
		if ($line =~ /ID=HaplotypeCaller/) {
			$caller = "HaplotypeCaller";
			last;
		}
		if ($line =~ /source=lofreq/) {
			$caller = "LoFreq";
			last;
		}
		if ($line =~ /muTector/) {
			# mutect actually uses the txt output instead of vcf (more data and better formatting)
			$caller = "MuTect114";
			last;
		}
		if ($line =~ /MuTect:1.1.6/) {
			# mutect actually uses the txt output instead of vcf (more data and better formatting)
			$caller = "MuTect116";
			last;
		}
		if ($line =~ /GATKCommandLine.MuTect2/) {
			$caller = "MuTect2";
			last;
		}
		if ($line =~ /source=strelka/) {
			$caller = "Strelka";
			last;
		}

	}
	close IN;

	return ($caller);
}

# format HaplotypeCaller line
sub format_haplotypecaller {
	my $line = $_[0];
	my $sample = $_[1];
	my $out;

	# split the needed columns (#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE])
	my @cols = split(/\t/, $line);
	# sample fields (GT:AD:DP:GQ:PL)
	my @sample_fields = split(':', $cols[9]);
	# allelic depth field (Allelic depths for the ref and alt alleles in the order listed)
	my @ad_cols = split(',', $sample_fields[1]);

	# output columns
	my $chr = $cols[0];
	my $pos = $cols[1];
	my $ref = $cols[3];
	my $alt = $cols[4];
	my $depth = $sample_fields[2];

	# do not report if
	# there are not any quality reads (reported depth 0)
	# less than 5 variant call supporting reads
	if ( ($depth > 0) && ($ad_cols[1] >= 5) ) {
		($pos, $ref, $alt) = adjust_indels($pos, $ref, $alt);
		my $qual = (sprintf "%.1f", $cols[5]);
		my $freq = sprintf("%.3f", ( $ad_cols[1] / $depth ));
		$out = join("\t", "${chr}:${pos}:${ref}:${alt}", $sample, $chr, $pos, $qual, $depth, $freq);
	}

	$out;
}

# format lofreq line
sub format_lofreq {
	my $line = $_[0];
	my $sample = $_[1];
	my $out;

	# split the needed columns (#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE])
	my @cols = split(/\t/, $line);

	# output columns
	my $chr = $cols[0];
	my $pos = $cols[1];
	my $ref = $cols[3];
	my $alt = $cols[4];
	my $qual = $cols[5];
	my $info = $cols[7];
	my $depth = $info;
	$depth =~ s/.*DP=(.+?);.*/$1/;
	my $freq = $info;
	$freq =~ s/.*AF=(.+?);.*/$1/;

	# do not report if frequency is less than 1%
	if ($freq > 0.01) {
		($pos, $ref, $alt) = adjust_indels($pos, $ref, $alt);
		$freq = sprintf("%.3f", $freq);
		$out = join("\t", "${chr}:${pos}:${ref}:${alt}", $sample, $chr, $pos, $qual, $depth, $freq);
	}

	$out;
}

# format MuTect 1.1.4 line
sub format_mutect_114 {
	my $line = $_[0];
	my $out;

	# split the needed columns (#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE])
	my @cols = split(/\t/, $line);

	# output columns
	my $chr = $cols[0];
	my $pos = $cols[1];
	my $ref = $cols[3];
	my $alt = $cols[4];
	my $sample_t = $cols[5];
	my $sample_n = $cols[6];
	my $t_lod_fstar = sprintf("%.1f", $cols[16]);
	my $t_ref_count = $cols[20];
	my $t_alt_count = $cols[21];
	my $n_ref_count = $cols[30];
	my $n_alt_count = $cols[31];

	my $depth_t = $t_alt_count + $t_ref_count;
	my $freq_t = 0;
	if ($depth_t > 0) {
		$freq_t = sprintf("%.3f", ( $t_alt_count / $depth_t ));
	}

	my $depth_n = $n_alt_count + $n_ref_count;
	my $freq_n = 0;
	if ($depth_n > 0) {
		$freq_n = sprintf("%.3f", ( $n_alt_count / $depth_n ));
	}

	my $judgement = $cols[34];

	if ($judgement eq "KEEP") {
		$out = join("\t", "${chr}:${pos}:${ref}:${alt}", $sample_n, $sample_t, $chr, $pos, $t_lod_fstar, $depth_n, $freq_n, $depth_t, $freq_t);
	}

	$out;
}

# format MuTect 1.1.6/1.1.7 line
sub format_mutect_116 {
	my $line = $_[0];
	my $out;

	# split the needed columns (#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE])
	my @cols = split(/\t/, $line);

	# output columns
	my $chr = $cols[0];
	my $pos = $cols[1];
	my $ref = $cols[3];
	my $alt = $cols[4];
	my $sample_t = $cols[5];
	my $sample_n = $cols[6];
	my $t_lod_fstar = sprintf("%.1f", $cols[18]);
	my $t_q20_count = $cols[24];
	my $t_ref_count = $cols[25];
	my $t_alt_count = $cols[26];
	my $tumor_f = sprintf("%.3f", $cols[21]);
	my $n_q20_count = $cols[36];
	my $n_ref_count = $cols[37];
	my $n_alt_count = $cols[38];
	my $normal_f = sprintf("%.3f", $cols[35]);

	# consider replacing raw depth with q20 depth (q20 depth and frequencies were not available in 1.1.4)

	my $depth_t = $t_alt_count + $t_ref_count;
	#my $freq_t = 0;
	#if ($depth_t > 0) {
	#	$freq_t = sprintf("%.3f", ( $t_alt_count / $depth_t ));
	#}

	my $depth_n = $n_alt_count + $n_ref_count;
	#my $freq_n = 0;
	#if ($depth_n > 0) {
	#	$freq_n = sprintf("%.3f", ( $n_alt_count / $depth_n ));
	#}

	my $judgement = $cols[50];

	if ($judgement eq "KEEP") {
		#$out = join("\t", "${chr}:${pos}:${ref}:${alt}", $sample_n, $sample_t, $chr, $pos, $t_lod_fstar, $depth_n, $freq_n, $depth_t, $freq_t);
		$out = join("\t", "${chr}:${pos}:${ref}:${alt}", $sample_n, $sample_t, $chr, $pos, $t_lod_fstar, $depth_n, $normal_f, $depth_t, $tumor_f);
	}

	$out;
}

# format MuTect2 line
sub format_mutect2 {
	my $line = $_[0];
	my $sample = $_[1];
	my $out;

	my ($sample_t, $sample_n) = split(':', $sample);

	# split the columns (#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT TUMOR NORMAL)
	my @cols = split(/\t/, $line);

	# T/N sample fields (GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:PGT:PID:QSS:REF_F1R2:REF_F2R1)
	my @sample_t_fields = split(':', $cols[9]);
	my @sample_n_fields = split(':', $cols[10]);

	# allelic depth field (Allelic depths for the ref and alt alleles in the order listed)
	# not using AF for now because number format is scientific for low fractions
	my @t_ad_cols = split(',', $sample_t_fields[1]);
	my @n_ad_cols = split(',', $sample_n_fields[1]);

	# output columns
	my $chr = $cols[0];
	my $pos = $cols[1];
	my $ref = $cols[3];
	my $alt = $cols[4];
	my $qual = $cols[7];
	$qual =~ s/.*TLOD=([0-9.]+).*/$1/;
	my $t_depth = $t_ad_cols[0] + $t_ad_cols[1];
	my $n_depth = $n_ad_cols[0] + $n_ad_cols[1];
	my $t_freq = $t_ad_cols[1] / $t_depth;
	my $n_freq = $n_ad_cols[1] / $n_depth;

	# report if
	# T frequency is more than 3%
	# N frequency is less than 5%
	# at least 5 variant call supporting reads
	# T frequency is sufficiently higher than N frequency
	# "we recommend applying post-processing filters, e.g. by hard-filtering calls with low minor allele frequencies"
	if ( ($t_freq > 0.03) && ($n_freq < 0.05) && ($t_ad_cols[1] >= 5) && ($t_freq > $n_freq * 5) ) {
		($pos, $ref, $alt) = adjust_indels($pos, $ref, $alt);
		$t_freq = sprintf("%.3f", $t_freq);
		$n_freq = sprintf("%.3f", $n_freq);
		my $mut_id = "${chr}:${pos}:${ref}:${alt}";
		$out = join("\t", $mut_id, $sample_t, $sample_n, $chr, $pos, $qual, $t_depth, $t_freq, $n_depth, $n_freq);
	}

	$out;
}

# format Strelka line
sub format_strelka {
	my $line = $_[0];
	my $sample = $_[1];
	my $out;

	my ($sample_t, $sample_n) = split(':', $sample);

	# split the columns (#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NORMAL TUMOR)
	my @cols = split(/\t/, $line);

	# T/N sample fields (DP:FDP:SDP:SUBDP:AU:CU:GU:TU)
	my @sample_t_fields = split(':', $cols[10]);
	my @sample_n_fields = split(':', $cols[9]);

	# output columns
	my $chr = $cols[0];
	my $pos = $cols[1];
	my $ref = $cols[3];
	my $alt = $cols[4];
	my $qual = $cols[7];
	my $t_depth = 0;
	my $n_depth = 0;
	my $t_freq = 1;
	my $n_freq = 1;
	my $t_ref_count = 0;
	my $t_alt_count = 999;
	my $n_ref_count = 999;
	my $n_alt_count = 0;

	# SNV lines
	if ($line =~ /AU:CU:GU:TU/) {
		# allele fields (AU:CU:GU:TU)
		my %t_alleles = (
			A => (split(',', $sample_t_fields[4]))[0],
			C => (split(',', $sample_t_fields[5]))[0],
			G => (split(',', $sample_t_fields[6]))[0],
			T => (split(',', $sample_t_fields[7]))[0]
		);
		my %n_alleles = (
			A => (split(',', $sample_n_fields[4]))[0],
			C => (split(',', $sample_n_fields[5]))[0],
			G => (split(',', $sample_n_fields[6]))[0],
			T => (split(',', $sample_n_fields[7]))[0]
		);

		# output columns
		$qual =~ s/.*;QSS=(.+?);.*/$1/;
		$t_alt_count = $t_alleles{$alt};
		$n_alt_count = $n_alleles{$alt};
		$t_depth = $sample_t_fields[0];
		$n_depth = $sample_n_fields[0];
		$t_freq = $t_alt_count / $t_depth;
		$n_freq = $n_alt_count / $n_depth;
	}

	# indel lines
	if ($line =~ /TAR:TIR:TOR/) {
		# from https://github.com/Illumina/strelka/issues/3
		# somatic indel allele frequency is: alt_t1count / (ref_t1count + alt_t1count)
		# ref_t1count = 1st value of FORMAT column value TAR
		# alt_t1count = 1st value of FORMAT column value TIR

		# extract ref/alt allele counts from TAR/TIR fields
		$t_ref_count = (split(',', $sample_t_fields[2]))[0];
		$t_alt_count = (split(',', $sample_t_fields[3]))[0];
		$n_ref_count = (split(',', $sample_n_fields[2]))[0];
		$n_alt_count = (split(',', $sample_n_fields[3]))[0];

		# output columns
		$qual =~ s/.*;QSI=(.+?);.*/$1/;
		$t_depth = $t_ref_count + $t_alt_count;
		$n_depth = $n_ref_count + $n_alt_count;
		$t_freq = $t_alt_count / $t_depth;
		$n_freq = $n_alt_count / $n_depth;
	}

	# report if
	# T frequency is more than 3%
	# at least 5 variant call supporting reads
	if ( ($t_freq > 0.03) && ($t_alt_count >= 5) ) {
		($pos, $ref, $alt) = adjust_indels($pos, $ref, $alt);
		$t_freq = sprintf("%.3f", $t_freq);
		$n_freq = sprintf("%.3f", $n_freq);
		my $mut_id = "${chr}:${pos}:${ref}:${alt}";
		$out = join("\t", $mut_id, $sample_t, $sample_n, $chr, $pos, $qual, $t_depth, $t_freq, $n_depth, $n_freq);
	}

	$out;
}

# adjust indels to match annovar output
sub adjust_indels {
	my $pos = $_[0];
	my $ref = $_[1];
	my $alt = $_[2];

	# deletion (pos has to be incremented to match annovar)
	if (length($ref) > length($alt)) {
		$pos = $pos + 1;
		$ref = substr($ref, 1);
		$alt = "-";
	}

	# insertion
	if (length($alt) > length($ref)) {
		$ref = "-";
		$alt = substr($alt, 1);
	}

	return ($pos, $ref, $alt);
}



# end
