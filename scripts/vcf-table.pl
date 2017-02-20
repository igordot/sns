#!/usr/bin/env perl

use strict;

my $HELP = <<HELP;

  Convert VCF file to tab-delimited table for merging with ANNOVAR output.

  usage: vcf-format.pl file.vcf sample_name > vcf_table.txt

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

	# print header (for somatic, T and N are separate columns to make it clear which is which)
	if ($caller eq "HaplotypeCaller" || $caller eq "LoFreq") {
		print "#MUT\tSAMPLE\tCHR\tPOS\tQUAL\tDEPTH\tFREQ\n";
	}
	if ($caller eq "MuTect114" || $caller eq "MuTect116") {
		print "#MUT\tSAMPLE N\tSAMPLE T\tCHR\tPOS\tt_lod_fstar ( log of ( likelihood event is real / likelihood event is seq error ) )\tN Depth\tN Freq\tT Depth\tT Freq\n";
	}
	if ($caller eq "MuTect2") {
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
		($pos, $ref, $alt) = fix_indels($pos, $ref, $alt);
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
		($pos, $ref, $alt) = fix_indels($pos, $ref, $alt);
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

	# split the needed columns (#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT TUMOR NORMAL)
	my @cols = split(/\t/, $line);
	# sample fields (GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:PGT:PID:QSS:REF_F1R2:REF_F2R1)
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

	# do not report if
	# T frequency is less than 3%
	# N frequency is more than 5%
	# less than 5 variant call supporting reads
	# T frequency is sufficiently higher than N frequency
	# "we recommend applying post-processing filters, e.g. by hard-filtering calls with low minor allele frequencies"
	if ( ($t_freq > 0.03) && ($n_freq < 0.05) && ($t_ad_cols[1] >= 5) && ($t_freq > $n_freq * 5) ) {
		($pos, $ref, $alt) = fix_indels($pos, $ref, $alt);
		$t_freq = sprintf("%.3f", $t_freq);
		$n_freq = sprintf("%.3f", $n_freq);
		$out = join("\t", "${chr}:${pos}:${ref}:${alt}", $sample_t, $sample_n, $chr, $pos, $qual, $t_depth, $t_freq, $n_depth, $n_freq);
	}

	$out;
}

# fix indels to match annovar output
sub fix_indels {
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
