#!/usr/bin/env Rscript


##
## Parse a VCF file and generate a table compatible with ANNOVAR output.
##
## usage: Rscript --vanilla vcf-table.R sample_name in.vcf out.txt
##


# increase output width
options(width = 120)
# print warnings as they occur
options(warn = 1)

# get scripts directory (directory of this file) and load relevant functions
args_all = commandArgs(trailingOnly = FALSE)
scripts_dir = normalizePath(dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
source(paste0(scripts_dir, "/load-install-packages.R"))

# relevant arguments
args = commandArgs(trailingOnly = TRUE)
sample_name = args[1]
in_vcf = args[2]
out_txt = args[3]

# check that inputs are valid
if (length(args) != 3) stop("usage: Rscript --vanilla vcf-table.R sample_name in.vcf out.txt")
if (!file.exists(in_vcf)) stop("file does not exist: ", in_vcf)

# load relevant packages
load_install_packages("dplyr")
load_install_packages("readr")
load_install_packages("tidyr")
load_install_packages("stringr")
load_install_packages("glue")
load_install_packages("vcfR")

# GATK HaplotypeCaller
parse_gatkhc = function(vcfr_obj, sample_name) {

  ##INFO DP Approximate read depth; some reads may have been filtered
  ##FORMAT AD Allelic depths for the ref and alt alleles in the order listed
  ##FORMAT DP Approximate read depth (reads with MQ=255 or with bad mates are filtered)
  ##FORMAT GQ Genotype Quality
  ##FORMAT GT Genotype

  # convert vcfR object to a tibble
  muts_tbl = vcfR2tidy(vcfr_obj,
                       info_fields = c("DP"),
                       format_fields = c("AD", "DP"),
                       single_frame = TRUE, verbose = FALSE)
  muts_tbl = muts_tbl$dat

  # extract relevant metrics
  muts_tbl =
    muts_tbl %>%
    mutate(mut_id = glue("{CHROM}:{POS}:{REF}:{ALT}")) %>%
    separate(gt_AD, into = c("ref_counts", "alt_counts"), sep = ",", convert = TRUE, extra = "drop") %>%
    mutate(AF = alt_counts / (ref_counts + alt_counts)) %>%
    mutate(
      QUAL = round(as.numeric(QUAL), 1),
      DEPTH = gt_DP,
      FREQ = round(as.numeric(AF), 3)
    )

  # adjust indels to be in ANNOVAR format
  muts_tbl = muts_tbl %>% adjust_indels()

  # manual filtering
  muts_tbl %>% filter(DEPTH >= 10 & alt_counts >= 5 & FREQ > 0.01)

}

# LoFreq (with FORMAT and SAMPLE fields added using lofreq2_add_fake_gt.py)
parse_lofreq = function(vcfr_obj, sample_name) {

  # convert vcfR object to a tibble
  muts_tbl = vcfR2tidy(vcfr_obj,
                       info_fields = c("DP", "AF", "DP4"),
                       format_fields = c("GT"),
                       single_frame = TRUE, verbose = FALSE)
  muts_tbl = muts_tbl$dat

  # extract relevant metrics
  muts_tbl =
    muts_tbl %>%
    mutate(mut_id = glue("{CHROM}:{POS}:{REF}:{ALT}")) %>%
    separate(DP4,
             into = c("ref_fwd_counts", "ref_rev_counts", "alt_fwd_counts", "alt_rev_counts"),
             sep = ",", convert = TRUE, extra = "drop") %>%
    mutate(alt_counts = alt_fwd_counts + alt_rev_counts) %>%
    mutate(
      QUAL = round(as.numeric(QUAL), 1),
      DEPTH = DP,
      FREQ = round(as.numeric(AF), 3)
    )

  # adjust indels to be in ANNOVAR format
  muts_tbl = muts_tbl %>% adjust_indels()

  # manual filtering
  muts_tbl %>% filter(DEPTH >= 10 & alt_counts >= 5 & FREQ > 0.01)

}

# Mutect 2.1 (GATK 4.0) and Mutect 2.2 (GATK 4.1)
parse_mutect2 = function(vcfr_obj, sample_T, sample_N) {

  # confirm mutect version (header changed from "Mutect Version" to "MutectVersion" in GATK 4.0.9.0)
  if (!any(str_detect(vcfr_obj@meta, "##MutectVersion=2."))) {
    stop("version mismatch: expecting Mutect 2.X")
  }

  # confirm sample names
  if (!any(str_detect(vcfr_obj@meta, glue("##tumor_sample={sample_T}")))) {
    stop("sample mismatch")
  }
  if (!any(str_detect(vcfr_obj@meta, glue("##normal_sample={sample_N}")))) {
    stop("sample mismatch")
  }

  # if a read is considered uninformative, it is counted towards the DP, but not the AD
  # an uninformative read is not reported in the AD, it is still used in calculations for genotyping
  # if uninformative reads are the only reads, we report the potential variant allele, but keep the AD values 0
  # AD is the number of reads that more likely than not support an allele
  # if you have 10 reads, each with 0.6 probability of having a certain alt allele, you get an AD of 10, whereas you get essentially 0.6 x 10 = 6 reads for the purpose of AF

  ##INFO DP Approximate read depth; some reads may have been filtered
  ##INFO TLOD Tumor LOD score
  ##FORMAT AD Allelic depths for the ref and alt alleles in the order listed
  ##FORMAT AF Allele fractions of alternate alleles in the tumor
  ##FORMAT DP Approximate read depth (reads with MQ=255 or with bad mates are filtered

  # convert vcfR object to a tibble
  muts_tbl = vcfR2tidy(vcfr_obj, single_frame = TRUE, verbose = FALSE,
                       info_fields = c("DP", "TLOD"),
                       format_fields = c("AD", "AF", "DP"))
  muts_tbl = muts_tbl$dat

  # extract relevant metrics
  muts_tbl =
    muts_tbl %>%
    # unique mutation identifier (for joining T and N)
    mutate(mut_id = glue("{CHROM}:{POS}:{REF}:{ALT}")) %>%
    filter(FILTER == "PASS") %>%
    separate(gt_AD, into = c("ref_counts", "alt_counts"), sep = ",", convert = TRUE, extra = "drop") %>%
    # manual AF calculation for comparison
    mutate(myAF = alt_counts / (ref_counts + alt_counts)) %>%
    mutate(
      QUAL = round(as.numeric(TLOD), 1),
      T_DEPTH = ref_counts + alt_counts,
      T_FREQ = round(as.numeric(gt_AF), 3)
    )

  # extract samples T and N to put side by side ("wide" format)
  snvs_n_tbl =
    muts_tbl %>%
    filter(Indiv == sample_N) %>%
    rename(N_DEPTH = T_DEPTH, N_FREQ = T_FREQ) %>%
    select(mut_id, N_DEPTH, N_FREQ)
  muts_tbl =
    muts_tbl %>%
    filter(Indiv == sample_T) %>%
    inner_join(snvs_n_tbl, by = "mut_id")

  # adjust indels to be in ANNOVAR format
  muts_tbl = muts_tbl %>% adjust_indels()

  # manual filtering
  muts_tbl %>% filter(T_DEPTH >= 10 & N_DEPTH >= 10 & alt_counts >= 5 & T_FREQ > 0.03 & T_FREQ > (N_FREQ * 5))

}

# Strelka 2
parse_strelka2 = function(vcfr_obj, sample_T, sample_N) {

  # confirm strelka version
  if (!any(str_detect(vcfr_obj@meta, "##source_version=2"))) {
    stop("version mismatch: expecting Strelka 2")
  }

  ##INFO QSS Quality score for any somatic snv
  ##INFO SOMATIC Somatic mutation">
  ##INFO QSI Quality score for the ALT haplotype to be present at a significantly different freq in the T and N
  ##FORMAT AU Number of 'A' alleles used in tiers 1,2
  ##FORMAT CU Number of 'C' alleles used in tiers 1,2
  ##FORMAT GU Number of 'G' alleles used in tiers 1,2
  ##FORMAT TU Number of 'T' alleles used in tiers 1,2
  ##FORMAT TAR Reads strongly supporting alternate allele for tiers 1,2
  ##FORMAT TIR Reads strongly supporting indel allele for tiers 1,2

  # convert vcfR object to a tibble
  muts_tbl = vcfR2tidy(vcfr_obj,
                       info_fields = c("SOMATIC", "QSS", "QSI"),
                       format_fields = c("DP", "AU", "CU", "GU", "TU", "TAR", "TIR"),
                       single_frame = TRUE, verbose = FALSE)
  muts_tbl = muts_tbl$dat
  colnames(muts_tbl)

  # extract relevant metrics
  muts_tbl =
    muts_tbl %>%
    mutate(mut_id = glue("{CHROM}:{POS}:{REF}:{ALT}")) %>%
    filter(FILTER == "PASS") %>%
    # extract tier1 counts for each nucleotide
    separate(gt_AU, into = "A_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    separate(gt_CU, into = "C_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    separate(gt_GU, into = "G_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    separate(gt_TU, into = "T_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    separate(gt_TAR, into = "indel_ref_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    separate(gt_TIR, into = "indel_alt_counts", sep = ",", convert = TRUE, extra = "drop") %>%
    # set ref/alt counts
    mutate(
      ref_counts = case_when(
        is.na(QSS) ~ indel_ref_counts,
        REF == "A" ~ A_counts,
        REF == "C" ~ C_counts,
        REF == "G" ~ G_counts,
        REF == "T" ~ T_counts
      )
    ) %>%
    mutate(
      alt_counts = case_when(
        is.na(QSS) ~ indel_alt_counts,
        ALT == "A" ~ A_counts,
        ALT == "C" ~ C_counts,
        ALT == "G" ~ G_counts,
        ALT == "T" ~ T_counts
      )
    ) %>%
    # extract quality
    mutate(
      QUAL = case_when(
        is.na(QSI) ~ QSS,
        is.na(QSS) ~ QSI
      )
    ) %>%
    mutate(T_DEPTH = gt_DP) %>%
    mutate(T_FREQ = round(alt_counts / (ref_counts + alt_counts), 3))

  # extract samples T and N to put side by side ("wide" format)
  snvs_n_tbl =
    muts_tbl %>%
    filter(Indiv == "NORMAL") %>%
    rename(N_DEPTH = T_DEPTH, N_FREQ = T_FREQ) %>%
    select(mut_id, N_DEPTH, N_FREQ)
  muts_tbl =
    muts_tbl %>%
    filter(Indiv == "TUMOR") %>%
    inner_join(snvs_n_tbl, by = "mut_id")

  # adjust indels to be in ANNOVAR format
  muts_tbl = muts_tbl %>% adjust_indels()

  # manual filtering
  muts_tbl %>% filter(T_DEPTH >= 10 & N_DEPTH >= 10 & alt_counts >= 5 & T_FREQ > 0.03 & T_FREQ > (N_FREQ * 5))

}

# adjust indels to match ANNOVAR output
adjust_indels = function(x) {
  x %>%
    # set mutation type
    mutate(
      mut_type = case_when(
        str_length(REF) > str_length(ALT) ~ "del",
        str_length(ALT) > str_length(REF) ~ "ins",
        TRUE ~ "pt"
      )
    ) %>%
    # deletion (pos has to be incremented)
    mutate(
      REF = if_else(mut_type == "del", str_sub(REF, 2), REF),
      ALT = if_else(mut_type == "del", "-", ALT),
      POS = if_else(mut_type == "del", as.integer(POS + 1), POS)
    ) %>%
    # insertion
    mutate(
      REF = if_else(mut_type == "ins", "-", REF),
      ALT = if_else(mut_type == "ins", str_sub(ALT, 2), ALT)
    )
}

# split sample name for somatic variants
if (str_detect(sample_name, ":")) {
  sample_T = str_split_fixed(sample_name, ":", 2)[1]
  sample_N = str_split_fixed(sample_name, ":", 2)[2]
}

# import VCF as a vcfR object
muts_vcfr = read.vcfR(in_vcf, verbose = FALSE)

# check if there are any variants
if (nrow(muts_vcfr@fix) == 0) stop("no variants in imported VCF")

# determine variant caller based on VCF contents and parse accordingly
if (any(str_detect(muts_vcfr@meta, "##GATKCommandLine.HaplotypeCaller"))) {

  # GATK HaplotypeCaller
  message("parsing GATK HaplotypeCaller VCF")
  caller_type = "germline"
  vcf_tbl = parse_gatkhc(vcfr_obj = muts_vcfr, sample_name = sample_name)

} else if (any(str_detect(muts_vcfr@meta, "##source=lofreq"))) {

  # LoFreq
  message("parsing LoFreq VCF")
  caller_type = "germline"
  vcf_tbl = parse_lofreq(vcfr_obj = muts_vcfr, sample_name = sample_name)

} else if (any(str_detect(muts_vcfr@meta, "##source=Mutect2"))) {

  # Mutect 2.1 (GATK 4.0) or Mutect 2.2 (GATK 4.1)
  message("parsing Mutect 2.X VCF")
  caller_type = "somatic"
  vcf_tbl = parse_mutect2(vcfr_obj = muts_vcfr, sample_T = sample_T, sample_N = sample_N)

} else if (any(str_detect(muts_vcfr@meta, "##source=strelka"))) {

  # Strelka 2
  message("parsing Strelka 2 VCF")
  caller_type = "somatic"
  vcf_tbl = parse_strelka2(vcfr_obj = muts_vcfr, sample_T = sample_T, sample_N = sample_N)

} else {

  stop("unknown variant caller")

}

# check if table is empty
if (nrow(vcf_tbl) == 0) stop("output table is empty after parsing")

# update locale for sorting (to match the rest of the pipeline)
Sys.setlocale(category = "LC_ALL", locale = "C")

# create and sort by #MUT
vcf_tbl =
  vcf_tbl %>%
  mutate(`#MUT` = glue("{CHROM}:{POS}:{REF}:{ALT}")) %>%
  rename(CHR = CHROM) %>%
  arrange(`#MUT`)

# determine output table columns
if (caller_type == "germline") {
  vcf_tbl = vcf_tbl %>% mutate(SAMPLE = sample_name)
  out_cols = c("#MUT", "SAMPLE", "CHR", "POS", "QUAL", "DEPTH", "FREQ")
}
if (caller_type == "somatic") {
  vcf_tbl = vcf_tbl %>% mutate(SAMPLE_T = sample_T, SAMPLE_N = sample_N)
  out_cols = c("#MUT", "SAMPLE_T", "SAMPLE_N", "CHR", "POS", "QUAL", "T_DEPTH", "T_FREQ", "N_DEPTH", "N_FREQ")
}

# keep only relevant columns
clean_vcf_tbl = vcf_tbl[, out_cols]

# export
write_tsv(clean_vcf_tbl, path = out_txt)



# end
