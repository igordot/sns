#!/usr/bin/env Rscript


##
## Extract and plot fragment size distribution from a BAM file.
##
## usage: Rscript --vanilla fragment-sizes.R sample_name file.bam
##


# increase output width
options(width = 150)

# relevent arguments
args = commandArgs(trailingOnly = TRUE)
sample_name = args[1]
bam_file = args[2]

# check that input files exist
if (!file.exists(bam_file)) stop("file does not exist: ", bam_file)

# load relevant packages
library(Rsamtools, quietly = TRUE)
library(readr, quietly = TRUE)
library(reshape2, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)

# import fragment size for only one mate of paired reads (other mate will be opposite)
scan_flag = scanBamFlag(isPaired = TRUE, isFirstMateRead = TRUE, isSecondMateRead = FALSE)
scan_params = ScanBamParam(what = c("isize"), flag = scan_flag)
scanned_sizes = scanBam(file = bam_file, param = scan_params)

# convert to vector
sizes = unlist(scanned_sizes[[1]]["isize"], use.names = FALSE)

# make all sizes positive (half should be negative)
sizes = abs(sizes)

# remove unreasonable values
sizes = sizes[!is.na(sizes)]
sizes = sizes[sizes > 0]
sizes = sizes[sizes < 5000]

# display stats
message("min: ", min(sizes))
message("max: ", max(sizes))
message("mean: ", round(mean(sizes), digits = 1))
message("median: ", median(sizes))

# stats table
stats_table = data.frame(SAMPLE = sample_name,
                         MEAN = round(mean(sizes), digits = 1),
                         MEDIAN = median(sizes),
                         SD = round(sd(sizes), digits = 1))
stats_table_file = paste0(sample_name, ".stats.csv")
write_csv(x = stats_table, path = stats_table_file)

# calculate frequency table
freq_table = as.data.frame(table(sizes), stringsAsFactors = FALSE)
colnames(freq_table) = c("#SIZE", sample_name)

# export frequency table
freq_table_file = paste0(sample_name, ".freq.csv")
write_csv(x = freq_table, path = freq_table_file)

# filter out fragments that will not be plotted (with bin width padding)
sizes = sizes[sizes < 1010]

# plot
plot_file = paste0(sample_name, ".png")
pdf(file = NULL)
ggplot(melt(sizes, value.name = "size"), aes(size)) +
geom_freqpoly(binwidth = 10, size = 1.5) +
scale_x_continuous(limits = c(min(sizes), 1000), breaks = 1:10*100) +
background_grid(major = "x", minor = "none") +
ggtitle(sample_name) +
ggsave(plot_file, width = 8, height = 5, units = "in")
dev.off()



# end
