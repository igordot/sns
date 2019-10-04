#!/usr/bin/env Rscript


##
## Extract and plot fragment size distribution from a BAM file.
##
## usage: Rscript --vanilla fragment-sizes.R sample_name file.bam
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
bam_file = args[2]

# check that input files exist
if (!file.exists(bam_file)) stop("file does not exist: ", bam_file)

# load relevant packages
load_install_packages("Rsamtools")
load_install_packages("tidyverse")
load_install_packages("glue")
load_install_packages("cowplot")

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
sizes = sizes[sizes < 2500]

# display stats
message("mean: ", round(mean(sizes), digits = 1))
message("median: ", median(sizes))
message("sd: ", round(sd(sizes), digits = 1))

# stats table
stats_table = data.frame(SAMPLE = sample_name,
                         MEAN = round(mean(sizes), digits = 1),
                         MEDIAN = median(sizes),
                         SD = round(sd(sizes), digits = 1))
stats_table_file = paste0(sample_name, ".stats.csv")
write_csv(x = stats_table, path = stats_table_file)

# create a table of fragment sizes
sizes_tbl = as_tibble(list(size = sizes)) %>%
  group_by(size) %>%
  summarize(n = n())

# change column names for export and export fragment sizes table
colnames(sizes_tbl) = c("#SIZE", sample_name)
write_csv(x = sizes_tbl, path = glue("{sample_name}.sizes.csv"))

# calculate size frequencies and remove rare occurrences (10/1M)
colnames(sizes_tbl) = c("size", "n")
sizes_tbl = sizes_tbl %>%
  mutate(freq = n/sum(n)) %>%
  filter(n > 10, freq > 0.00001, size < 1000)

# plot (all points with loess smoothing)
plot_file = paste0(sample_name, ".png")
freq_plot =
  ggplot(sizes_tbl, aes(x = size, y = freq)) +
  geom_point(size = 1.2, shape = 16, color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "black", size = 1.5) +
  scale_x_continuous(limits = c(0, 1020), expand = c(0, 0), breaks = 0:10 * 100) +
  scale_y_continuous(limits = c(0, max(sizes_tbl$freq) * 1.02), expand = c(0, 0)) +
  labs(title = sample_name, x = "Fragment Size", y = "Frequency") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  background_grid(major = "x", minor = "none")
save_plot(filename = plot_file, plot = freq_plot, base_height = 5, base_width = 8)

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")



# end
