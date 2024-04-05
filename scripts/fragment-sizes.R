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
message("num read pairs: ", length(sizes))

# make all sizes positive (half should be negative)
sizes = abs(sizes)

# remove unreasonable values
sizes = sizes[!is.na(sizes)]
sizes = sizes[sizes > 0]
sizes = sizes[sizes < 5000]
message("num fragments 0-5kb: ", length(sizes))

# display stats
message("mean size: ", round(mean(sizes), digits = 1))
message("median size: ", round(median(sizes)), digits = 1)
message("size sd: ", round(sd(sizes), digits = 1))

# save stats table
stats_table =
  tibble(
    SAMPLE = sample_name,
    MEAN = round(mean(sizes), digits = 1),
    MEDIAN = round(median(sizes), digits = 1),
    SD = round(sd(sizes), digits = 1)
  )
stats_table_file = paste0(sample_name, ".stats.csv")
write_csv(stats_table, file = stats_table_file)

# create a table of fragment sizes
sizes_tbl = tibble(size = sizes) %>% count(size, name = "n")

# change column names for export and export fragment sizes table
colnames(sizes_tbl) = c("#SIZE", sample_name)
write_csv(sizes_tbl, file = glue("{sample_name}.sizes.csv"))
Sys.sleep(1)

# calculate size frequencies and remove rare occurrences (1/1M)
colnames(sizes_tbl) = c("size", "n")
sizes_tbl = sizes_tbl %>%
  mutate(freq = n/sum(n)) %>%
  filter(n > 10, freq > (1 / 1000000))

# plot size frequency distribution (all points with loess smoothing)
plot_file = paste0(sample_name, ".sizes.1200.png")
freq_plot =
  filter(sizes_tbl, size <= 1200) %>%
  ggplot(aes(x = size, y = freq)) +
  geom_point(size = 1.2, shape = 16, color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "black", linewidth = 1.5) +
  scale_x_continuous(limits = c(0, 1200), expand = expansion(mult = c(0, 0.02)), breaks = 0:12 * 100) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = sample_name, x = "Fragment Size", y = "Frequency") +
  theme_cowplot() +
  theme(plot.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5)) +
  background_grid(major = "x", minor = "none")
ggsave(filename = plot_file, plot = freq_plot, width = 9, height = 5)
Sys.sleep(1)

# plot size frequency distribution (all points with loess smoothing)
plot_file = paste0(sample_name, ".sizes.600.png")
freq_plot =
  filter(sizes_tbl, size <= 600) %>%
  ggplot(aes(x = size, y = freq)) +
  geom_point(size = 1.2, shape = 16, color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "black", linewidth = 1.5) +
  scale_x_continuous(limits = c(0, 600), expand = expansion(mult = c(0, 0.02)), breaks = 0:6 * 100) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = sample_name, x = "Fragment Size", y = "Frequency") +
  theme_cowplot() +
  theme(plot.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5)) +
  background_grid(major = "x", minor = "none")
ggsave(filename = plot_file, plot = freq_plot, width = 8, height = 5)
Sys.sleep(1)

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")



# end
