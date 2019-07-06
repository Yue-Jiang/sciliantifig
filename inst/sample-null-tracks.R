suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
  library(BiocParallel)
})

parser <- ArgumentParser()
parser$add_argument("-i", "--observed_tracks", required=TRUE, help="observed tracks: haploid.no0.hb.filtered.tab")
parser$add_argument("-r", "--reference_all_snps", required=TRUE, help="reference bed containing all possible track ends: spret.snp.pos.includeLOW.gz")
parser$add_argument("-p", "--n_proc", required=FALSE, default=0, type="integer", help="number of parallel processes to use")
parser$add_argument("-n", "--n_tracks_to_sample", required=FALSE, default=-1, type="integer", help="number of tracks to simulate, will be set equal to observed tracks if not provided")
parser$add_argument("-m", "--repeat_times", required=FALSE, default=1, type="integer", help="repeat the sampling m times, resulting in m files")
parser$add_argument("-o", "--output_dir", required=TRUE, help="output directory for sampled null tracks")
args <- parser$parse_args()

if (args$n_proc > 0) register(MulticoreParam(args$n_proc))
if (!dir.exists(args$output_dir)) dir.create(args$output_dir)
n_id <- floor(log10(args$repeat_times)) + 1

# obs_df <- read_tsv("~/Desktop/sci-lianti-files/haploid.no0.hb.filtered.tab")
# ref <- read_tsv("~/Desktop/sci-lianti-files/spret.snp.pos.includeLOW.gz", comment="##")

obs_df <- read_tsv(args$observed_tracks)
ref <- read_tsv(args$reference_all_snps, comment="##")
colnames(ref) <- gsub("#", "", colnames(ref))

chr_size <- ref %>%
  group_by(CHROM) %>%
  summarise(MIN_POS=min(POS), MAX_POS=max(POS)) %>%
  ungroup() %>%
  mutate(CHROM_SIZE= MAX_POS - MIN_POS) %>%
  filter(!CHROM %in% c("chrX", "chrY"))

chr_prob <- chr_size$CHROM_SIZE / sum(chr_size$CHROM_SIZE)
names(chr_prob) <- chr_size$CHROM

ref_per_chr <- map(set_names(chr_size$CHROM), function(chr) {
  filter(ref, CHROM == chr)$POS
})

remove(ref)

n_samples <- ifelse(args$n_tracks_to_sample > 0, args$n_tracks_to_sample, nrow(obs_df))

chr_levels <- chr_size$CHROM[order(as.numeric(gsub("chr", "", chr_size$CHROM)))]

draw_chr <- function() {
  sample(x=names(chr_prob), size=1, prob=chr_prob)
}

draw_segment <- function(chr) {
  this_chr <- filter(chr_size, CHROM == chr)
  s <- sample(this_chr$MIN_POS:this_chr$MAX_POS, 1)
  len <- sample(obs_df$DISTANCE, 1)
  e <- s + len
  return(c("s"=s, "e"=e))
}

coerce_to_closest_snp <- function(chr, pos) {
  all_pos <- ref_per_chr[[chr]]
  new_pos <- all_pos[which.min(abs(all_pos - pos))]
  return(new_pos)
}

for(i in 1:args$repeat_times) {
  samples <- bplapply(1:n_samples, function(x) {
    n_retry <- 0
    while (n_retry < 1000) {
      chr <- draw_chr()
      seg <- draw_segment(chr)
      ret <- tibble(
        CHROM=chr,
        START=coerce_to_closest_snp(chr, seg["s"]),
        END=coerce_to_closest_snp(chr, seg["e"])
      )
      if (ret$START != ret$END) return(ret)
      n_retry <- n_retry + 1
    }
    if (n_retry == 1000) stop("Having trouble getting tracks with length > 0")
  })
  out_df <- bind_rows(samples) %>%
    mutate(CHROM=factor(CHROM, levels=chr_levels)) %>%
    arrange(CHROM, START, END)
  out_file <- file.path(args$output_dir, paste0("sampled-track_", str_pad(i, n_id, pad="0"), ".tab"))
  print(out_file)
  write_tsv(out_df, out_file)
}
