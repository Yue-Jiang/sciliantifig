suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
  library(BiocParallel)
})

parser <- ArgumentParser()
parser$add_argument("-i", "--input_rds", required=TRUE, help="input rds")
parser$add_argument("-c", "--cut_off", required=TRUE, type="double", help="threshold to call a chromosome UPD")
parser$add_argument("-o", "--output_upd", required=TRUE, help="per chromosome upd status table, CELL_ID CHROM UPD_STATUS")
args <- parser$parse_args()

cells <- readRDS(args$input_rds)
upd_df_list <- bplapply(cells, function(cell) {
  cell$alleleFrequency %>%
    group_by(CHROM) %>%
    summarise(FREQ_ALT=mean(FREQ_ALT, na.rm=TRUE)) %>%
    mutate(UPD_STATUS=case_when(FREQ_ALT > args$cut_off ~ "ALT",
                                FREQ_ALT < 1 - args$cut_off ~ "REF",
                                TRUE ~ "OTHER"))
})
upd_df <- bind_rows(upd_df_list, .id="CELL_ID")
write_tsv(upd_df, args$output_upd)

