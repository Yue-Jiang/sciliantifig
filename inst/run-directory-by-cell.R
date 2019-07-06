# TODO: make this into a library so don't have to source another file this way
#
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }
  
  if (!is.null(this.file)) return(dirname(this.file))
  
  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
  
  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))
  
  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}
current.dir = LocationOfThisScript()
##
source(file.path(current.dir, "sci-lianti.R"))

library(argparse)
library(tidyverse)
library(BiocParallel)
library(yaml)

parser <- ArgumentParser()

parser$add_argument("-i", "--in_dir", required=TRUE, help="directory containing input vcf files")
parser$add_argument("-or", "--out_rds_dir", default="", help="directory of rds containing fitted individual cells")
parser$add_argument("-hb", "--out_hap_bp_dir", default="", help="directory of break point bed for haploid cells")
parser$add_argument("-mb", "--out_m2d_bp_dir", default="", help="directory of break point bed for m2 haploid cells")
parser$add_argument("-hp", "--out_hap_plot_dir", default="", help="directory of plot for haploid cells")
parser$add_argument("-mp", "--out_m2d_plot_dir", default="", help="directory of plot for m2 haploid cells")
parser$add_argument("-c", "--out_cell_status_dir", default="", help="directory of cell status summary")
args <- parser$parse_args()

dir.create(args$out_rds_dir)
dir.create(args$out_hap_bp_dir)
dir.create(args$out_m2d_bp_dir)
dir.create(args$out_hap_plot_dir)
dir.create(args$out_m2d_plot_dir)
dir.create(args$out_cell_status_dir)

files <- list.files(args$in_dir, full.names=TRUE)
names(files) <- list.files(args$in_dir)

bplapply(names(files), function(file_name) {
  file <- files[[file_name]]
  base_name <- gsub("\\.vcf|\\.gz", "", file_name)
  cell <- sciLianti(file, bin_window_size=40, filter='rslow')
  if (args$out_rds_dir != "") {
    saveRDS(cell, file.path(args$out_rds_dir, paste0(base_name, ".rds")))
  }
  if (args$out_hap_bp_dir != "") {
    cell$haploidBreakPoints %>%
      add_column(CELL_ID=basename(cell$vcfName), .before=1) %>%
      write_tsv(file.path(args$out_hap_bp_dir, paste0(base_name, ".hb.tab")))
  }
  if (args$out_hap_plot_dir != "") {
    pdf(file.path(args$out_hap_plot_dir, paste0(base_name, ".hp.pdf")), 10, 12)
    print(plot_haploid_state(cell))
    dev.off()
  }
  if (args$out_m2d_bp_dir != "") {
    cell$m2DiploidBreakPoints %>%
      add_column(CELL_ID=basename(cell$vcfName), .before=1) %>%
      write_tsv(file.path(args$out_m2d_bp_dir, paste0(base_name, ".mb.tab")))
  }
  if (args$out_m2d_plot_dir != "") {
    pdf(file.path(args$out_m2d_plot_dir, paste0(base_name, ".mp.pdf")), 10, 12)
    print(plot_m2diploid_state(cell))
    dev.off()
  }
  if (args$out_cell_status_dir != "") {
    status_df <- tibble(
      CELL_ID=basename(cell$vcfName),
      CELL_STATUS=paste(cell$cellState, collapse="|"),
      HAPLOID_CHROMOSOMES=paste(cell$haploidChromosomes, collapse="|"),
      M2DIPLOID_CHROMOSOMES=paste(cell$m2DiploidChromosomes, collapse="|")
    )
    write_tsv(status_df, file.path(args$out_cell_status_dir, paste0(base_name, ".cellstatus.tab")))
    # yml <- yaml::as.yaml(cell[c("vcfName", "UID", "cellState", "haploidChromosomes", "m2DiploidChromosomes")])
    # writeLines(yml, con=file.path(args$out_cell_status_dir, paste0(base_name, ".cell-status.yaml")))
  }
})