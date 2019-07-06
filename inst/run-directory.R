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

parser <- ArgumentParser()

parser$add_argument("-i", "--in_dir", required=TRUE, help="directory containing input vcf files")
parser$add_argument("-or", "--out_rds", default="", help="rds containing fitted individual cells")
parser$add_argument("-hb", "--out_hap_bp", default="", help="break point bed for haploid cells")
parser$add_argument("-mb", "--out_m2d_bp", default="", help="break point bed for m2 haploid cells")
parser$add_argument("-hp", "--out_hap_plot", default="", help="plot for haploid cells")
parser$add_argument("-mp", "--out_m2d_plot", default="", help="plot for m2 haploid cells")
parser$add_argument("-c", "--out_cell_status", default="", help="cell status summary")
args <- parser$parse_args()

files <- list.files(args$in_dir, full.names=TRUE)
names(files) <- list.files(args$in_dir)

cells <- map(files, function(c) {
  cell <- sciLianti(c, bin_window_size=40, filter='rslow')
  cell <- trimSciLianti(cell)
  cell
})

if (args$out_rds != "") {
  saveRDS(cells, args$out_rds)
}
if (args$out_hap_bp != "") {
  hap_bp <- map_df(cells, ~.x$haploidBreakPoints, .id="CELL_ID")
  write_tsv(hap_bp, args$out_hap_bp)
}
if (args$out_hap_plot != "") {
  ps <- map(cells, plot_haploid_state)
  pdf(args$out_hap_plot, 10, 12)
  walk(ps, print)
  dev.off()
}
if (args$out_m2d_bp != "") {
  m2d_bp <- map_df(cells, ~.x$m2DiploidBreakPoints, .id="CELL_ID")
  write_tsv(m2d_bp, args$out_m2d_bp)
}
if (args$out_m2d_plot != "") {
  ps <- map(cells, plot_m2diploid_state)
  pdf(args$out_m2d_plot, 10, 12)
  walk(ps, print)
  dev.off()
}
if (args$out_cell_status != "") {
  status_df <- map_df(cells, function(x) {
    tibble(CELL_STATUS=paste(x$cellState, collapse="|"),
           HAPLOID_CHROMOSOMES=paste(x$haploidChromosomes, collapse="|"),
           M2DIPLOID_CHROMOSOMES=paste(x$m2DiploidChromosomes, collapse="|"))
  }, .id="CELL_ID")
  write_tsv(status_df, args$out_cell_status)
}