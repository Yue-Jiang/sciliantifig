context("sciLianti pipeline, from vcf to called break points")

test_that("sciLianti pipeline can run from vcf to break points", {
  cell_file <- system.file("testdata", "yi190_ACGCGA.ATCGGGACCGGTCTT.filtered.vcf.gz", package="sciliantifig")
  test_cell <- sciLianti(cell_file, bin_window_size=40, filter='rs')
  trimmed_cell <- trimSciLianti(test_cell)
  restored_cell <- restoreSciLianti(trimmed_cell)
  expect_identical(test_cell, restored_cell)
})
