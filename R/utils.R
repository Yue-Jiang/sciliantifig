#' Locate a file installed with this library
#'
#' @param f The name of a file installed with this library.
#'
#' @return The reconstructed full file path of said file.
#'
#' @export
slfile <- function(f) {
  ret <- system.file("extdata", f, package = "sciliantifig")
  if (!file.exists(ret)) {
    warning("Cannot locate %s in sciliantifig. It's not distributed with the library in `inst/extdata/`.", f)
  }
  return(ret)
}

#' Calculate the distance between to break points
#'
#' @param s1 Start point of break point 1
#' @param e1 End point of break point 1
#' @param s2 Start point of break point 2
#' @param e2 End point of break point 2
#'
#' @return A numeric value for the distance between the two break points.
#'
#' @export
bp_distance <- function(s1, e1, s2, e2) {
  # if (!((s1 < s2 & e1 < s2) | (s1 > e2 & s1 > s2))) warning(str(c(s1, e1, s2, e2)))
  assertthat::assert_that((s1 < s2 & e1 < s2) | (s1 > e2 & s1 > s2))
  if (e1 < s2) return(s2 - e1)
  if (s1 > e2) return(s1 - e2)
}

#' Check that the two break points are sane, in that one does not fall exclusively within the other
#'
#' @param s1 Start point of break point 1
#' @param e1 End point of break point 1
#' @param s2 Start point of break point 2
#' @param e2 End point of break point 2
#'
#' @return Logical. TRUE if the the two break points are sane, FALSE if not.
#'
#' @export
bp_distance_check <- function(s1, e1, s2, e2) {
  # return false if problematic
  (s1 < s2 & e1 < s2) | (s1 > e2 & s1 > s2)
}

#' Beautified version of summary plot for BAS
#'
#' @param x A BAS model object.
#' @param top.models The number of top models to include in the plot.
#'
#' @return A gtable.
#'
#' @import tidyverse gridExtra
#' @export
bas_plot <- function (x, top.models = 20) {
  # require(tidyverse)
  xcoef <- coef(x)
  incprob <- tibble(Predictor=xcoef$namesx, InclusionProb=xcoef$probne0)

  postprob <- x$postprobs
  top.models <- min(top.models, x$n.models)
  best <- order(-x$postprobs)[1:top.models]
  postprob <- postprob[best]/sum(postprob[best])
  which.mat <- list2matrix.which(x, best)
  nvar <- ncol(which.mat)
  subset <- 1:nvar
  which.mat <- which.mat[, subset, drop = FALSE]
  df <- as.data.frame(which.mat) %>%
    mutate(Rank=row_number()) %>%
    gather(Predictor, Included, 1:nvar) %>%
    left_join(tibble(Rank=1:top.models, PosteriorProb=postprob)) %>%
    mutate(PosteriorProb=case_when(Included == 0 ~ -1,
                                   TRUE ~ PosteriorProb)) %>%
    left_join(incprob) %>%
    group_by(Predictor) %>%
    mutate(PredIncludedMean=mean(Included == 1)) %>%
    ungroup() %>%
    arrange(desc(InclusionProb)) %>%
    mutate(Predictor=factor(Predictor, levels=rev(unique(.$Predictor))))
  p_left <- ggplot(df, aes(Rank, Predictor, fill=PosteriorProb)) +
    geom_tile() +
    scale_fill_gradient2(low="grey50", mid="white", high="red", limits=c(0, 1)) +
    xlim(1, top.models) +
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.justification = c(1, 0),
          legend.position = c(0.99, 0.01),
          legend.box.background = element_rect(colour="black"))
  p_right <- distinct(df, Predictor, InclusionProb) %>%
    ggplot(aes(Predictor, InclusionProb)) +
    geom_bar(stat='identity', fill='dodgerblue') +
    coord_flip() +
    theme_bw() +
    xlab("") +
    ylab("Inclusion probability") +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.text.y=element_blank(),
          axis.line=element_line(size=0.2)
    )

  p <- gridExtra::grid.arrange(p_left, p_right, nrow=1, widths=c(80, 25))
  return(p)
}

#' Color cell attributes on scatter plot of the first two PCs
#'
#' @param cell_pca The PCA object.
#' @param cell_df The data frame of cell attributes.
#'
#' @return a list of ggplots.
#'
#' @export
plot_cell_pca <- function(cell_pca, cell_df) {
  # require(viridis)
  # require(gridExtra)
  s <- summary(cell_pca)
  s1 <- s$importance["Proportion of Variance", "PC1"] * 100
  s2 <- s$importance["Proportion of Variance", "PC2"] * 100
  cell_pca_df <- cell_df %>%
    mutate(PC1=cell_pca$x[, "PC1"],
           PC2=cell_pca$x[, "PC2"])
  cols <- setdiff(names(cell_pca_df)[-1], c("PC1", "PC2"))
  plts <- map(cols, function(x) {
    if (is.numeric(cell_pca_df[[x]])) {
      this_df <- cell_pca_df
      this_df[[x]] <- log1p(cell_pca_df[[x]])
    } else {
      this_df <- cell_pca_df
    }
    ggplot(this_df, aes(PC1, PC2)) +
      geom_point(aes_string(color=x), alpha=0.6) +
      scale_color_viridis(discrete=!is.numeric(this_df[[x]])) +
      xlab(sprintf("PC1 (%1.2f%%)", s1)) +
      ylab(sprintf("PC2 (%1.2f%%)", s2)) +
      theme_classic() +
      theme(legend.position="bottom")
  })
  return(plts)
}
