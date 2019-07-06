#####################################
### functions for sci-lianti analysis
#####################################

#' The Cell object
#'
#' The Cell object that stores information about one cell.
#'
#' @param vcf_file The path to the input vcf file.
#'
#' @return A Cell object that is a list.
Cell <- function(vcf_file) {
  uid <- gsub("\\.vcf|\\.gz$", "", basename(vcf_file))
  object <- list(
    vcfName               = vcf_file,
    UID                   = uid,
    alleleFrequency       = tibble(),
    haploidState          = tibble(),
    haploidChromosomes    = character(0),
    haploidStateBlock     = tibble(),
    haploidBreakPoints    = tibble(),
    binnedAlleleFrequency = tibble(),
    binSize               = NA,
    m2DiploidState        = tibble(),
    m2DiploidStateBlock   = tibble(),
    m2DiploidChromosomes  = character(0),
    m2DiploidBreakPoints  = tibble(),
    cellState             = character(0)
  )
  class(object) <- "Cell"
  return(object)
}

#####################################
### use observed snps without binning
### works for haploid cells
#####################################

#' Calculate allele frequencies at each SNP for one cell
#'
#' Calculate allele frequencies at each SNP for one cell. Optionaly filter based on quality of the
#'   alternative allele, mapping quality etc.
#'
#' @import tidyverse
#' @keywords internal
#' @noRd
call_sc_allele_frequency <- function(cell, min_qual=0, min_amq=0, annotated_only=FALSE, named_only=FALSE) {
  vcf_df <- read_tsv(cell$vcfName, comment="##")
  colnames(vcf_df) <- gsub(".*/|#", "", colnames(vcf_df))
  chrs <- unique(vcf_df$CHROM)[!grepl("_|chrM", unique(vcf_df$CHROM))]
  chr_levels <- chrs[order(as.numeric(gsub("chr", "", chrs)))]
  sc_col <- colnames(vcf_df)[11]
  sc_df <- vcf_df %>%
    separate_(sc_col, c("GT", "ADF", "ADR", "LTDROP"), sep=":")
  if (annotated_only) sc_df <- filter(sc_df, grepl("rs", ID))
  if (named_only) sc_df <- filter(sc_df, grepl("^rs", ID))
  sc_df <- sc_df %>%
    mutate(AMQ=gsub("AMQ=", "", INFO)) %>%
    mutate(AMQ=gsub("\\.", "0", AMQ)) %>%
    separate(AMQ, c("AMQ_REF", "AMQ_ALT"), sep=",", remove=FALSE) %>%
    mutate(AMQ_REF=as.numeric(AMQ_REF), AMQ_ALT=as.numeric(AMQ_ALT)) %>%
    filter(AMQ_REF > min_amq | AMQ_ALT > min_amq) %>%
    filter(GT %in% c("./.", "0/0", "0/1", "1/0", "1/1")) %>%
    filter(QUAL >= min_qual) %>%
    separate(ADF, c("ADF1", "ADF2"), sep=",") %>%
    separate(ADR, c("ADR1", "ADR2"), sep=",") %>%
    mutate(ADF1=as.numeric(ADF1),
           ADF2=as.numeric(ADF2),
           ADR1=as.numeric(ADR1),
           ADR2=as.numeric(ADR2)) %>%
    mutate(AD=ADF1 + ADF2 + ADR1 + ADR2,
           ADFR1=ADF1 + ADR1,
           ADFR2=ADF2 + ADR2) %>%
    filter(AD > 0) %>%
    mutate(FREQ_ALT=case_when(.$ADFR1 > .$ADFR2 ~ 0,
                              .$ADFR1 < .$ADFR2 ~ 1,
                              TRUE ~ -1)) %>%
    filter(FREQ_ALT >= 0) %>% # due to loss currently impossible to cover the same snp on both chromosomes
    filter(!grepl("_", CHROM)) %>% # get rid of decoy
    filter(CHROM != "chrM") %>%
    mutate(CHROM=factor(CHROM, levels=chr_levels))
  cell$alleleFrequency <- sc_df
  return(cell)
}

#' Calculate allele frequencies at each SNP for one fly cell
#'
#' Calculate allele frequencies at each SNP for one fly cell. Optionaly filter based on quality of
#'  the alternative allele, mapping quality etc. Note that vcf format for fly is a little different
#'  than mouse, July 7 2018. Also per cell vcfs are not annotated, need to join with known vcf to
#'  filter for known snps.
#'
#' @import tidyverse
#' @keywords internal
#' @noRd
call_sc_allele_frequency_fly <- function(cell, min_qual=0, min_amq=0, annotated_only=FALSE, named_only=FALSE, known_vcf=NA) {
  full_vcf_df <- read_tsv(cell$vcfName, comment="##")
  colnames(full_vcf_df) <- gsub(".*/|#", "", colnames(full_vcf_df))
  # grab known snps only
  vcf_df <- known_vcf %>%
    mutate(CHROM=paste0("chr", .$CHROM)) %>%
    left_join(full_vcf_df, by=c("CHROM", "POS"), suffix=c("_HET", "")) %>%
    mutate(ALT_HET=toupper(ALT_HET), ALT=toupper(ALT)) %>%
    filter(ALT_HET == ALT)

  chrs <- unique(vcf_df$CHROM)[!grepl("_|chrM", unique(vcf_df$CHROM))]
  chr_levels <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
  sc_col <- colnames(vcf_df)[ncol(vcf_df)]

  sc_df <- vcf_df %>%
    separate_(sc_col, c("GT", "ADF", "ADR", "NOUSE"), sep=":") %>%
    mutate(AMQ=gsub("AMQ=", "", INFO)) %>%
    mutate(AMQ=gsub("\\.", "0", AMQ)) %>%
    separate(AMQ, c("AMQ_REF", "AMQ_ALT"), sep=",", remove=FALSE) %>%
    mutate(AMQ_REF=as.numeric(AMQ_REF), AMQ_ALT=as.numeric(AMQ_ALT)) %>%
    filter(AMQ_REF > min_amq | AMQ_ALT > min_amq) %>%
    filter(GT %in% c("./.", "0/0", "0/1", "1/0", "1/1")) %>%
    filter(QUAL >= min_qual) %>%
    separate(ADF, c("ADF1", "ADF2"), sep=",") %>%
    separate(ADR, c("ADR1", "ADR2"), sep=",") %>%
    mutate(ADF1=as.numeric(ADF1),
           ADF2=as.numeric(ADF2),
           ADR1=as.numeric(ADR1),
           ADR2=as.numeric(ADR2)) %>%
    mutate(AD=ADF1 + ADF2 + ADR1 + ADR2,
           ADFR1=ADF1 + ADR1,
           ADFR2=ADF2 + ADR2) %>%
    filter(AD > 0) %>%
    mutate(FREQ336=case_when(.$L336 == "1/1" & .$L382 == "0/0" ~ ADFR1 / AD,
                             .$L336 == "0/0" & .$L382 == "1/1" ~ ADFR2 / AD,
                             TRUE ~ -1)) %>%
    mutate(FREQ_ALT=case_when(.$FREQ336 > 0.7 ~ 0, # abuse of nomenclature to be consistent with downstream, ALT means 382
                              .$FREQ336 < 0.3 ~ 1,
                              TRUE ~ -1)) %>%
    filter(FREQ_ALT >= 0) %>%
    mutate(CHROM=factor(CHROM, levels=chr_levels))

  cell$alleleFrequency <- sc_df
  return(cell)
}

#' Run HMM to infer reference / alternative states
#'
#' @import tidyverse
#' @importFrom HMM initHMM viterbi
#' @keywords internal
#' @noRd
call_sc_state_hmm <- function(cell) {
  set.seed(12345)
  config_hmm <- function(transprob=1e-16) {
    states <- c("REF", "ALT", "HET")
    symbols <- c('0', '1')
    transprobs <- matrix(c(1 - transprob, transprob * 0.3, transprob * 0.7, transprob * 0.3, 1 - transprob, transprob * 0.7, transprob * 0.5, transprob * 0.5, 1 - transprob), nrow=3, byrow=TRUE)
    # ### mus sp het
    # mus
    # sp
    # het
    emissionprobs <- matrix(c(0.9, 0.1, 0.1, 0.9, 0.5, 0.5), nrow=3, byrow=TRUE)
    # ### 0 1
    # mus
    # sp
    # het
    list("States"=states,
         "Symbols"=symbols,
         "transProbs"=transprobs,
         "emissionProbs"=emissionprobs)
  }

  sc_df <- cell$alleleFrequency
  sc_state_df <- map_df(unique(sc_df$CHROM), function(chr) {
    chr_df <- sc_df %>%
      filter(CHROM == chr)
    hmm_config <- config_hmm(transprob=1e-10 / nrow(chr_df))
    this_hmm <- initHMM(States=hmm_config$States,
                        Symbols=hmm_config$Symbols,
                        transProbs=hmm_config$transProbs,
                        emissionProbs=hmm_config$emissionProbs)
    if (nrow(chr_df) < 10) {
      chr_df$STATE <- "UNDETERMINED"
    } else {
      chr_df$STATE <- HMM::viterbi(this_hmm, as.character(chr_df$FREQ_ALT))
    }
    chr_df %>%
      filter(STATE != "UNDETERMINED") %>%
      mutate(ALT_STATE=case_when(.$STATE == "REF" ~ 0,
                                 .$STATE == "ALT" ~ 1,
                                 .$STATE == "HET" ~ 0.5))
  })
  cell$haploidState <- sc_state_df
  return(cell)
}

#' Call haploid chromosomes
#'
#' @keywords internal
#' @noRd
call_haploid <- function(cell, min_hap_frac=0.8) {
  sc_state_df <- cell$haploidState
  chr_status_df <- sc_state_df %>%
    mutate(HAPLOID_FLAG=ALT_STATE != 0.5) %>%
    group_by(CHROM) %>%
    summarise(HAPLOID_FLAG=mean(HAPLOID_FLAG))
  cell$haploidChromosomes <- as.character(chr_status_df$CHROM)[chr_status_df$HAPLOID_FLAG >= min_hap_frac]
  return(cell)
}

#' Call blocks of chromosomes that belong to the same states
#'
#' @importFrom dplyr lag
#' @keywords internal
#' @noRd
call_state_block <- function(cell, min_block_size=1e5) {
  # chrom chromStart(block start) chromEnd(block end) name(cell id, add in map_df call) score(0,1,0.5) strand(+)
  hap_chrs_df <- filter(cell$haploidState, CHROM %in% cell$haploidChromosomes)
  if (nrow(hap_chrs_df) == 0) {
    cell$haploidStateBlock <- tibble()
    return(cell)
  }
  chrs <- unique(hap_chrs_df$CHROM)[!grepl("_", unique(hap_chrs_df$CHROM))]
  chr_levels <- chrs[order(as.numeric(gsub("chr", "", chrs)))]
  ret <- map_df(chrs, function(chr) {
    chr_df <- filter(hap_chrs_df, CHROM == chr)
    break_points <- which(chr_df$ALT_STATE != lag(chr_df$ALT_STATE))
    block_start_idxs <- c(1, break_points)
    block_end_idxs <- c(break_points - 1, nrow(chr_df))
    map2_df(block_start_idxs, block_end_idxs, function(s, e) {
      line <- tibble(chrom=unique(chr_df$CHROM),
                     chromStart=chr_df$POS[s],
                     chromEnd=chr_df$POS[e],
                     score=unique(chr_df$ALT_STATE[s:e]),
                     strand="+",
                     distance=chromEnd - chromStart) %>%
        mutate(flag=case_when(.$score == 0.5 ~ "DISCARD",
                              .$distance < min_block_size ~ "DISCARD",
                              TRUE ~ "KEEP"))
      line
    })
  })
  if (nrow(ret) > 0) {
    ret <- ret %>%
      mutate(chrom=factor(chrom, levels=chr_levels)) %>%
      arrange(chrom)
  }
  cell$haploidStateBlock <- ret
  return(cell)
}

#' Call break points on chromosomes
#'
#' Call break points on chromosomes, where the chromosome states change from one to another.
#'
#' @import tidyverse
#' @keywords internal
#' @noRd
call_break_points <- function(cell) {
  use_block_df <- filter(cell$haploidStateBlock, flag == "KEEP")
  if (nrow(use_block_df) == 0) {
    cell$haploidBreakPoints <- tibble()
    return(cell)
  }
  chrs <- unique(use_block_df$chrom)
  chr_levels <- chrs[order(as.numeric(gsub("chr", "", chrs)))]
  ret <- map_df(chrs, function(chr) {
    chr_df <- filter(use_block_df, chrom == chr) %>%
      arrange(chromStart)
    if (nrow(chr_df) < 2) return(NULL)
    tibble(chrom=chr,
           chromStart=chr_df$chromEnd[1:(nrow(chr_df) - 1)],
           chromEnd=chr_df$chromStart[2:(nrow(chr_df))],
           score=chr_df$score[2:(nrow(chr_df))] - chr_df$score[1:(nrow(chr_df) - 1)],
           strand='+',
           distance=chromEnd - chromStart)
  })
  if (nrow(ret) > 0) {
    ret <- ret %>%
      mutate(chrom=factor(chrom, levels=chr_levels)) %>%
      arrange(chrom)
  }
  cell$haploidBreakPoints <- ret
  return(cell)
}

#' Plot haploid state for a cell
#'
#' @import tidyverse
#' @keywords internal
#' @noRd
plot_haploid_state <- function(cell, ncol=3) {
  if (!any(class(cell$haploidState) == 'data.frame')) cell <- restoreSciLianti(cell)
  df <- cell$haploidState
  if (length(cell$haploidChromosomes) > 0) {
    df <- df %>%
      mutate(CHROM_STATUS=case_when(.$CHROM %in% cell$haploidChromosomes ~ "HAPLOID",
                                    TRUE ~ "NOT HAPLOID"))
  }
  p <- ggplot(df, aes(POS, FREQ_ALT)) +
    geom_point(alpha=0.1, size=0.4) +
    geom_line(aes(POS, ALT_STATE - (ALT_STATE - 0.5) / 4), color="red") +
    theme_bw() +
    ylim(c(0, 1))
  if (length(cell$haploidChromosomes) > 0) {
    p <- p +
      facet_wrap(~CHROM + CHROM_STATUS, scales="free_x", ncol=ncol)
  } else {
    p <- p +
      facet_wrap(~CHROM, scales="free_x", ncol=ncol)
  }
  return(p)
}

#####################################
### bin snps
### works for M2 diploid cells
#####################################

#' Calculate allele frequencies for binned SNPs for selected window size for one fly cell
#'
#' @import tidyverse
#' @keywords internal
#' @noRd
sc_allele_freq_binned <- function(cell, window_size) {
  binned_allele_freq <- cell$alleleFrequency %>%
    group_by(CHROM) %>%
    arrange(POS) %>%
    mutate(window=round((row_number() - 1)/ window_size)) %>%
    ungroup()%>%
    group_by(CHROM, window) %>%
    summarise(POS=median(POS), FREQ_ALT=mean(FREQ_ALT)) %>%
    ungroup()
  cell$binnedAlleleFrequency <- binned_allele_freq
  cell$binSize <- window_size
  return(cell)
}

#' Run HMM to infer reference / alternative states from binned allele frequencies
#'
#' @import tidyverse
#' @importFrom depmixS4 depmix setpars fit posterior
#' @importFrom stats gaussian
#' @keywords internal
#' @noRd
call_sc_state_hmm_binned <- function(cell) {
  set.seed(12345)
  sc_binned_df <- cell$binnedAlleleFrequency
  sc_state_df <- map_df(unique(sc_binned_df$CHROM), function(chr) {
    chr_df <- sc_binned_df %>%
      filter(CHROM == chr)
    hmm_init <- depmix(list(FREQ_ALT~1), data=chr_df, nstates=3, family=list(gaussian()))
    transprob <- 1e-16 #/ nrow(chr_df)
    pars <- c(1/3, 1/3, 1/3,
              1 - transprob, transprob, 0, transprob / 2, 1 - transprob, transprob / 2, 0, transprob, 1 - transprob,
              0.05, 0.1, 0.5, 0.1, 0.95, 0.1)
    hmm_init <- setpars(hmm_init, pars)
    fixed <- c(0, 0, 0,
               rep(1, 9),
               1, 0, 1, 0, 1, 0)
    # # 0 0.5 1
    # 0
    # 0.5
    # 1
    chr_df$STATE <- "UNDETERMINED"
    try({
      # this is pretty aggressive error handling, since depmixS4 sometimes just fail to fit
      suppressMessages(hmm_fit<-fit(hmm_init, fixed=fixed, verbose=FALSE))
      hmm_post<-depmixS4::posterior(hmm_fit)
      chr_df$STATE <- c("REF", "HET", "ALT")[hmm_post$state]
    })
    chr_df %>%
      mutate(ALT_STATE=case_when(.$STATE == "REF" ~ 0,
                                 .$STATE == "ALT" ~ 1,
                                 .$STATE == "HET" ~ 0.5))
  })
  cell$m2DiploidState <- sc_state_df
  return(cell)
}

#' Call blocks of chromosomes that belong to the same states for binned allele frequencies
#'
#' @importFrom dplyr lag
#' @keywords internal
#' @noRd
call_state_block_binned <- function(cell, min_block_size=1e5, add=FALSE) {
  # Note this is done on all chromosomes
  sc_status_df <- cell$m2DiploidState
  chrs <- unique(sc_status_df$CHROM)[!grepl("_", unique(sc_status_df$CHROM))]
  chr_levels <- chrs[order(as.numeric(gsub("chr", "", chrs)))]
  ret <- map_df(chrs, function(chr) {
    chr_df <- filter(sc_status_df, CHROM == chr)
    break_points <- which(chr_df$ALT_STATE != lag(chr_df$ALT_STATE))
    block_start_idxs <- c(1, break_points)
    block_end_idxs <- c(break_points - 1, nrow(chr_df))
    block_df <- map2_df(block_start_idxs, block_end_idxs, function(s, e) {
      line <- tibble(chrom=unique(chr_df$CHROM),
                     chromStart=chr_df$POS[s],
                     chromEnd=chr_df$POS[e],
                     score=unique(chr_df$ALT_STATE[s:e]),
                     strand="+",
                     distance=chromEnd - chromStart) %>%
        mutate(flag=case_when(.$distance < min_block_size ~ "DISCARD",
                              TRUE ~ "KEEP"))
      line
    })
    if (add) {
      chr_df$BLOCK <- 0
      chr_df$BLOCK_FLAG <- "DISCARD"
      for (i in seq_len(nrow(block_df))) {
        if (block_df$flag[i] == "KEEP") {
          chr_df[which(chr_df$POS >= block_df$chromStart[i] & chr_df$POS <= block_df$chromEnd[i]), "BLOCK"] <- i
          chr_df[which(chr_df$POS >= block_df$chromStart[i] & chr_df$POS <= block_df$chromEnd[i]), "BLOCK_FLAG"] <- "KEEP"
        }
      }
      chr_df
    } else {
      block_df
    }
  })
  if (nrow(ret) > 0) {
    ret <- ret %>%
      mutate(chrom=factor(chrom, levels=chr_levels)) %>%
      arrange(chrom)
  }
  cell$m2DiploidStateBlock <- ret
  return(cell)
}

#' Call m2 diploid chromosomes
#'
#' @import tidyverse
#' @keywords internal
#' @noRd
call_m2diploid <- function(cell, min_het_length=1e7, min_homo_length=1e7) {
  # Note this is called after removing spiky changes
  # because it requires absolute block length other than fraction
  sc_stateblock_df <- cell$m2DiploidStateBlock
  df <- sc_stateblock_df %>%
    filter(flag == "KEEP") %>%
    drop_na(score) %>%
    arrange(chrom, chromStart)
  prev_row <- filter(df, chrom == "to get an empty row")
  df_list <- list()
  for (i in seq_len(nrow(df))) {
    if (nrow(prev_row) == 0) {
      prev_row <- df[i, ]
      next
    }
    if (prev_row$chrom == df$chrom[i] && prev_row$score == df$score[i]) {
      prev_row[1, "chromEnd"] <- df$chromEnd[i]
      prev_row[1, "distance"] <- prev_row$chromEnd[1] - prev_row$chromStart[1]
      next
    }
    if (prev_row$chrom != df$chrom[i] || prev_row$score != df$score[i]) {
      df_list[[length(df_list) + 1]] <- prev_row
      prev_row <- df[i, ]
      next
    }
  }
  df_list[[length(df_list) + 1]] <- prev_row
  cleaned_df <- bind_rows(df_list)
  homo_df <- cleaned_df %>%
    filter(score %in% c(0, 1)) %>%
    group_by(chrom) %>%
    summarise(HOMO_FLAG=max(distance >= min_homo_length)) %>%
    filter(HOMO_FLAG == TRUE)
  het_df <- cleaned_df %>%
    filter(score == 0.5) %>%
    group_by(chrom) %>%
    summarise(HET_FLAG=max(distance >= min_het_length)) %>%
    filter(HET_FLAG == TRUE)
  cell$m2DiploidChromosomes <- intersect(homo_df$chrom, het_df$chrom)
  return(cell)
}

#' A likelihood-ratio test used to call chromosome break points for m2 diploid chromosomes
#'
#' @import tidyverse
#' @keywords internal
#' @noRd
breakpoint_by_lr <- function(obs, states, error_prob=1e-3, lr_threshold=2) {
  # obs can be 0 or 1
  # states (vector of 2, can be 0, 0.5, 1)
  # error_prob: prob snp is called wrong
  prob <- matrix(c(1 - error_prob, error_prob, 0.5, 0.5, error_prob, 1 - error_prob), nrow=3, byrow=TRUE)
  rownames(prob) <- as.character(c(0, 0.5, 1))
  colnames(prob) <- as.character(c(0, 1))
  len <- length(obs)
  state_left <- states[1]
  state_right <- states[2]
  lls <- c()
  for (i in seq_len(len - 1)) {
    # break after i
    obs_left <- obs[1:i]
    obs_right <- obs[(i + 1):len]
    ll <- log10(prob[as.character(state_left), "0"]) * sum(obs_left == 0) +
      log10(prob[as.character(state_left), "1"]) * sum(obs_left == 1) +
      log10(prob[as.character(state_right), "0"]) * sum(obs_right == 0) +
      log10(prob[as.character(state_right), "1"]) * sum(obs_right == 1)
    lls <- c(lls, ll)
  }
  return(which(lls >= max(lls) - lr_threshold))
}

#' Call chromosome break points for m2 diploid chromosomes
#'
#' @import tidyverse
#' @keywords internal
#' @noRd
call_break_points_binned <- function(cell, check_range=20) {
  sc_df <- cell$alleleFrequency
  state_block_df <- cell$m2DiploidStateBlock
  # given the break points found by binning snps, go back to per snp table and refine where the break point actually is
  use_block_df <- filter(state_block_df, flag == "KEEP")
  if (nrow(use_block_df) == 0) {
    cell$m2DiploidBreakPoints <- tibble()
    return(cell)
  }
  chrs <- unique(use_block_df$chrom)
  chr_levels <- chrs[order(as.numeric(gsub("chr", "", chrs)))]

  ret <- map_df(chrs, function(chr) {
    chr <- as.character(chr)
    chr_block_df <- filter(use_block_df, chrom == chr) %>%
      arrange(chromStart)
    chr_sc_df <- filter(sc_df, CHROM == chr) %>%
      arrange(POS)
    if (nrow(chr_block_df) < 2) return(NULL)
    map_df(2:nrow(chr_block_df), function(i) {
      mid_pos <- chr_block_df$chromStart[i]
      check_left <- max(1, which(chr_sc_df$POS == mid_pos) - check_range)
      check_right <- min(nrow(chr_sc_df), which(chr_sc_df$POS == mid_pos) + check_range)
      break_idx <- breakpoint_by_lr(obs=chr_sc_df$FREQ_ALT[check_left:check_right],
                                    states=c(chr_block_df$score[i - 1], chr_block_df$score[i]),
                                    error_prob=1e-3,
                                    lr_threshold=3)
      left <- chr_sc_df$POS[c(check_left:check_right)[min(break_idx)]]
      right <- chr_sc_df$POS[c(check_left:check_right)[max(break_idx)]]
      tibble(chrom=as.character(chr),
             chromStart=left,
             chromEnd=right,
             score=chr_block_df$score[i] - chr_block_df$score[i - 1],
             strand='+',
             distance=chromEnd - chromStart)
    })
  })
  if (nrow(ret) > 0) {
    ret <- ret %>%
      mutate(chrom=factor(chrom, levels=chr_levels)) %>%
      arrange(chrom)
  }
  cell$m2DiploidBreakPoints <- ret
  return(cell)
}

#' Plot m2-diploid state for a cell
#'
#' @import tidyverse
#' @keywords internal
#' @noRd
plot_m2diploid_state <- function(cell, ncol=3) {
  if (!any(class(cell$m2DiploidState) == 'data.frame')) cell <- restoreSciLianti(cell)
  df <- cell$m2DiploidState
  if (length(cell$m2DiploidChromosomes) > 0) {
    df <- df %>%
      mutate(CHROM_STATUS=case_when(.$CHROM %in% cell$m2DiploidChromosomes ~ "M2DIPLOID",
                                    TRUE ~ "NOT M2DIPLOID"))
  }
  p <- ggplot(df, aes(POS, FREQ_ALT)) +
    geom_point(alpha=0.2) +
    geom_line(aes(POS, ALT_STATE), color="red") +
    theme_bw()
  if (length(cell$m2DiploidChromosomes) > 0) {
    p <- p +
      facet_wrap(~CHROM + CHROM_STATUS, scales="free_x", ncol=ncol)
  } else {
    p <- p +
      facet_wrap(~CHROM, scales="free_x", ncol=ncol)
  }
  return(p)
}

#' Call the state for a cell
#'
#' @import tidyverse
#' @keywords internal
#' @noRd
call_cell_state <- function(cell) {
  cell$cellState <- character(0)
  if (length(cell$haploidChromosomes) >= 15) cell$cellState <- c(cell$cellState, "haploid")
  if (length(cell$m2DiploidChromosomes) >= 4) cell$cellState <- c(cell$cellState, "m2diploid")
  return(cell)
}

#####################################
### wrapper, runs one cell
#####################################

#' Analyze break points for one cell
#'
#' Analyze break points for one cell.
#'
#' @param single_cell_file Path to the input vcf file.
#' @param bin_window_size Window size for calculating binned allele frequencies.
#' @param filter Filter SNPs based on annotation. One of c('none', 'rs', 'rslow').
#' @param organism One of c('mouse', 'fly').
#' @param known_vcf Path to a vcf file containing known SNPs. Only used when organism is 'fly' and when we
#'   want to filter by known SNPs.
#'
#' @return A Cell object.
#'
#' @examples
#' \dontrun{
#' cell <- sciLianti(single_cell_file="/path/to/input/vcf", bin_window_size=40, filter='rslow',
#'  organism='mouse', known_vcf=NA)
#' }
#'
#' @export
sciLianti <- function(single_cell_file, bin_window_size, filter='none', organism='mouse', known_vcf=NA) {
  cell <- Cell(single_cell_file)
  # filter
  # if rs, filter named snps only e.g. rs012345
  # if rslow, filter named snps and .rsLOW
  # if anything else, don't filter for known snps
  if (organism == 'fly') {
    cell <- call_sc_allele_frequency_fly(cell, known_vcf=known_vcf)
  } else {
    if (tolower(filter) == 'rs') {
      cell <- call_sc_allele_frequency(cell, named_only=TRUE)
    } else if (tolower(filter) == 'rslow') {
      cell <- call_sc_allele_frequency(cell, annotated_only=TRUE)
    } else {
      cell <- call_sc_allele_frequency(cell, annotated_only=FALSE, named_only=FALSE)
    }
  }
  cell <- call_sc_state_hmm(cell)
  cell <- call_haploid(cell)
  cell <- call_state_block(cell)
  cell <- call_break_points(cell)
  cell <- sc_allele_freq_binned(cell, bin_window_size)
  cell <- call_sc_state_hmm_binned(cell)
  cell <- call_state_block_binned(cell)
  cell <- call_m2diploid(cell)
  cell <- call_break_points_binned(cell)
  cell <- call_cell_state(cell)
  return(cell)
}

#' Make a Cell object smaller in size by removing data that can be relatively easily regenerated,
#' for storage purpose
#'
#' @param cell A processed full-size Cell object.
#'
#' @return A trimmed Cell object.
#'
#' @examples
#' \dontrun{
#' cell <- trimSciLianti(cell)
#' }
#'
#' @export
trimSciLianti <- function(cell) {
  # remove data that can be relatively easily regenerated, for storage purpose
  cell$haploidState <- NA
  cell$binnedAlleleFrequency <- NA
  cell$m2DiploidState <- NA
  return(cell)
}

#' Restore a trimmed Cell object to full size
#'
#' @param cell A trimmed Cell object.
#' @param bin_window_size The window size to calculate binned allele frequencies. Will use stored
#'   old value if left as NA.
#'
#' @return A full-size Cell object.
#'
#' @examples
#' \dontrun{
#' cell <- restoreSciLianti(cell, bin_window_size=NA)
#' }
#'
#' @export
restoreSciLianti <- function(cell, bin_window_size=NA) {
  if (!is.na(bin_window_size)) {
    cell$binSize <- bin_window_size
  }
  must_have <- c("vcfName", "UID", "alleleFrequency", "binSize")
  missing <- map_lgl(cell[must_have], is.null)
  if (any(missing)) {
    stop(sprintf("%s missing, cannot restore cell.", paste(must_have[missing], collapse=", ")))
  }
  cell <- call_sc_state_hmm(cell)
  cell <- call_haploid(cell)
  cell <- call_state_block(cell)
  cell <- call_break_points(cell)
  cell <- sc_allele_freq_binned(cell, cell$binSize)
  cell <- call_sc_state_hmm_binned(cell)
  cell <- call_state_block_binned(cell)
  cell <- call_m2diploid(cell)
  cell <- call_break_points_binned(cell)
  cell <- call_cell_state(cell)
  return(cell)
}
