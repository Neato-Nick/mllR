
#' Filter a vcfR object
#'
#' @param x a vcfR object
#' @param min_depth Minimum depth a genotype call must be to not be censored.
#' @param max_depth Maximum depth a genotype call can be or else it wil be censored.
#' @param min_mq Minimum mapping quality a genotype call can be or else it will be censored.
#' @param samples_miss1 Maximum percent missing data a sample can have before it gets omitted. Occurs before variant omission.
#' @param samples_miss2 Maximum percent missing data a sample can have before it gets omitted. Occurs after variant omission.
#' @param vars_miss Maximum percent missing data a variant can have before it gets omitted. Occurs after first wave of sample omission.
#' @param top_95_quantile_filter Censor the 5th and 95th percentiles of genotype call depth or not.
#' @param remove_nonpolymorphs After removing samples, some variants will no longer be polymorphic. These should be removed to not inflate the true number of variants.
#' @return A filtered \code{vcfR} object.
#' @keywords vcf filter
#' @export
#' @examples
#' filtered_vcf <- filter_vcf(unfiltered_vcf)
#' filtered_vcf <- filter_vcf(unfiltered_vcf, samples_miss1 = 0.8, samples_miss2 = 0.4, vars_miss = 0.001)
#'
filter_vcf <- function (x, min_depth = 4, max_depth = 800, min_mq = 20, samples_miss1 = 0.9999, samples_miss2 = 0.7, vars_miss = 0.2, top_dp_quantile_filter = TRUE, remove_nonpolymorphs = TRUE) {
  dp <- vcfR::extract.gt(x, element = "DP", as.numeric = TRUE)
  # Censor variants depth < 8 or > 100
  dp[dp < min_depth] <- NA
  dp[dp > max_depth] <- NA
  x.2 <- x
  x.2@gt[,-1][ is.na(dp) == TRUE] <- NA

  # Censor variants in top 95% quantile
  if (top_dp_quantile_filter == TRUE) {
    dp <- vcfR::extract.gt(x.2, element = "DP", as.numeric = TRUE)
    quants <- apply(dp, MARGIN=2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)
    # Sweep up the file to get rid of the quartiles
    dp2 <- sweep(dp, MARGIN=2, FUN="-", quants[1,])
    dp[dp2 < 0] <- NA
    dp2 <- sweep(dp, MARGIN=2, FUN="-", quants[2,])
    dp[dp2 > 0] <- NA
    x.3 <- x.2
    x.3@gt[,-1][ is.na(dp) == TRUE] <- NA
  } else {
    x.3 <- x.2
  }

  # Censor variants with mappint quality < min_mq
  mq <- extract.info(x, element = "MQ", as.numeric = TRUE)
  # bypass mapping quality filter if there are missing values in MQ
  # any NA values here would indicate denovo assembly
  if (is.na(sum(mq)) == TRUE) {
    x.4 <- x.3
  } else {
    mq[mq < min_mq] <- NA
    x.4 <- x.3
    x.4@gt[,-1][ is.na(mq) == TRUE] <- NA
  }

  # Omit samples with > samples_miss1
  dp <- extract.gt(x.4, element = "DP", as.numeric = TRUE)
  # Now omitting samples
  myMiss <- apply(dp, MARGIN = 2, function(y){ sum( is.na(y) ) } )
  myMiss <- myMiss / nrow(dp)
  x.5 <- x.4
  x.5@gt <- x.5@gt[, c(TRUE, myMiss < samples_miss1)]

  # Omit variants with > vars_miss
  dp <- extract.gt(x.5, element = "DP", as.numeric = TRUE)
  # Now omitting vars
  myMiss <- apply(dp, MARGIN = 1, function(y) {sum(is.na(y))})
  myMiss <- myMiss / ncol(dp)
  x.6 <- x.5
  x.6 <- x.6[myMiss < vars_miss, ]

  # Omit samples with > samples_miss2
  dp <- extract.gt(x.6, element = "DP", as.numeric = TRUE)
  # Omit samples
  myMiss <- apply(dp, MARGIN = 2, function(y) {sum(is.na(y))})
  myMiss <- myMiss / ncol(dp)
  x.7 <- x.6
  x.7@gt <- x.7@gt[, c(TRUE, myMiss < samples_miss2)]

  # Now remove all variants that are no longer present
  # Measured by a minor allele count > 0
  if (remove_nonpolymorphs == FALSE) {
    return(x.7)
  } else {
    x.8 <- x.7
    x.8 <- x.8[is.polymorphic(x.8, na.omit = TRUE)]
    return(x.8)
  }
}


#' Filter out samples of a vcfR object based on depth
#'
#' @param x A vcfR object
#' @param lower_quant lower quantile to filter by. The top quantile will be \code{1-lower_quant}
#' @param stat Statistic to filter by. Right now, supports "median" to filter on median depth, "mean" to filter on mean depth, or "missingness" to filter on percent missing data.
#' @param filter_bot Omit samples in the bottom quantile
#' @param filter_top Omit samples in the top quantile
#' @param remove_nonpolymorphs After omitting samples, some variants will be nonpolymorphic. These should be removed.
#' @keywords vcf filter
#' @return A filtered vcfR object.
#' @export
#' @examples
#' ultrafiltered_vcf <- samples_filter(filtered_vcf)
samples_filter <- function(x, lower_quant = 0.05, stat = 'median', filter_bot = TRUE, filter_top = FALSE, remove_nonpolymorphs = TRUE) {
  upper_quant <- 1-lower_quant
  dp <- vcfR::extract.gt(x, element = "DP", as.numeric = TRUE)
  if (identical(stat, 'median') == TRUE) {
    median_dps <- apply(dp, MARGIN = 2, median, na.rm = TRUE)
    quants_median <- quantile(median_dps, probs=c(lower_quant, upper_quant))
    # upper quantile
    bot <- median_dps[median_dps < quants_median[1]]
    keep_bot <- !(median_dps %in% bot)
    # lower quantile
    top <- median_dps[median_dps > quants_median[2]]
    keep_top <- !(median_dps %in% top)

    keep_both <- !(median_dps %in% top) & !(median_dps %in% bot)
  } else if (identical(stat, 'mean') == TRUE) {
    mean_dps <- apply(dp, MARGIN = 2, mean, na.rm = TRUE)
    quants_mean <- quantile(mean_dps, probs=c(lower_quant, upper_quant))
    # upper quantile
    bot <- mean_dps[mean_dps < quants_mean[1]]
    keep_bot <- !(mean_dps %in% bot)
    # lower quantile
    top <- mean_dps[mean_dps > quants_mean[2]]
    keep_top <- !(mean_dps %in% top)

    keep_both <- !(mean_dps %in% top) & !(mean_dps %in% bot)
  } else if (identical(stat, 'missingness') == TRUE) {
    myMiss <- apply(dp, MARGIN = 2, function(y){ sum( is.na(y) ) } )
    myMiss <- myMiss / nrow(dp)
    quants_miss <- quantile(myMiss, probs=c(lower_quant, upper_quant))
    # upper quantile
    bot <- myMiss[myMiss < quants_miss[1]]
    keep_bot <- !(myMiss %in% bot)
    # lower quantile
    top <- myMiss[myMiss > quants_miss[2]]
    keep_top <- !(myMiss %in% top)
    keep_both <- !(myMiss %in% top) & !(myMiss %in% bot)
    filter_bot = FALSE
    filter_top = TRUE
  }

  if (filter_bot == TRUE & filter_top == FALSE) {
    x.samp_quant_filtered_bot <- x
    x.samp_quant_filtered_bot@gt <- x.samp_quant_filtered_bot@gt[, c(TRUE, keep_bot)]
    #return(x.samp_quant_filtered_bot)
    x.filtered <- x.samp_quant_filtered_bot
  } else if (filter_bot == FALSE & filter_top == TRUE) {
    x.samp_quant_filtered_top <- x
    x.samp_quant_filtered_top@gt <- x.samp_quant_filtered_top@gt[, c(TRUE, keep_top)]
    #return(x.samp_quant_filtered_top)
    x.filtered <- x.samp_quant_filtered_top
  } else if (filter_bot == TRUE & filter_top == TRUE) {
    x.samp_quant_filtered_both <- x
    x.samp_quant_filtered_both@gt <- x.samp_quant_filtered_both@gt[, c(TRUE, keep_both)]
    #return(x.samp_quant_filtered_both)
    x.filtered <- x.samp_quant_filtered_both
  }
  if (remove_nonpolymorphs == FALSE) {
    return(x.filtered)
  } else {
    x.filtered.polymorphic <- x.filtered
    x.filtered.polymorphic <- x.filtered.polymorphic[vcfR::is.polymorphic(x.filtered.polymorphic, na.omit = TRUE)]
    return(x.filtered.polymorphic)
  }
}

#' Make violin plots for all samples in a vcfR object
#'
#' @param x A vcfR object
#' @param samples_per_row Number of violins per row
#' @param fill_group_aesthetic Fills violin plots based on a vcf statistic. Currently not supported.
#' @keywords vcf filter
#' @return Violin plots of all samples
#' @export
#' @examples
#' plot_violins(ultra_filtered_vcf)
plot_violins <- function(x, samples_per_row = 10, fill_group_aesthetic = NULL) {
  dp <- vcfR::extract.gt(x, element = "DP", as.numeric = TRUE)
  dpf <- reshape2::melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
  dpf_2 <- dpf[ dpf$Depth > 0,]
  # Now make violin plots by sample
  samples <- colnames(x@gt)[-1]
  samps_per_row <- samples_per_row
  myRows <- ceiling(length(levels(dpf_2$Sample))/samps_per_row)
  myList <- vector(mode = "list", length = myRows)

  # Make list of necessary plots
  for(i in 1:myRows){
    myIndex <- c(i*samps_per_row - samps_per_row + 1):c(i*samps_per_row)
    myIndex <- myIndex[myIndex <= length(levels(dpf_2$Sample))]
    myLevels <- levels(dpf_2$Sample)[myIndex]
    myRegex <- paste(myLevels, collapse = "$|^")
    myRegex <- paste("^", myRegex, "$", sep = "")
    myList[[i]] <- dpf_2[grep(myRegex, dpf_2$Sample),]
    myList[[i]]$Sample <- factor(myList[[i]]$Sample)
  }

  # If necessary, calculate heterozygosity
  if (is.null(fill_group_aesthetic) == FALSE) {
    all_het <- vcfR::is_het(extract.gt(x), na_is_false = TRUE)
    myHet <- colSums(all_het)
    nvar <- nrow(x@gt)
    myHet_prop <- myHet/nvar
  }

  # Create the plot.
  myPlots <- vector(mode = "list", length = myRows)
  for(i in 1:myRows){
    if (is.null(fill_group_aesthetic) == TRUE) {
      myPlots[[i]] <- ggplot2::ggplot(myList[[i]], ggplot2::aes(x=Sample, y=Depth)) +
        ggplot2::geom_violin(fill="#8dd3c7", adjust=1.0, scale = "count", trim=TRUE)
    } else {
      myPlots[[i]] <- ggplot2::ggplot(myList[[i]], ggplot2::aes(x=Sample, y=Depth)) +
        ggplot2::geom_violin(fill = myHet_prop, adjust=1.0, scale = "count", trim=TRUE)
    }

    myPlots[[i]] <- myPlots[[i]] + ggplot2::theme_bw()
    myPlots[[i]] <- myPlots[[i]] + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                         axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    myPlots[[i]] <- myPlots[[i]] + ggplot2::scale_y_continuous(trans=scales::log2_trans(),
                                                      breaks=c(1, 10, 100, 800),
                                                      minor_breaks=c(1:10, 2:10*10, 2:8*100))
    myPlots[[i]] <- myPlots[[i]] + ggplot2::theme( panel.grid.major.y=ggplot2::element_line(color = "#A9A9A9", size=0.6) )
    myPlots[[i]] <- myPlots[[i]] + ggplot2::theme( panel.grid.minor.y=ggplot2::element_line(color = "#C0C0C0", size=0.2) )
  }
  cowplot::plot_grid(plotlist = myPlots, nrow = myRows)
}
