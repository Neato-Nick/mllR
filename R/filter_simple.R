filter_vcf <- function (x, min_dp = 4, max_dp = 1000, min_mq = 20, samples_miss1 = 0.999, samples_miss2 = 0.7, vars_miss = 0.2, top_dp_quantile_filter = TRUE, remove_nonpolymorphs = TRUE) {
  dp <- extract.gt(x, element = "DP", as.numeric = TRUE)

  # Censor variants depth < min_dp or > max_dp
  dp[dp < min_dp] <- NA
  dp[dp > max_dp] <- NA
  x.2 <- x
  x.2@gt[,-1][ is.na(dp) == TRUE] <- NA

  # Censor variants in top 95% quantile of depth
  if (top_dp_quantile_filter == TRUE) {
    dp <- extract.gt(x.2, element = "DP", as.numeric = TRUE)
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
  mq <- extract.info(genotype_calls_missing, element = "MQ", as.numeric = TRUE)
  # bypass mapping quality filter if there are missing values in MQ
  # indicates denovo assembly
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
    # } else {
    #   x.maf <- as.data.frame(maf(x.6))
    #   x.maf$Variant <- seq(1:length(x.maf$Count))
    #   x.maf <- subset(x.maf, subset = x.maf$Count > 0)
    #   x.7 <- x.6
    #   x.7 <- x.7[x.maf$Variant,]
    #   x.7
  } else {
    x.8 <- x.7
    x.8 <- x.8[is.polymorphic(x.8, na.omit = TRUE)]
    return(x.8)
  }
}
