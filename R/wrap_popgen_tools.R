#' Make DAPC plots
#'
#' @param x A genclone object
#' @param formula Hierarchical formula to group the populations by in plots
#' @param name Used to label plot axes and titles
#' @keywords genclone genind dapc
#' @return Three DAPC plots
#' @export
#' @examples
#' dapc_workflow(ultra_filtered_genclone, ~Nursery, "Nurseries as Pops")
dapc_workflow <- function(x, formula, name) {
  adegenet::setPop(x) <- formula

  # First the scatter plot
  # Get optimum number of PC axes
  dapc.genclone <- adegenet::dapc(x, var.contrib = TRUE, center = FALSE, n.pca = 40,  n.da = nPop(x) - 1)
  optimization <- adegenet::optim.a.score(dapc.genclone)
  # Repeat DAPC calculation with correct number of axes
  dapc.genclone <- adegenet::dapc(x, var.contrib = TRUE, center = FALSE, n.pca = optimization$best,  n.da = nPop(x) - 1)
  p1 <- ade4::scatter(dapc.genclone, mstree = TRUE, lwd = 2, lty = 2, main = paste("DAPC of", name), col = rainbow(4))

  # Now compoplot
  # Melt DAPC results into a data frame
  dapc.results <- as.data.frame(dapc.genclone$posterior)
  dapc.results$pop <- adegenet::pop(x)
  dapc.results$indNames <- rownames(dapc.results)
  dapc.results <- reshape2::melt(dapc.results)
  colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

  # Employ sorting method used in Structure plots
  dapc_sorted <- dapc.results %>%
    dplyr::group_by(Original_Pop, Assigned_Pop) %>%
    dplyr::summarise(avg = mean(Posterior_membership_probability)) %>%
    dplyr::right_join(dapc.results, by = c("Original_Pop", "Assigned_Pop")) %>%
    dplyr::arrange(Original_Pop, desc(avg), desc(Posterior_membership_probability)) %>%
    dplyr::select(-avg) %>%
    dplyr::ungroup
  #compoplot.dapc(dapc_sorted)

  # Plot the plot
  p <- ggplot2::ggplot(dapc_sorted, aes(x=fct_inorder(Sample), y=Posterior_membership_probability, fill=Assigned_Pop))
  p <- p + ggplot2::geom_bar(stat='identity', width = 1)
  p <- p + ggplot2::scale_fill_manual(values = rainbow(4))
  p <- p + ggplot2::facet_grid(~Original_Pop, scales = "free_x")
  p <- p + ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
  #p <- p + xlab(paste("Samples of", name))
  p <- p + ggplot2::ylab("Posterior membership probability")
  p
  p3 <- ade4::scatter(dapc.genclone, 1, 1, bg = "white", scree.da=FALSE, legend = TRUE, solid = .4, col = rainbow(4))
  return(list(p1, p3, p))
}

#' Make PCA plots
#'
#' @param x A genlight object
#' @param pop_def Define what populations are to be used
#' @param pop_map Reference population map/strata
#' @param filter_strata Filter down your strata to only samples used
#' @keywords genlight pca
#' @return PCA plots
#' @export
#' @examples
#' pca_workflow(ultra_filtered_genlight, 'Nursery', strata_defs)
pca_workflow <- function(x, pop_def, pop_map, filter_strata = TRUE) {

  # first determine how to set populations
  if ((is.null(pop_def) == FALSE) & (is.null(pop_map) == FALSE) & (filter_strata == TRUE)) {
    pop_map_filtered <- subset(pop_map[,], pop_map$Sample %in% (x@ind.names))
    adegenet::pop(x) <- pop_map_filtered[[pop_def]]
    x
  } else if ((is.null(pop_def) == FALSE) & (is.null(pop_map) == FALSE) & (filter_strata == FALSE)) {
    adegenet::pop(x) <- pop_map[[pop_def]]
  }
  # the else condition is that populations are already supplied to the genind...
  x.pca <- adegenet::glPca(x, nf = 3)
  barplot(x.pca$eig, col = heat.colors(50), main = "PCA Eigenvalues")

  iPst.pca.scores <- as.data.frame(x.pca$scores)
  iPst.pca.scores$pop <- adegenet::pop(x)
  iPst.pca.scores$sample <- rownames(iPst.pca.scores)
  #return(iPst.pca.scores)
  set.seed(9)
  p12 <- ggplot2::ggplot(iPst.pca.scores, aes(x=PC1, y=PC2, colour=pop))
  p12 <- p12 + ggplot2::geom_point(size=2)
  p12 <- p12 + ggplot2::geom_text(aes(label=sample), hjust = 0, vjust = 0)
  p12 <- p12 + ggplot2::stat_ellipse(level = 0.95, size = 1)
  p12 <- p12 + ggplot2::scale_color_manual(values = rainbow(nPop(x)))
  p12 <- p12 + ggplot2::geom_hline(yintercept = 0)
  p12 <- p12 + ggplot2::geom_vline(xintercept = 0)
  p12 <- p12 + ggplot2::theme_bw()
  p12

  p13 <- ggplot2::ggplot(iPst.pca.scores, aes(x=PC1, y=PC3, colour=pop))
  p13 <- p13 + ggplot2::geom_point(size=2)
  p13 <- p13 + ggplot2::geom_text(aes(label=sample), hjust = 0, vjust = 0)
  p13 <- p13 + ggplot2::stat_ellipse(level = 0.95, size = 1)
  p13 <- p13 + ggplot2::scale_color_manual(values = rainbow(nPop(x)))
  p13 <- p13 + ggplot2::geom_hline(yintercept = 0)
  p13 <- p13 + ggplot2::geom_vline(xintercept = 0)
  p13 <- p13 + ggplot2::theme_bw()
  p13
  return(list(p12, p13))
}


#' Make phylogenetic trees colored according to strata
#'
#' @param x Preferably a genclone object, genlights sometimes break it. Must have strata
#' @param formula Hierarchical formula to group the populations by in plots. Do not wrap in quotes... I think.
#' @param name Used to label plot axes and titles.
#' @param scale_quality What to color the tree tips by. Use 'population' to color by any "formula" you set. Other options are: 'missing', 'median_depth', 'mean_depth', 'heterozygosity', and 'allele_depth' 
#' @param vcf Original vcfR object, used to calculate statistics on. Be sure it's the same exact individuals and variants as the data in `x`.
#' @param geo_tree Optionally supply your own tree instead of running poppr::aboot
#' @param treeType Algorithm to make tree by. Make sure the package to use is loaded.
#' @param dist Algorithm to calculate distance. Make sure the package to use is loaded.
#' @param input_type Specify if a genclone or genlight is inputted
#' @keywords genclone genind phylogeny ape
#' @return Phylogenetic tree plot
#' @export
#' @examples
#' tree_workflow(ultra_filtered_genlight, ~Nursery, "Nursery of origin", "population", ultra_filtered_vcf)
tree_workflow <- function(x, formula = NULL, name, scale_quality = 'population', vcf = NULL, geo_tree = NULL, treeType = 'upgma', dist = 'provesti.dist', input_type = 'genclone') {
  if (is.null(formula) == FALSE) {
    setPop(x) <- formula
  }
  
  # if bitwise.dist is used, object must all be one type of ploid
  # so remove all samples that aren't that one ploid
  if (identical(dist, 'bitwise.dist') == TRUE) {
    filter_ploidy <- function(x) {
      # only run if ploidy really is variable
      if (identical(min(adegenet::ploidy(x)), max(adegenet::ploidy(x))) == FALSE) {
        # find most common ploid
        x.ploidy <- adegenet::ploidy(x)
        ploidy_mode <- as.integer(names(sort(summary(as.factor(x.ploidy)), decreasing = T)[1]))
        x.remove <- x.ploidy[x.ploidy != ploidy_mode]
        # Remove samples that don't match most common ploidy
        x <- x[adegenet::indNames(x) != names(x.remove)]
      }
    }
    x <- filter_ploidy(x)
  }
  # Set legend and color properties based on scale type
  if (identical(scale_quality, 'population') == FALSE) {
    # Until I get actual continuous legends figured out...
    discrete_legend <- FALSE
    measure <- adegenet::pop(x)
    legend_measure <- levels(measure)
    # Make legend
    #discrete_legend <- FALSE
    #measure <- pop(x)
    #legend_measure <- levels(measure)
    # Make color palette and
    # assign numbers to colors assuming input is raw vector of numbers
    col_scale_num <- 20
    min_color <- '#fee0d2'
    max_color <- '#8b0000'
    scale_to_cols <- function(raw_vector) {
      col_scale_num <- 20
      min_color <- '#fee0d2'
      max_color <- '#8b0000'
      color_palette <- colorRampPalette(colors=c(min_color, max_color))(col_scale_num)
      # myMiss_rescale <- round(rescale(myMiss, to = c(0, col_scale_num)))
      # cols <- color_palette[myMiss_rescale]
      bins <- cut(raw_vector, breaks = seq(min(raw_vector), max(raw_vector), len = col_scale_num), include.lowest = TRUE)
      color_palette[bins]
    }
  } else if (identical(scale_quality, 'population') == TRUE) {
    discrete_legend <- TRUE
    measure <- adegenet::pop(x)
    legend_measure <- levels(measure)
    
  }
  
  # Create actual scales based on inputted parameter
  if (identical(scale_quality, 'population') == TRUE) {
    #dark_color_palette <- brewer.pal(n = nPop(x), name = "Dark2")
    #cols <- dark_color_palette[pop(x)]
    rainbow_pal <- rainbow(length(levels(adegenet::pop(x))))
    cols <- rainbow_pal[adegenet::pop(x)]
  } else if (identical(scale_quality, 'missing') == TRUE) {
    # Assign sample to color based on missingness
    dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)
    myMiss <- apply(dp, MARGIN = 2, function(y){ sum( is.na(y) ) } )
    myMiss <- myMiss / nrow(dp)
    print(myMiss)
    cols <- scale_to_cols(myMiss)
    explanation <- "Proportion Missing Data"
    legend2_min <- min(myMiss)
    legend2_max <- max(myMiss)
  } else if (identical(scale_quality, 'median_depth') == TRUE) {
    dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)
    myDepth <- apply(dp, MARGIN = 2, function(y){ median(y, na.rm = TRUE) })
    print(myDepth)
    cols <- scale_to_cols(myDepth)
    explanation <- "Median variant depth"
    legend2_min <- min(myDepth)
    legend2_max <- max(myDepth)
  } else if (identical(scale_quality, 'mean_depth') == TRUE) {
    dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)
    myDepth <- apply(dp, MARGIN = 2, function(y){ mean(y, na.rm = TRUE) })
    print(myDepth)
    cols <- scale_to_cols(myDepth)
    explanation <- "Mean variant depth"
    legend2_min <- min(myDepth)
    legend2_max <- max(myDepth)
  } else if (identical(scale_quality, 'heterozygosity') == TRUE) {
    myHet <- colSums(vcfR::is_het(extract.gt(vcf), na_is_false = TRUE))
    print(myHet)
    cols <- scale_to_cols(myHet)
    explanation <- "Heterozygous loci"
    legend2_min <- min(myHet)
    legend2_max <- max(myHet)
  } else if (identical(scale_quality, 'allele_depth') == TRUE) {
    ad <- vcfR::extract.gt(vcf, element = 'AD', as.numeric = TRUE)
    myAlleleDepth <- apply(ad, MARGIN = 2, function(x){ median(x, na.rm = TRUE)})
    print(myAlleleDepth)
    cols <- scale_to_cols(myAlleleDepth)
    explanation <- "Median allele depth"
    legend2_min <- min(myAlleleDepth)
    legend2_max <- max(myAlleleDepth)
  } else {
    warnings(paste(c("scale_quality", scale_quality, "not found.")))
    stop()
  }
  
  # Make tree
  # works for genclone
  # geo_tree <- aboot(x, dist = provesti.dist, sample = 50, tree = treeType, missing = "ignore",
  #                   quiet = TRUE, showtree = FALSE)
  # must use bitwise.dist for genlight. Let user specify distance using elipsis?
  geo_tree <- poppr::aboot(x, sample = 100, tree = treeType, missing = "ignore",
                    quiet = TRUE, showtree = FALSE, dist = dist)
  #edge_color_palette <- brewer.pal(n = nPop(x), name = "Dark2")
  edge_color_palette <- rainbow_pal <- rainbow(length(levels(pop(x))))
  ape::plot.phylo(ape::ladderize(phangorn::midpoint(geo_tree)), cex = 0.8, font = 2, adj = 0,
             label.offset = 0.0125, underscore = TRUE, #align.tip.label = 2,
             #edge.color = rainbow(nPop(x))[pop(x)], tip.color = cols)
             #edge.color = edge.cols, tip.color = cols)
             #edge.color = edge_color_palette[pop(x)],
             tip.color = cols)
  
  nodelabels(geo_tree$node.label,  adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
  axis(side = 1)
  
  # Plot legend
  # since I can't get edge colors to be accurate...
  if (discrete_legend == TRUE) {
    legend1_fill_measure <- edge_color_palette
    legend1_name <- deparse(substitute((formula)))
    legend('right', legend = legend_measure, fill = legend1_fill_measure, title = legend1_name,
           border = FALSE, bty = "n", cex = 1)
  }
  # Add second legend if continuous node color scaling
  if (discrete_legend == FALSE) {
    # Continuous legend
    #x.coord = c(-20, -15, -15, -20)
    #y.coord = c(0, 60, 60, 0)
    #x.coord = c(0, 0, 0, 0)
    #y.coord = c(0, 0, 0, 0)
    #pnts = cbind(x =c(0.8,0.9,0.9,0.8), y =c(1.0,1.0,0.8,0.8))
    #legend.gradient(cbind(x.coord = x.coord - 150, y.coord = y.coord - 30),
    #legend.gradient(pnts, cols = c('#fee0d2', '#8b0000'), title = "TITLE", limits = c(min(myMiss), max(myMiss)))
    #legend('topleft', legend = legend_measure, fill = rainbow(nPop(x)), title = "Population (branches)",
    #       border = FALSE, bty = "n", cex = 1)
    # For now, just put on another legend with the discrete upper and lower bounds of color
    legend('left', title = paste(explanation, "(tips)"),
           legend = c(round(legend2_min, 3), round(legend2_max, 3)),
           fill = c(min_color, max_color),
           border = FALSE, bty = "n", cex = 1, ncol = 2)
  }
  
  # Label tree title
  if(missing(name)){
    
    title(xlab = paste("Genetic distance, tips colored by", scale_quality))
  } else {
    title(xlab = paste("Genetic distance using", name, "; tips colored by", scale_quality))
  }
}
