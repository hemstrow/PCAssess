#' Run PCA permutation pipeline
#'
#' Runs the full PCA boostrap/permutation testing pipeline using genetic data.
#'
#' @param x Genetic data formatted as single numbers indicated genotype, where
#'   0 and 2 are homozygoytes and 1 is a heterozygote. Rows are SNPs and columns
#'   are individuals.
#' @param facet character. Categorical variable by which the data was broken up e.g., population, treatment.
#' @param n numeric. Number of permutations to run.
#' @param fst_cut numeric, default 0.95. Cut-off used to select highest fst loci.
#' @param par numeric. Number of parallel threads to use.
#' @param store_pca logical, default FALSE. If TRUE, returns the raw PCA.
#'
#'
#' @return a list containing: null distribution of the F-statistic observed values, the observed F-statistic, and the p-values from PCA permutation testing
#'
#' @references Jombart, T., Devillard, S. & Balloux, F. Discriminant analysis of
#'  principal components: a new method for the analysis of genetically
#'  structured populations. BMC Genet 11, 94 (2010).
#'  \url{https://doi.org/10.1186}
#'
#' @author William Hemstrom
#' @author Andy Lee
#'
#' @export
#' @examples
#' run_permutation(snps_dat, "pop", 1000, fst_cut, par = FALSE, store_pca = FALSE)

run_permutation <- function(x, facet, n, fst_cut = .95, par = FALSE, store_pca = FALSE){
  #===============sanity checks===================
  msg <- character()
  warn <- character()


  if(fst_cut > 1 | fst_cut < 0){
    msg <- c(msg, "Fst cut-off must be between 0 and 1.\n")
  }


  if(length(msg) > 0){
    stop(msg)
  }

  if(length(warn) > 0){
    warning(warn)
  }

  #===============run the permutation=============
  x_as <- .prep_boots(x, facet, n)
  real_x <- generate_summary_stats(as = x_as$real_as, genotypes = x, facet = facet, fst_cut = fst_cut, store_pca = store_pca)
  boots_x <- .do_boots(x_as$boot_as, genotypes = x, par = par, fst_cut = fst_cut, store_pca = store_pca)

  if(store_pca){
    p_x <- get_p_values(real_x[which(names(real_x) != "pca")], boots_x$res[,-1],
                        c("greater", "greater", "greater", "greater", "greater", "greater"))
    return(list(null_distribution = boots_x$res, observed_values = real_x[which(names(real_x) != "pca")], pvalues = p_x,
                real_pca = real_x$pca,
                null_pca = boots_x$pca))
  }
  else{
    p_x <- get_p_values(real_x, boots_x[,-1], c("greater", "greater", "greater", "greater", "greater", "greater"))
    return(list(null_distribution = boots_x, observed_values = real_x, pvalues = p_x))
  }
}

#' Plot PCA permutation results
#'
#' Plots results from PCA permutation testing
#'
#' @param permute_res a list containing the null distribution, observed values, and p-vaules from PCA permutation testing
#' @param n_permute_pcas number of permutation to plot
#'
#' @return a ggplot of the null distribution of F-statistic and observed F-statistic from PCA permutation testing.
#'  if do_pcas is \code{TRUE}, also returns a ggplot of the observed PCA and n number of permuted PCAs
#'
#'
#' @author William Hemstrom
#' @author Andy Lee
#'
#' @export
#'
#' @examples
#' # plot observed PCA and 100 permutations
#' plot_permutation_res(permute_res, 100)
#'
#' # plot only the observed PCA
#' plot_permutation_res(permute_res, 0)
#'
plot_permutation_res <- function(permute_res, n_boot_pcas = 1){

  #=======result distribution plot, always make this=======
  dist_plot <- ggplot2::ggplot(permute_res$null_distribution, ggplot2::aes(x = delta_Fstat)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = permute_res$observed_values$delta_Fstat, color = "red") +
    ggplot2::geom_label(data = data.frame(x = permute_res$observed_values$delta_Fstat,
                                         lab = paste0("p = ", permute_res$pvalues["delta_Fstat"]),
                                         y = max(density(permute_res$null_distribution$delta_Fstat)$y)*1.1),
                       ggplot2::aes(x = x, y = y, label = lab), hjust = "inward") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw() +
    ggplot2::xlab(bquote(Delta * "F (change in clustering)")) +
    ggplot2::ylab("Density")


  do_pcas <- "real_pca" %in% names(permute_res)
  if(!do_pcas){
    return(dist_plot)
  }

  #==========PCA plots if present=========
  plot_pca <- function(all, selected){
    cols <- c("PC1", "PC2", colnames(all)[-grep("PC[1-9]+", colnames(all))])
    pcad <- rbindlist(list(all = .fix..call(all[,..cols]), selected = .fix..call(selected[,..cols])), idcol = "type")

    ccol <- cols[3]
    colnames(pcad)[4] <- "ccol"
    p <- ggplot2::ggplot(pcad, ggplot2::aes(x = PC1, y = PC2, color = ccol)) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(type~., ) +
      ggplot2::theme_bw() +
      ggplot2::guides(color = ggplot2::guide_legend(title = ccol))  +
      khroma::scale_color_batlow(discrete = TRUE) +
      ggplot2::theme(axis.text = element_blank(),
                     axis.ticks = element_blank())

    return(p)
  }

  # real
  pca_real <- plot_pca(permute_res$real_pca$all_vars, permute_res$real_pca$selected)

  # boots
  if(n_boot_pcas == 0){
    return(cowplot::plot_grid(pca_real, dist_plot, nrow = 2, align = "hv", axis = "trbl"))
  }

  bp <- vector("list", n_boot_pcas)
  use <- sample(length(permute_res$null_pca), n_boot_pcas, replace = FALSE)
  for(i in 1:n_boot_pcas){
    bp[[i]] <- plot_pca(permute_res$null_pca[[use[i]]]$all_vars, permute_res$null_pca[[use[i]]]$selected) + ggplot2::guides(color = "none")
  }
  top <- cowplot::plot_grid(plotlist = c(list(pca_real +
                                                ggplot2::guides(color = "none") +
                                                ggplot2::theme(plot.background = element_rect(color = "lightgrey", fill = "lightgrey"))),
                                         bp),
                            nrow = length(bp) + 1)
  top <- cowplot::plot_grid(top, ggpubr::get_legend(pca_real), nrow = 1, rel_widths = c(.9, .1))
  return(cowplot::plot_grid(top, dist_plot, nrow = 2, rel_heights = c(1*((n_boot_pcas + 1)/(n_boot_pcas + 2)),
                                                                      1*(1/(n_boot_pcas + 2)))))

}
