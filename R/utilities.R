#' Calculate global FST
#'
#' Calculates FST globally across all subsets
#'
#' @param x allele-count and observed heterozygosity data. Required columns are
#'   \code{variable} containing locus IDs, \code{vf} containing
#'   population/facet IDs, two containing allele counts for
#'   the major/ref and minor/alt loci, respectively, and \code{ho} for observed
#'   heterozygosity. Note that this object can be created from 0/1/2 formatted
#'   genotypic data with \code{prep_as_from_sn}.
#'
#' @param ac_cols character, default c("0", "1"). Column names for the two
#'   alleles in x.
#'
#' @author William Hemstrom
#' @author Andy Lee
#'
#' @examples
#' #global_fst(mon_sn)
#'
#' @export
global_fst <- function(x, ac_cols = c("0", "1")){
  ac_cols <- ac_cols[which(ac_cols %in% colnames(x))]
  nt <- data.table::dcast(x[,rowSums(.SD), .SDcols = ac_cols, by = .(vf, variable)], variable ~ vf, value.var = "V1")
  psm <- x[,.SD/rowSums(.SD), .SDcols = ac_cols, by = .(vf, variable)]
  ntotm <- nt[,-1]/2
  ntotm <- cbind(nt[,1], ntotm)

  hom <- data.table::dcast(x, variable ~ vf, value.var = "ho")

  pops <- unique(x[,vf])
  r <- length(pops) # number of comps
  nbar <- rowMeans(ntotm[,-1]) #average sample size in individuals
  CV <- matrixStats::rowSds(as.matrix(ntotm[,-1]))/nbar # coefficient of variation in sample size
  nc <- nbar*(1-(CV^2)/r)
  parts <- vector("list", length(ac_cols))

  out <- data.table::data.table(variable = sort(unique(x$variable)),
                                fst = 0,
                                a = 0,
                                b = 0,
                                c = 0)

  for(k in 1:length(ac_cols)){
    psf_m <- dcast(psm, variable ~ vf, value.var = colnames(psm)[k + 2])

    # need to determine per-allele hom
    # if(!bi_allelic){
    #   het_cols_containing_k <- which(xor(gc_colnames_1 == colnames(ps1_f)[k], gc_colnames_2 == colnames(ps1_f)[k])) # hets (xor) with this allele
    #   tiho <- .fix..call(rowSums(idat[,..gc_cols][,..het_cols_containing_k])/intot)
    #   tjho <- .fix..call(rowSums(jdat[,..gc_cols][,..het_cols_containing_k])/jntot)
    # }

    # otherwise we have the ho already
    absent <- which(rowSums(psf_m[,-1]) == 0)
    thom <- data.table::copy(hom)
    if(length(absent) != 0){
      data.table::set(thom, absent, 2:ncol(thom), value = 0)
    }

    parts[[k]] <- .per_all_f_stat_components(ntotm = ntotm[,-1], psm = psf_m[,-1], r = r, nbar = nbar, nc = nc, hom = thom[,-1])
  }

  a <- matrix(unlist(purrr::map(parts, "a")), ncol = length(parts))
  b <- matrix(unlist(purrr::map(parts, "b")), ncol = length(parts))
  c <- matrix(unlist(purrr::map(parts, "c")), ncol = length(parts))

  data.table::set(out, j = "a", value = rowSums(a)) # write a
  data.table::set(out, j = "b", value = rowSums(b)) # write b
  data.table::set(out, j = "c", value = rowSums(c)) # write v
  # Fst <- rowSums(a)/rowSums(a + b + c)

  out[,fst := a/(a + b + c)]
  meta_colnames <- colnames(x)[1:(which(colnames(x) %in% ac_cols)[1] - 1)]
  out$nk <- rowSums(nt[,-1])


  means <- out[,list(weighted_fst = stats::weighted.mean(a, w = nk, na.rm = T)/
                       stats::weighted.mean(a + b + c, w = nk, na.rm = T),
                     fst = mean(a, na.rm = TRUE)/mean(a + b + c, na.rm = TRUE))]

  return(list(means = means, pairwise = out))
}

#' smartPCAs
#' Generate a raw smartPCA

#' @param sn allele-count and observed heterozygosity data. Required columns are
#'   \code{variable} containing locus IDs, \code{vf} containing
#'   population/facet IDs, two containing allele counts for
#'   the major/ref and minor/alt loci, respectively, and \code{ho} for observed
#'   heterozygosity. Note that this object can be created from 0/1/2 formatted
#'   genotypic data with \code{prep_as_from_sn}.
#'
#' @author William Hemstrom
#' @author Andy Lee
#'
#' @export
#'
#' @examples
#' # example code, generate a smartPCA
#' # quick_smartPCA(mon_sn)

quick_smartPCA <- function(sn){
  ms <- colMeans(sn, na.rm = TRUE)
  ms <- ms[col(sn)]
  afs <- (1+ms)/(2 + 2*nrow(sn)) # adjusted eqn, from Price et al 2006
  # afs <- ms/2 # according to patterson
  sn <- (sn - ms)/sqrt(afs*(1-afs)) # eqn 3, patterson et al 2006
  sn[is.na(sn)] <- 0

  pca_r <- svd(sn)
  colnames(pca_r$u) <- paste0("PC", 1:ncol(pca_r$u))

  return(pca_r$u)
}

#' Generate summary statistics

#' @param as allele-count and observed heterozygosity data. Required columns are
#'   \code{variable} containing locus IDs, \code{vf} containing
#'   population/facet IDs, two containing allele counts for
#'   the major/ref and minor/alt loci, respectively, and \code{ho} for observed
#'   heterozygosity. Note that this object can be created from 0/1/2 formatted
#'   genotypic data with \code{prep_as_from_sn}.
#' @param genotypes genotypic data, formatted with genotypes as 0/1/2 for the
#'   homozygous ref/major, heterozygous, and homozygous alt/minor, respectively.
#'   Each row must be one locus, each column one individual. Column names must
#'   be unique.
#' @param facet character vector noting population informtion for each
#'   inidivudal sample.
#' @param ac_cols character, default c("0", "1"). Column names for the two
#'   alleles in 'as' format.
#' @param fst_cut numeric, default 0.95. The Fst quantile above which to label
#'   loci as high-Fst.
#' @param store_pca logical, default FALSE. If TRUE, PCAs for both all and high-
#'   Fst loci will be retained and returned.
#' @references Jombart, T., Devillard, S. & Balloux, F. Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. BMC Genet 11, 94 (2010). https://doi.org/10.1186
#'
#' @author William Hemstrom
#' @author Andy Lee
#'
#' @export
#'
#' @examples
#' #Generate summary statistics
#' generate_summary_stats(mon_sn, genotypes, facet="population", fst_cut=.95, store_pca = FALSE)
generate_summary_stats <- function(as, genotypes, facet, ac_cols = c("0", "1"), fst_cut = .95,
                                   store_pca = FALSE){
  fst <- global_fst(as, ac_cols)
  ofst <- fst$means$fst

  opca <- quick_smartPCA(t(genotypes))
  opca <- as.data.table(opca)
  opca$pop <- facet
  omav <- manova(cbind(PC1, PC2) ~ pop, opca)
  omav <- summary(omav)$stats[1,"approx F"]


  high_fst <- which(fst$pairwise$fst >= quantile(fst$pairwise$fst, fst_cut, na.rm = TRUE))
  high_fst <- fst$pairwise$variable[high_fst]

  high_fst_as <- as[which(as$variable %in% high_fst),]
  fst_high <- global_fst(high_fst_as, ac_cols)
  tfst <- fst_high$means$fst


  tpca <- quick_smartPCA(t(genotypes[which(colnames(genotypes) %in% high_fst),]))

  tpca <- as.data.table(tpca)
  tpca$pop <- facet

  mav <- manova(cbind(PC1, PC2) ~ pop, tpca)
  mav <- summary(mav)$stats[1,"approx F"]

  if(store_pca){
    return(list(Fstat = mav, fst = tfst, delta_Fstat = mav - omav, delta_fst = tfst - ofst,
                init_Fstat = omav, init_fst = ofst,
                pca = list(all_vars = opca, selected = tpca)))
  }
  return(list(Fstat = mav, fst = tfst, delta_Fstat = mav - omav, delta_fst = tfst - ofst, init_Fstat = omav, init_fst = ofst))
}

#' Calculate p-value for the observed PCA
#'
#' @param observed the observed dataset
#' @param null the null dataset
#'
#' @author William Hemstrom
#' @author Andy Lee
#'
#' @export
#' @examples # get_p_values(observed, null)
#'
get_p_values <- function(observed, null, h0 = c("less", "greater")){
  ec_dists <- lapply(null, ecdf)

  p <- numeric(length(observed))
  names(p) <- names(observed)
  for(i in 1:length(ec_dists)){
    p[i] <- ec_dists[[i]](observed[[i]])
    if(h0[i] == "greater"){
      p[i] <- 1-p[i]
    }
  }
  return(p)
}


#' Prepare data in 'as' format from data in 'sn' format
#' @param x genotypic data, formatted with genotypes as 0/1/2 for the
#'   homozygous ref/major, heterozygous, and homozygous alt/minor, respectively.
#'   Each row must be one locus, each column one individual. Column names must
#'   be unique.
#' @param facet character vector noting population inform for each
#'   individual sample.
#' @author William Hemstrom
#' @author Andy Lee
#'
#' @export
#' @examples # prep_as_from_sn(mon_sn, facet="populations" )

prep_as_from_sn <- function(x, facet){
  browser()
  .tab_func <- function(x){
    x <- data.table::melt(data.table::transpose(x, keep.names = "samp"), id.vars = "samp") # transpose and melt
    data.table::set(x, j = "vf", value = facet[as.numeric(x$samp)])

    gmat <- data.table::dcast(data.table::setDT(x), variable + vf ~ value, value.var='value', length) # cast

    ho <- gmat$`1`/(gmat$`0` + gmat$`2`)

    amat <- gmat[,.(variable, vf)]
    amat$`0` <- gmat$`0`*2 + gmat$`1`
    amat$`1` <- gmat$`2`*2 + gmat$`1`

    gmat$`NA` <- NULL
  }


  # loop to save memory if large
  max_rows <- 100000000 # about 1.5G of storage for the big melt
  max_snps <- ceiling(max_rows/ncol(x)) # max snps tolerable at once
  n_iters <- ceiling(nrow(x)/max_snps)
  titer <- 1


  colnames(x) <- as.character(1:ncol(x))





  return(list(gmat = gmat, amat = amat, ho = ho))
}
