global_fst <- function(x, ac_cols = c("A", "C", "G", "T")){
  x <- as.data.table(x)
  ac_cols <- ac_cols[which(ac_cols %in% colnames(x))]
  nt <- data.table::dcast(x[,rowSums(.SD), .SDcols = ac_cols, by = .(subfacet, .snp.id)], .snp.id ~ subfacet, value.var = "V1")
  psm <- x[,.SD/rowSums(.SD), .SDcols = ac_cols, by = .(subfacet, .snp.id)]
  ntotm <- nt/2
  ntotm[,1] <- nt[,1]

  hom <- data.table::dcast(x, .snp.id ~ subfacet, value.var = "ho")

  pops <- unique(x[,subfacet])
  r <- length(pops) # number of comps
  nbar <- rowMeans(ntotm[,-1]) #average sample size in individuals
  CV <- matrixStats::rowSds(as.matrix(ntotm[,-1]))/nbar # coefficient of variation in sample size
  nc <- nbar*(1-(CV^2)/r)
  parts <- vector("list", length(ac_cols))

  out <- data.table::data.table(.snp.id = sort(unique(x$.snp.id)),
                                fst = 0,
                                a = 0,
                                b = 0,
                                c = 0)

  for(k in 1:length(ac_cols)){
    psf_m <- dcast(psm, .snp.id ~ subfacet, value.var = colnames(psm)[k + 2])

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

    parts[[k]] <- snpR:::.per_all_f_stat_components(ntotm = ntotm[,-1], psm = psf_m[,-1], r = r, nbar = nbar, nc = nc, hom = thom[,-1])
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
  out <- merge(out, snpR:::.fix..call(x[subfacet == pops[1],..meta_colnames]),
               by = ".snp.id")


  out <- cbind(comparison = ".GLOBAL", out)
  out$nk <- rowSums(nt[,-1])

  col.ord <- c("comparison", meta_colnames, "fst", "a", "b", "c", "nk")
  out <- snpR:::.fix..call(out[,..col.ord])

  means <- out[,list(weighted_fst = stats::weighted.mean(a, w = nk, na.rm = T)/
                       stats::weighted.mean(a + b + c, w = nk, na.rm = T),
                     fst = mean(a, na.rm = TRUE)/mean(a + b + c, na.rm = TRUE)),
               by = subfacet]

  return(list(means = means, pairwise = out))
}

boot_as <- function(x, n, facet, ret_gs = FALSE, by = "sample"){
  ..tm <- ..ord <-  NULL

  out <- vector("list", n)
  opts <- snpR:::.get.task.list(x, facet)
  shuff_ord <- out
  for(i in 1:n){

    # get the maf identities for the base facet
    shuff_ord[[i]] <- sample(ncol(x), ncol(x), replace = TRUE)
    shuff <- x
    sample.meta(shuff)[,facet] <- sample.meta(shuff)[shuff_ord[[i]],facet]
    tas <- vector("list", nrow(opts))

    # do for each facet level.
    for(j in 1:nrow(opts)){
      tm <- snpR:::.fetch.sample.meta.matching.task.list(shuff, opts[j,])
      ttab <- snpR:::.tabulate_genotypes(snpR:::.fix..call(as.data.table(genotypes(shuff))[,..tm]), x@mDat)
      tas[[j]] <- data.table::as.data.table(as.matrix(ttab$as))

      # fix missing columns (no genotypes with this allele in this subfacet)
      missing <- which(!colnames(x@geno.tables$as) %in% colnames(tas[[j]]))
      if(length(missing) != 0){
        fill <- matrix(0, nrow(tas[[j]]), length(missing))
        fill <- as.data.table(fill)
        colnames(fill) <- colnames(x@geno.tables$as)[missing]
        tas[[j]] <- cbind(tas[[j]], fill)
      }
      ord <- colnames(x@geno.tables$as)
      tas[[j]] <- snpR:::.fix..call(tas[[j]][,..ord])

      if(ret_gs){
        gst <- as.data.table(as.matrix(ttab$gs))

        # fix missing columns (no genotypes with this allele in this subfacet)
        missing <- which(!colnames(x@geno.tables$gs) %in% colnames(gst))
        if(length(missing) != 0){
          fill <- matrix(0, nrow(gst), length(missing))
          fill <- as.data.table(fill)
          colnames(fill) <- colnames(x@geno.tables$gs)[missing]
          gst <- cbind(gst, fill)
        }
        ord <- colnames(x@geno.tables$gs)
        gst<- snpR:::.fix..call(gst[,..ord])

        tas[[j]] <- cbind(tas[[j]], gst)
      }

      tas[[j]]$ho <- snpR:::.ho_func(ttab, x@snp.form)
      tas[[j]]$.snp.id <- x@snp.meta$.snp.id
      tas[[j]]$subfacet <- opts[j,2]
      tas[[j]]$facet <- facet

    }
    tas <- data.table::rbindlist(tas)
    tas[is.na(tas)] <- 0
    ord <- c((ncol(tas) - 3):ncol(tas), 1:(ncol(tas) - 4))
    out[[i]] <- snpR:::.fix..call(tas[,..ord])
    out[[i]] <- snpR:::.suppress_specific_warning(cbind(snp.meta(x)[,which(colnames(snp.meta(x)) != ".snp.id"),drop=FALSE], out[[i]]), "row names were found from a short variable")
  }

  return(list(shuff_ord = shuff_ord, out = out))
}

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

generate_summary_stats <- function(as, genotypes, facet, ac_cols = c("A", "C", "G", "T"), fst_cut = .95,
                                   store_pca = FALSE){
  fst <- global_fst(as, ac_cols)
  ofst <- fst$means$fst

  opca <- quick_smartPCA(t(genotypes[,-which(colnames(genotypes) == ".snp.id")]))
  opca <- as.data.table(opca)
  opca$pop <- facet
  omav <- manova(cbind(PC1, PC2) ~ pop, opca)
  omav <- summary(omav)$stats[1,"approx F"]


  high_fst <- which(fst$pairwise$fst >= quantile(fst$pairwise$fst, fst_cut, na.rm = TRUE))
  high_fst <- fst$pairwise$.snp.id[high_fst]

  high_fst_as <- as[which(as$.snp.id %in% high_fst),]
  fst_high <- global_fst(high_fst_as, ac_cols)
  tfst <- fst_high$means$fst


  tpca <- quick_smartPCA(t(genotypes[which(genotypes$.snp.id %in% high_fst),-which(colnames(genotypes) == ".snp.id")]))

  tpca <- as.data.table(tpca)
  tpca$pop <- facet

  mav <- manova(cbind(PC1, PC2) ~ pop, tpca)
  mav <- summary(mav)$stats[1,"approx F"]

  # cols <- c("PC1", "PC2")
  # tpca[,c("pop_means_PC1", "pop_means_PC2")  := lapply(.SD, mean), by = pop, .SDcols = cols]
  # tpca[,c("global_means_PC1", "global_means_PC2") := lapply(.SD, mean), .SDcols = cols]
  #
  # tpca[,global_dist := sqrt((PC1 - global_means_PC1)^2  + (PC2 - global_means_PC2)^2)]
  # tpca[,pop_dist := sqrt((PC1 - pop_means_PC1)^2  + (PC2 - pop_means_PC2)^2)]
  # sum_global <- sum(tpca$global_dist)
  # sum_pop <- sum(tpca$pop_dist)
  #
  # test_statistic <- sum_pop/sum_global

  if(store_pca){
    return(list(Fstat = mav, fst = tfst, delta_Fstat = mav - omav, delta_fst = tfst - ofst,
                init_Fstat = omav, init_fst = ofst,
                pca = list(all_vars = opca, selected = tpca)))
  }
  return(list(Fstat = mav, fst = tfst, delta_Fstat = mav - omav, delta_fst = tfst - ofst, init_Fstat = omav, init_fst = ofst))
}

do_boots <- function(x, genotypes, facet, nboots, par = FALSE, fst_cut = .95, store_pca){

  if(par > 2){
    cl <- parallel::makePSOCKcluster(par)
    doParallel::registerDoParallel(cl)

    out <- foreach::foreach(q = 1:nboots, .export = c("generate_summary_stats", "quick_smartPCA", "global_fst"),
                            .packages = "data.table", .inorder = FALSE) %dopar% {
                              tdata <- x$out[[q]]

                              tout <- generate_summary_stats(as = tdata, genotypes = genotypes, facet = facet[x$shuff_ord[[q]]], fst_cut = fst_cut, store_pca = store_pca)

                              if(store_pca){
                                tout <- list(res = data.table(boot = q, Fstat = tout$Fstat, fst = tout$fst,
                                                              delta_Fstat = tout$delta_Fstat, delta_fst = tout$delta_fst,
                                                              init_Fstat = tout$init_Fstat, init_fst = tout$init_fst),
                                             pca = tout$pca)
                              }
                              else{
                                tout <- data.table(boot = q, Fstat = tout$Fstat, fst = tout$fst, delta_Fstat = tout$delta_Fstat, delta_fst = tout$delta_fst, init_Fstat = tout$init_Fstat, init_fst = tout$init_fst)
                              }
                              tout
                            }

    parallel::stopCluster(cl)
    if(store_pca){
      pca <- purrr::map(out, "pca")
      out <- data.table::rbindlist(purrr::map(out, "res"))
      out <- dplyr::arrange(out, boot)

      out <- list(res = out, pca = pca)
    }
    else{
      out <- data.table::rbindlist(out)
      out <- dplyr::arrange(out, boot)
    }
  }
  else{
    out <- data.table(boot = 1:nboots)
    out$Fstat <- 0
    out$fst <- 0
    out$delta_Fstat <- 0
    out$delta_fst <- 0
    out$init_Fstat <- 0
    out$init_fst <- 0

    if(store_pca){
      pca <- vector("list", nboots)
    }

    for(i in 1:nboots){
      cat(i, "\n")
      tdata <- x$out[[i]]

      tout <- generate_summary_stats(as = tdata, genotypes = genotypes, facet = facet[x$shuff_ord[[i]]], fst_cut = fst_cut)
      if(store_pca){
        pca[[i]] <- tout$pca
      }

      out[i, Fstat := tout$Fstat]
      out[i, fst := tout$fst]
      out[i, delta_Fstat := delta_tout$Fstat]
      out[i, delta_fst := delta_tout$fst]
      out[i, init_Fstat := tout$init_Fstat]
      out[i, init_fst := tout$init_fst]
    }

    if(store_pca){
      out <- list(res = out, pca = pca)
    }
  }
  return(out)
}

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

prep_boots <- function(x, facet, n){
  x <- calc_ho(x, facet)

  as_real <- cbind(x@facet.meta, as.data.table(as.matrix(x@geno.tables$as)))
  as_real <- as_real[which(as_real$facet == facet),]
  ho_real <- get.snpR.stats(x, facet, "ho")$single
  ac_cols <- colnames(x@geno.tables$as)
  as_real <- merge(as_real, ho_real, sort = FALSE, all.x = TRUE, all.y = TRUE)

  boots <- boot_as(x, n, facet)

  return(list(real_as = as_real, boot_as = boots))
}

prep_sn <- function(x){
  sn <- format_snps(x, "sn", interpolate = FALSE, sn_remove_empty = FALSE)
  sn <- cbind(.snp.id = snp.meta(x)$.snp.id, sn[,-c(1:(ncol(snp.meta(x)) - 1))])
  return(sn)
}


test_fun <- function(x, y){
  if(all(c("cat", "dog") %in% c(x, y))){
    return("These aren't numbers.")
  }
  if(x == 3){stop("x cannot be 3, you swine.")}
  return(x + y)
}

#' Run PCA boostrapping pipeline
#'
#' Runs the full PCA boostrap/permutation testing pipeline using genetic data.
#'
#' @param x Genetic data formatted as single numbers indicated genotype, where
#'   0 and 2 are homozygoytes and 1 is a heterozygote. Rows are SNPs and columns
#'   are individuals
#'
#' @return Describe what the function returns to the user.
#'
#' @references
#' # list the citations
#'
#' @author William Hemstrom
#'
#' @author Andy Lee
#'
#' @examples
#' # provide example code here
#'
#'
#' @export
run_bootstrapping <- function(x, facet, n, fst_cut = .95, par = FALSE, store_pca = FALSE){
  x_as <- prep_boots(x, "pop", n)
  x_sn <- prep_sn(x)
  real_x <- generate_summary_stats(as = x_as$real_as, genotypes = x_sn, facet = sample.meta(x)$pop, fst_cut = fst_cut, store_pca = store_pca)
  boots_x <- do_boots(x_as$boot_as, x_sn, sample.meta(x)$pop, n, par = par, fst_cut = fst_cut, store_pca = store_pca)

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

basic_pca_one_rep <- function(x, rsp, quant, store_pca = FALSE){
  # real_pca
  dagg <- zoo::na.aggregate(d)
  pca <- as.data.frame(prcomp(x = dagg)$x)
  pca$pop <- rsp
  omav <- manova(cbind(PC1, PC2) ~ pop, pca)
  omav <- summary(omav)$stats[1,"approx F"]

  # selected pca
  select_fun <- function(x) anova(lm(as.numeric(x) ~ rsp))$`Pr(>F)`[1]
  fits <- unlist(lapply(d, select_fun))
  best_vars <- which(fits <= quantile(fits, quant))
  spca <- as.data.frame(prcomp(x = zoo::na.aggregate(d[,best_vars]))$x)
  spca$pop <- rsp
  smav <- manova(cbind(PC1, PC2) ~ pop, spca)
  smav <- summary(smav)$stats[1,"approx F"]

  if(store_pca){
    return(list(res = list(omav = omav, smav = smav, retained = best_vars),
                pca = list(all_vars = pca, selected = spca)))
  }
  else{
    return(list(omav = omav, smav = smav, retained = best_vars))
  }
}

basic_pca_bootstrapping <- function(x, rsp, quant = 0.05, nboots, store_pca = FALSE){
  real_res <- basic_pca_one_rep(x, rsp, quant, store_pca)
  # boots
  boots <- vector("list", nboots)
  for(i in 1:nboots){
    trsp <- sample(rsp, length(rsp), TRUE)
    while(length(unique(trsp)) == 1){
      trsp <- sample(rsp, length(rsp), TRUE)
    }
    boots[[i]] <- basic_pca_one_rep(x, trsp, 0.05, store_pca)
  }

  if(store_pca){
    res <- data.frame(omav = unlist(purrr::map(purrr::map(boots, "res"), "omav")),
                      smav = unlist(purrr::map(purrr::map(boots, "res"), "smav")))

    res$delta <- res$smav - res$omav

    p <- get_p_values(list(omav = real_res$res$omav, smav = real_res$res$smav,
                           delta = real_res$res$smav - real_res$res$omav), res,
                      c("greater", "greater", "greater"))
  }
  else{
    res <- data.frame(omav = unlist(purrr::map(boots, "omav")),
                      smav = unlist(purrr::map(boots, "smav")))

    res$delta <- res$smav - res$omav

    p <- get_p_values(list(omav = real_res$omav, smav = real_res$smav,
                           delta = real_res$smav - real_res$omav), res,
                      c("greater", "greater", "greater"))
  }


  if(store_pca){
    return(list(p = p, observed = list(init_Fstat = real_res$res$omav, final_Fstat = real_res$res$smav,
                                       delta_Fstat = real_res$res$smav - real_res$res$omav,
                                       pca = real_res$pca,
                                       retained = real_res$res$retained),
                null = list(values = res,
                            pca = purrr::map(boots, "pca"),
                            retained = purrr::map(purrr::map(boots, "res"), "retained"))))
  }
  else{
    return(list(p = p, observed = list(init_Fstat = real_res$omav, final_Fstat = real_res$smav,
                                       delta_Fstat = real_res$smav - real_res$omav,
                                       retained = real_res$retained),
                null = list(values = res)))
  }
}


parse_boot_res <- function(files){
  files <- lapply(files, readRDS)

  null <- purrr::map(files, "null_distribution")
  null <- rbindlist(null)

  real <- files[[1]]$observed_values

  p <- get_p_values(real, null[,-1], c("greater", "greater", "greater", "greater", "greater", "greater"))

  plot_res <- function(real, null){
    null$boot <- 1:nrow(null)
    mnull <- data.table::melt(null, id.vars = "boot")

    mobs <- as.data.table(real)
    suppressWarnings(mobs <- melt(mobs))

    p <- ggplot(mnull, aes(x = value)) +
      geom_density() +
      geom_vline(aes(xintercept = value), color = "red", data = mobs) +
      facet_wrap(~variable, ncol = 1, scales = "free") +
      theme_bw() +
      xlab("Value") +
      ylab("Density")

    return(p)

  }

  return(list(null = null, observed = real, p = p, plot = plot_res(real, null)))
}
