.boot_as <- function(x, n, facet){
  ..tm <- ..ord <-  NULL

  out <- vector("list", n)

  shuff_ord <- out

  # we could do all of this vectorized... but it would get very memory intense. Loop instead.
  for(i in 1:n){
    fshuf <- sample(facet, length(facet), FALSE)
    out[[i]] <- prep_as_from_sn(x, fshuf)
    out[[i]]$amat[,ho := out[[i]]$ho]
    out[[i]]$facet <- fshuf
  }

  return(out)
}

.fix..call <- function(fun){
  return(.suppress_specific_warning(fun, "variable in calling scope for clarity"))
}

.per_all_f_stat_components <- function(intot = NULL, jntot = NULL, ps1 = NULL, ps2 = NULL, r, nbar, nc, iho = NULL, jho = NULL, ntotm = NULL, psm = NULL, hom = NULL){
  if(r == 2 & !is.null(intot)){
    pbar <- ((intot*ps1) + (jntot*ps2))/(r*nbar) #average sample allele frequency
    ssq <- (((intot)*(ps1-pbar)^2) + ((jntot)*(ps2-pbar)^2))/((r-1)*nbar) #sample variance of allele frequencies
    hbar <- ((intot*iho) + (jntot*jho))/(r*nbar) #average heterozygote frequencies
  }
  else if (r == 1){
    pbar <- ps1 #average sample allele frequency
    ssq <- 0 #sample variance of allele frequencies
    hbar <- iho #average heterozygote frequencies
  }
  else if(r > 2 | (r == 2 & is.null(intot))){
    pbar <- rowSums(ntotm * psm)/(r*nbar) #average sample allele frequency
    ssq <- rowSums(ntotm*(psm-pbar)^2)/((r - 1)*nbar)
    hbar <- rowSums(ntotm*hom)/(r*nbar)
  }



  #equation parts used in both
  inner1 <- pbar*(1-pbar)
  inner2 <- ((r-1)/r)*ssq
  inner3 <- .25*hbar

  inner4 <- ((2*nbar - 1)/(4*nbar))*hbar
  a <- (nbar/nc) * (ssq - (1/(nbar - 1))*(inner1 - inner2 - inner3))
  b <- (nbar/(nbar-1))*(inner1 - inner2 - inner4)
  c <- .5*hbar

  # # browser()
  # if(.print & !is.null(intot)){
  #   saveRDS(data.frame(pbar = pbar, ssq = ssq, hbar = hbar, iho = iho, intot = intot, jntot = jntot, jho = jho, inner1 = inner1, inner2 = inner2, inner3 = inner3, inner4 = inner4, a = a, b = b, c = c), "test.RDS")
  # }
  # else if(.print){browser()}
  return(list(a = a, b = b, c = c))

  # weir--exactly the same
  # else{
  #   S1 <- ssq - (1/(nbar-1))*(inner1 - inner2 - inner3)
  #   S2i1 <- ((r*(nbar - nc))/nbar)*inner1
  #   S2i2 <- (1/nbar)*((nbar-1)+(r-1)*(nbar-nc))*ssq
  #   S2i3 <- ((nbar-nc)/(4*nc^2))*hbar
  #   S2 <- inner1 - (nbar/(r*(nbar-1)))*(S2i1 -S2i2 - S2i3)
  #   return(list(S1 = S1, S2 = S2))
  # }


}




.do_boots <- function(x, genotypes, par = FALSE, fst_cut = .95, store_pca){

  nboots <- length(x)
  if(par > 2){
    cl <- parallel::makePSOCKcluster(par)
    doParallel::registerDoParallel(cl)
    out <- foreach::foreach(q = 1:nboots, .export = c("generate_summary_stats", "global_fst", "quick_smartPCA"),
                            .packages = c("data.table"), .inorder = FALSE) %dopar% {
                              tdata <- x[[q]]

                              tout <- generate_summary_stats(as = tdata$amat, genotypes = genotypes, facet = tdata$facet, fst_cut = fst_cut, store_pca = store_pca)

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


      tout <- generate_summary_stats(as = x[[i]]$amat, genotypes = genotypes, facet = x[[i]]$facet, fst_cut = fst_cut)
      if(store_pca){
        pca[[i]] <- tout$pca
      }

      out[i, Fstat := tout$Fstat]
      out[i, fst := tout$fst]
      out[i, delta_Fstat := tout$delta_Fstat]
      out[i, delta_fst := tout$delta_fst]
      out[i, init_Fstat := tout$init_Fstat]
      out[i, init_fst := tout$init_fst]
    }

    if(store_pca){
      out <- list(res = out, pca = pca)
    }
  }
  return(out)
}

.prep_boots <- function(x, facet, n){
  x <- as.data.table(x)
  as <- prep_as_from_sn(x, facet)

  as_real <- as$amat
  as_real[,ho := as$ho]

  boots <- .boot_as(x, n, facet)

  return(list(real_as = as_real, boot_as = boots))
}

.suppress_specific_warning <- function(fun, warn_to_suppress){
  warn_handle <- function(warn){
    if(any(grepl(warn_to_suppress, warn))){
      invokeRestart("muffleWarning")
    }
  }

  withCallingHandlers(res <- fun,
                      warning = warn_handle)


  return(res)
}

