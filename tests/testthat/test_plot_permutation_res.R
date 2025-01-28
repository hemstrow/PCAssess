?PCAssess::plot_permutation_res

test_that("plot_permutation_res()", {
  set.seed(12332)

  # basic format OK
  res <- run_permutation(mon_sn[1:100,],
                         facet = sample(LETTERS[1:4], ncol(mon_sn), TRUE),
                         n = 10, fst_cut = 0.95, store)

  p <- plot_permutation_res(res, 10)

  expect_true(ggplot2::is.ggplot(p))
  expect_true(!any(is.na(p$data)))

  p <- plot_permutation_res(res, 0)
  p$data

  # basic format with PCA
  res2 <- run_permutation(mon_sn[1:100,],
                         facet = sample(LETTERS[1:4], ncol(mon_sn), TRUE),
                         n = 10, fst_cut = 0.95, store_pca = TRUE)
  expect_true(any(names(res2) == "real_pca"))

  p2 <- plot_permutation_res(res2, 10)

})

