test_that("plot_permutation_res()", {
  set.seed(12332)

  # basic format OK
  res2 <- run_permutation(mon_sn[1:100,],
                         facet = sample(LETTERS[1:4], ncol(mon_sn), TRUE),
                         n = 10, fst_cut = 0.95, store_pca = TRUE)
  res <- res2
  res$real_pca <- NULL
  res$null_pca <- NULL

  p <- plot_permutation_res(res)

  expect_true(ggplot2::is.ggplot(p))
  expect_true(!any(is.na(p$data)))
  expect_true(length(p$layers) == 3) # since this isn't cowplot, it'll be each layer of plot: point, vline, and label.

  # basic format with PCA
  expect_true(any(names(res2) == "real_pca"))

  # check that we are returning the correct number of panels
  p2 <- plot_permutation_res(res2, 10)
  r2 <- capture.output(p2$layers[[1]])
  expect_true(length(grep("axis.title.y.left", r2)) == 22) # two per rep + real

  p3 <- plot_permutation_res(res2, 2)
  r3 <- capture.output(p3$layers[[1]])
  expect_true(length(grep("axis.title.y.left", r3)) == 6)

  # check that both have the bottom panel
  expect_true(length(p2$layers) == 2)
  expect_true(length(p2$layers) == 2)

  # plot_observed_pcas handler set to FALSE
  p4 <- plot_permutation_res(res2, 0, plot_observed_pcas = FALSE)
  expect_true(length(p4$layers) == 3) # since this isn't cowplot, it'll be each layer of plot: point, vline, and label.
  ## with perms
  p5 <- plot_permutation_res(res2, 10, plot_observed_pcas = FALSE)
  r5 <- capture.output(p5$layers[[1]])
  expect_true(length(grep("axis.title.y.left", r5)) == 20) # two per rep, no real
  expect_true(length(p5$layers) == 2) # still has bottom panel
})

