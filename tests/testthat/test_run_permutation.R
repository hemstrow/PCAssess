test_that("run_permutation", {

  set.seed(12332)
  # basic format OK
  res <- run_permutation(mon_sn[1:100,],
                         facet = sample(LETTERS[1:4], ncol(mon_sn), TRUE),
                         n = 10, fst_cut = 0.95)

  expect_identical(names(res), c("null_distribution", "observed_values", "pvalues"))

  ## pvals
  expect_true(is.numeric(res$pvalues))
  expect_identical(names(res$pvalues), c("Fstat", "fst", "delta_Fstat", "delta_fst", "init_Fstat", "init_fst"))
  expect_true(all(res$pvalues >= 0))
  expect_true(all(res$pvalues <= 1))

  ## null
  expect_true(is.data.table(res$null_distribution))
  skip_if_not(is.data.table(res$null_distribution))
  expect_identical(colnames(res$null_distribution), c("boot", "Fstat", "fst", "delta_Fstat", "delta_fst", "init_Fstat", "init_fst"))
  expect_true(all(unlist(lapply(res$null_distribution, class)) %in% c("integer", "numeric")))
  expect_true(all(is.finite(unlist(res$null_distribution)))) # check that everything is a actual number

  ## obs
  expect_true(is.numeric(unlist(res$observed_values)))
  expect_identical(names(res$observed_values), c("Fstat", "fst", "delta_Fstat", "delta_fst", "init_Fstat", "init_fst"))
  expect_true(all(is.finite(unlist(res$observed_values)))) # check that everything is a actual number

  # does this actually do what we expect.
  skip_on_cran()
  set.seed(12332)
  res <- run_permutation(mon_sn,
                         facet = sample(LETTERS[1:4], ncol(mon_sn), TRUE),
                         n = 10, fst_cut = 0.95, par = 4)
  expect_true(res$observed_values$Fstat > res$observed_values$init_Fstat) # there should be more clustering!
  expect_true(all(res$null_distribution$delta_Fstat > 0)) # in all boots!
  expect_false(res$pvalues["delta_Fstat"] <= 0.1) # but it shouldn't be a significant increase.
})

test_that("check errors and options", {
  # test fst cut-off options
  expect_error(run_permutation(mon_sn[1:100,],
                               facet = sample(LETTERS[1:4], ncol(mon_sn), TRUE),
                               n = 10, fst_cut = 1.05),
               regexp = "Fst cut-off must be between 0 and 1")
  # TODO: add reporting for the number of SNPs maintained (verbose argument to generate_summary_stats?), test that it is correct
  # TODO: add a sanity check to run_permutation that checks that the fst_cut will actually retain some SNPs (cut = .9999 with 10 SNPs will eliminate all?)
  # run_permutation(mon_sn[1:100,],
  #                 facet = sample(LETTERS[1:4], ncol(mon_sn), TRUE),
  #                 n = 10, fst_cut = .9)
})


