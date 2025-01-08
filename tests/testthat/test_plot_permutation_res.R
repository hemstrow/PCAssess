?PCAssess::plot_permutation_res

test_that("Check input files", {
  expect_error(run_permutation(x, facet))
})

test_that("n_boots", {
  expect_true(is.numeric(n_boots))
})
