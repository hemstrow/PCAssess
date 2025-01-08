library(testthat)
library(PCAssess)

?PCAssess::run_permutation_res

test_that("Check input files", (
  expect_error(run_permutation(x, facet)))
)

test_that("check fst cut-offs", {
  expect_error(run_permutation(fst_cut < 0))
  expect_error(run_permutation(fst_cut > 1))
})

# expect_output(x, y)
# expect_message(x, y)
# expect_warning(x, y)
# expect_error(x, y)
