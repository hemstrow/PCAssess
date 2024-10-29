test_that("this is an example test",{
  expect_true(all(c(4, 5) < test_fun(4, 5)))
  expect_true(4 + 5 == test_fun(4, 5))
  expect_error(test_fun(3, 5), "swine")

  expect_true(test_fun("cat", "dog") == "These aren't numbers.")
  expect_true(test_fun("dog", "cat") == "These aren't numbers.")
  expect_true(test_fun(c("dog", 2), "cat") == "These aren't numbers.")
})
