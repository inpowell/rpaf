context("test-modifiers")

test_that("categorical NAs handled correctly", {
  df <- data.frame(f = gl(2, 2, labels = letters[1:2]))
  df$f[4] <- NA
  edf <- data.frame(f = factor(c('a', 'a', 'a', NA), levels = c('a', 'b')))

  expect_equal(ModifyCategorical('f', 'a')(df), edf)
  expect_equal(ModifyCategorical('f', 'a', 'b')(df), edf)
})
