context("test-utils")

test_that("vec2mat works as expected", {
  dat <- c(1, 2, 3, 4, 5, 6, 7, 8)
  id <-  c(4, 2, 3, 1, 2, 3, 1, 4)
  per <- c(2, 2, 2, 2, 1, 1, 1, 1)

  expect_equal(
    vec2mat(dat, id, per),
    matrix(c(7, 5, 6, 8, 4, 2, 3, 1), nrow = 4, dimnames = list(1:4, 1:2))
  )
})

test_that("mat2arr works as expected", {
  mat <- matrix(1:12, nrow = 6, ncol = 2, dimnames = list(NULL, c('a', 'b')))
  id <- c(1, 2, 3, 3, 2, 1)
  per <- c(1, 2, 1, 2, 1, 2)

  exp <- array(
    c(1,  5, 3,    6, 2,  4,
      7, 11, 9,   12, 8, 10),
    dim = c(3, 2, 2),
    dimnames = list(1:3, 1:2, c('a', 'b'))
  )

  expect_equal(mat2arr(mat, id, per), exp)
})
