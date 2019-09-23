context("test-paf-missing")

test_that("missing data handling in paf functions", {
  df <- data.frame(
    resp = gl(2, 2, 4, labels = c('-', '+')),
    pred = gl(2, 1, 4, labels = c('-', '+')),
    count = c(65, 20, 5, 10)
  )
  model <- glm(resp ~ pred, family = binomial, data = df, weights = count)

  risk <- risk_glm(model)
  pop_df <- data.frame(
    pred = c('-', '+')
  )[rep.int(1:2, c(30, 70)), , drop = FALSE]
  pop_df_na <- data.frame(
    pred = c('-', '+', NA)
  )[rep.int(1:3, c(30, 70, 10)), , drop = FALSE]

  expect_equal(
    paf(risk, pop_df_na, mod_df = ModifyCategorical('pred', '-')),
    paf(risk, pop_df, mod_df = ModifyCategorical('pred', '-'))
  )
  expect_equal(
    paf_ci(risk, pop_df_na, mod_df = ModifyCategorical('pred', '-')),
    paf_ci(risk, pop_df, mod_df = ModifyCategorical('pred', '-'))
  )
})
