context("test-risk-hazard")

library("survival")

test_that("risk_survreg gives same hazard as 1/predict in exponential", {
  sr <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian,
                dist="exponential")
  risk <- risk_survreg(sr)

  expect_equal(
    risk$hazardfn(ovarian, 1),
    1 / predict(sr, newdata = ovarian, type = 'response')
  )
})

test_that("risk_survreg gives correct hazard and grad-hazard in exponential", {
  sr <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian,
                dist="exponential")
  risk <- risk_survreg(sr)

  X <- model.matrix(sr, data = ovarian)
  ehaz <- 1 / predict(sr, newdata = ovarian, type = 'response')
  edhaz <- - ehaz * X

  expect_equal(risk$hazardfn(ovarian, 1), ehaz)
  expect_equal(risk$dhazardfn(ovarian, 1), edhaz)
})
