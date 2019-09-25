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

  expect_equal(risk$hazardfn(ovarian, 2), ehaz)
  expect_equal(risk$dhazardfn(ovarian, 2), edhaz)
})

test_that("risk_pch gives correct hazard", {
  period_base <- factor(1:3, levels = 1:3,
                        labels = c("(0,60]", "(60,120]", "(120,1000]"))
  vet2 <- survSplit(Surv(time, status) ~ ., veteran,
                    cut = c(60, 120), episode = ".period", id = ".ID")
  vet2$.period <- period_base[vet2$.period]

  sr <- survreg(Surv(time - tstart, status) ~ karno * .period + age + trt,
                data = vet2, dist = "exponential")

  newvet <- merge(cbind(veteran, ".ID" = seq_len(nrow(veteran))),
                  data.frame(".period" = period_base))

  X <- model.matrix(sr, data = newvet)
  ehaz <- 1 / predict(sr, newdata = newvet, type = 'response')
  edhaz <- - ehaz * X

  ehaz_array <- vec2mat(ehaz, newvet$.ID, newvet$.period)
  edhaz_array <- mat2arr(edhaz, newvet$.ID, newvet$.period)

  risk <- risk_pch(Surv(time, status) ~ karno * .period + age + trt,
                   data = veteran, breaks = c(0, 60, 120, 1000))

  expect_equal(risk$hazardfn(veteran, 1000), ehaz_array)
  expect_equal(risk$dhazardfn(veteran, 1000), edhaz_array)
})
