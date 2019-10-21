# See Laaksonen et al. (2011) for these risk estimates
RR <- list(
  SEX = c("Male" = 1, "Female" = 1.0302),
  SMOKE = c("Never" = 1, "Former" = 1.4179, "<30/day" = 1.1446, ">=30/day" = 2.7868)
)

ci <- list(
  SEX = rbind(
    "lwr" = c("Male" = 1, "Female" = 0.7504),
    "upr" = c("Male" = 1, "Female" = 1.4144)
  ),
  SMOKE = rbind(
    "lwr" = c("Never" = 1, "Former" = 0.9863, "<30/day" = 0.7724, ">=30/day" = 1.2521),
    "upr" = c("Never" = 1, "Former" = 2.0385, "<30/day" = 1.6962, ">=30/day" = 6.2023)
  )
)

se <- lapply(ci, ci2se, trans = log)
se

risk <- risk_relative(RR, se, eta = 1)

paf_ci(risk, minifhs, ModifyCategorical("SMOKE", "Never", "<30/day", ">=30/day"), na.action = na.pass)
