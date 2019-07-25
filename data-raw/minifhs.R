library(readr)
minifhs <- read_csv(
  "data-raw/example.csv",
  col_types = cols(
    AGEGRP = col_factor(levels = c("1", "2", "3", "4")),
    BMI_2 = col_factor(levels = c("1", "2")),
    BP = col_factor(levels = c("1", "2")),
    B_COHORT = col_factor(levels = c("1890", "1910", "1930")),
    DEATH = col_logical(),
    DIAB = col_logical(),
    SEX = col_factor(levels = c("1", "2")),
    SMOKE = col_factor(levels = c("1", "2", "3", "4"))
  ),
  na = "B"
)

levels(minifhs$AGEGRP) <- c("1" = "40-49", "2" = "50-59", "3" = "60-69",
                            "4" = "70-79")
levels(minifhs$SEX) <- c('1' = "Male", '2' = "Female")
levels(minifhs$BMI_2) <- c("1" = '<25.0', '2' = '>=25.0')
levels(minifhs$BP) <- c("1" = "Normal", "2" = "Elevated")
levels(minifhs$SMOKE) <- c('1' = 'Never', '2' = 'Former', '3' = '<30/day',
                           '4' = '>=30/day')

minifhs <- as.data.frame(minifhs)
attr(minifhs, "spec") <- NULL

usethis::use_data(minifhs, overwrite = TRUE)
