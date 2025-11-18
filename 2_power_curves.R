library(ggplot2)
library(tidyverse)

### Description: this notebook produces power curves equivalent to manuscript
###   Figures 3,4 from the output of `1_calculate_power.R`

setwd("~/aou_phenotyping/example_analysis/")

df_power <- read_csv("power_output.csv")

# Set names for case_status sets:

names_key <- data.frame(
  case_status = unique(df_power$case_status),
  case_criteria = c("eMERGE", "1+", "2+")
)

df_power <- df_power %>% left_join(
  names_key,
  by = "case_status"
)

df_power <- df_power %>% 
  mutate(
    freq = factor(
      case_when(
        freq == 0.01 ~ "Carrier freq. = 0.01",
        freq == 0.003 ~ "Carrier freq. = 0.003",
        freq == 0.001 ~ "Carrier freq. = 0.001",
        freq == 5e-4 ~ "Carrier freq. = 0.0005"
      )
    )
  ) %>%
  mutate(
    freq = factor(
      freq, 
      levels = c(
        "Carrier freq. = 0.01",
        "Carrier freq. = 0.003",
        "Carrier freq. = 0.001",
        "Carrier freq. = 0.0005"
        )
      )
    )

ggsave("power.png", plot = p_power, dpi = 300)