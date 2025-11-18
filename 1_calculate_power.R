library(ggplot2)
library(tidyverse)
library(plyr)

setwd("~/Desktop/aou_phenotyping/example_analysis")
source("power_functions.R")

### Set analysis-specific parameters ###########################################

# Example pathogenic variant count data for HS, LS, and control sets
df_counts <- read_csv("pathogenic_counts.csv")

K <- 0.005 # Disease prevalence

alpha <- 0.0005 # Significance level of hypothesis tests
freq_to_test <- c(0.01, 0.003, 0.001, 0.0005) # Target variant frequencies
or_min <- 1.1; or_max <- 15 # Range of target variant odds ratios

################################################################################


df_counts <- df_counts %>%
  mutate(
    freq = n_pathogenic / n
  )

# Form freq confidence intervals, if desired

score_ci <- function(pi_hat, n, alpha) {
  z_sq <- qnorm(alpha / 2, lower.tail = FALSE)^2
  
  center <- pi_hat * (n / (n + z_sq)) + 0.5 * (z_sq / (n + z_sq))
  shift <- sqrt(z_sq) * sqrt(
    1 / (n + z_sq) * (
      pi_hat * (1 - pi_hat) * (n / (n + z_sq)) +
        0.25 * (z_sq / (n + z_sq))
    )
  )
  
  return(
    c(center - shift, center + shift)
    )
}

ci_wrapper <- function(row, alpha) {
  pi_hat <- as.numeric(row["freq"])
  n <- as.numeric(row["n"])
  
  return(
    score_ci(pi_hat = pi_hat, n = n, alpha = alpha)
    )
}


ci_table <- data.frame(
  t(
    apply(
      df_counts, MARGIN = 1,
      FUN = function(x) ci_wrapper(x, alpha = 0.09)
    )
  )
)
colnames(ci_table) <- c("l_bound", "u_bound")

df_counts <- cbind( df_counts, ci_table )

f_1 <- df_counts %>% filter( case_status == "HS" ) %>% pull( freq )
f_0 <- df_counts %>% filter( case_status == "control" ) %>% pull( freq )

df_counts <- df_counts %>%
  mutate(
    ppv = 1 - (f_1 - freq) / (f_1 - f_0)
  )


# Define grid of values for power calculation

or_grid <- seq( from = or_min, to = or_max, by = 0.1 )

or <- rep( or_grid, length(freq_to_test) )
freq <- rep( freq_to_test, each = length(or_grid) )

df_grid <- data.frame(or, freq)

# Compute power

df_power <- data.frame()

sets <- df_counts$case_status
sets <- sets[sets != "control"]

n_controls <- df_counts %>% filter( case_status == "control" ) %>% pull( n )
n_hs <- df_counts %>% filter( case_status == "HS" ) %>% pull( n )

for ( s in sets ) {
  
  n_s <- df_counts %>% filter( case_status == s ) %>% pull( n )
  ppv_s <- df_counts %>% filter( case_status == s ) %>% pull( ppv )
  phi <- n_s * (1 - ppv_s) / (n_controls + n_s * (1 - ppv_s))
  
  df_s <- df_grid

  df_s$power <- apply(
    X = df_s, 
    MARGIN = 1, 
    FUN = function(x) calc_power(
      freq = as.numeric(x[2]), 
      effect_size = as.numeric(x[1]), 
      prev = K, 
      theta = 0, 
      phi = phi,
      n_cases = ifelse(
        startsWith(s, "LS"),
        n_s + n_hs,
        n_hs
      ),
      n_controls = n_controls, 
      alpha = alpha
    )
    )
  
  df_s$case_status <- s
  
  df_power <- rbind.fill( df_power, df_s )
  
}

df_power %>% write_csv("power_output.csv")