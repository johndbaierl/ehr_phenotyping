### Description: this notebook contains the core functions used in computing
### power for association testing under estimated outcome misclassification


### Main function:

calc_power <- function(
    freq, effect_size, prev, theta, phi,
    n_cases, n_controls, alpha
) {
  #############################################################################################
  # Description: computes power of association tests under phenotype misclassification based on
  #     Edwards et al 2005
  #
  # Parameters:
  #    freq: carrier frequency
  #    effect_size: variant effect size
  #    prev: prevalence
  #    theta: Pr(affected misclasssified as control) = 1 - Sensitivity
  #    phi: Pr(unaffected misclassified as a case) = 1 - Specificity
  #    n_cases: number of cases
  #    n_controls: number of controls
  #    alpha: significance level
  #############################################################################################
  
  roots <- quad_roots(
    a = (1 - prev) * (1 - effect_size),
    b = (effect_size * (freq - prev) - 1 + prev - freq),
    c = freq
  )
  p1j <- roots[roots >= 0]
  
  if (identical(p1j, numeric(0))) {
    return("ERROR: frequency computation returned no positive roots")
  }
  
  if (length(p1j) > 1) {
    return("ERROR: frequency computation returned multiple positive roots")
  }
  
  p0j <- (freq - (1 - prev) * p1j) / prev
  
  genotype <- c("hom wildtype", "het")
  p0j <- c(1 - p0j, p0j)
  p1j <- c(1 - p1j, p1j)
  df_geno_sum <- data.frame(genotype, p0j, p1j)
  
  p0j_star_vec <- apply(
    X = df_geno_sum, MARGIN = 1, 
    function(x) p0j_star(
      as.numeric(x[2]), 
      as.numeric(x[3]), 
      theta, phi, K
    )
  )
  
  p1j_star_vec <- apply(
    X = df_geno_sum, MARGIN = 1, 
    function(x) p1j_star(
      as.numeric(x[2]), 
      as.numeric(x[3]), 
      theta, phi, K
    )
  )
  
  df_geno_sum$p0j_star <- p0j_star_vec
  df_geno_sum$p1j_star <- p1j_star_vec
  
  Na_star <- n_cases
  Nu_star <- n_controls
  ncp <- 0
  
  for (i in seq(nrow(df_geno_sum))) {
    p0j_star <- df_geno_sum$p0j_star[i]
    p1j_star <- df_geno_sum$p1j_star[i]
    
    next_term <- Na_star * Nu_star * (
      (p0j_star - p1j_star)^2 / (Na_star * p0j_star + Nu_star * p1j_star)
    )
    
    ncp <- ncp + next_term
  }
  
  power <- pchisq(
    qchisq(alpha, df = 1, lower = F), 
    df = 1,
    ncp = ncp, 
    lower = F
  )
  
  return(power)
}


### Helper functions:

quad_roots <- function(
    a, b, c
    ) {
  center <- -b / (2 * a)
  shift <- sqrt(b^2 - 4 * a * c) / (2 * a)
  
  return(
    sort(c(center - shift, center + shift))
    )
}


# For adjusting frequencies in misclassification calculation

p0j_star <- function(
    p0j, p1j, theta, phi, K
    ) {
  return(
    (p0j * (1 - theta) * K + p1j * phi * (1 - K)) /
      ((1 - theta) * K + phi * (1 - K))
  )
}

p1j_star <- function(
    p0j, p1j, theta, phi, K
    ) {
  return(
    (p0j * theta * K + p1j * (1 - phi) * (1 - K)) /
      (theta * K + (1 - phi) * (1 - K))
  )
}

### Other

compute_ppv <- function(
    df,
    confident_group,
    tentative_groups,
    control_group = 0,
    max_one = TRUE,
    overall_ppv = FALSE # indicates whether PPV should include confident cases or just tentative cases
) {
  ### Description: takes pooled df of counts and frequencies from phenotyping subgroups and 
  ###     computes PPV among tentative case groups
  
  # Preliminaries:
  
  require(tidyverse)
  
  cols_to_require <- c("group_label", "freq")
  
  if ( !all( cols_to_require %in% names(df) ) ) {
    stop("expected columns in df not found")
  }
  
  n_groups <- length( tentative_groups )
  
  df_out <- data.frame(
    group_label = tentative_groups,
    ppv = numeric( length(tentative_groups) )
  )
  
  # Compute PPV
  
  f_1 <- df %>%
    filter( group_label == confident_group ) %>%
    pull( freq )
  
  f_0 <- df %>%
    filter( group_label == control_group ) %>%
    pull( freq )
  
  f <- df %>%
    filter( group_label %in% tentative_groups ) %>%
    pull( freq )
  
  df_out$ppv <- 1 - (f_1 - f) / (f_1 - f_0)
  
  if ( max_one ) {
    df_out$ppv[df_out$ppv > 1] <- 1
  }
  
  if ( overall_ppv ) {
    
    n_tent <- df %>%
      filter(
        group_label %in% tentative_groups # retains same group_label ordring as df_out$ppv
      ) %>%
      pull( n )
    
    n_conf <- df %>%
      filter(
        group_label == confident_group
      ) %>%
      pull( n )
    
    df_out$ppv <- (df_out$ppv * n_tent + n_conf) / (n_tent + n_conf)
    
  }
  
  return( df_out )
}

score_ci <- function(pi_hat, n, alpha) {
  
  ### Description: this produces score-based confidence intervals for sample proportions
  ###     at specified confidence level
  
  z_sq <- qnorm(alpha / 2, lower.tail = FALSE)^2
  
  center <- pi_hat * (n / (n + z_sq)) + 0.5 * (z_sq / (n + z_sq))
  shift <- sqrt(z_sq) * sqrt(
    1 / (n + z_sq) * (
      pi_hat * (1 - pi_hat) * (n / (n + z_sq)) +
        0.25 * (z_sq / (n + z_sq))
    )
  )
  
  return(c(center - shift, center + shift))
}

ci_wrapper <- function(row, alpha) {
  pi_hat <- as.numeric(row["freq"])
  n <- as.numeric(row["n"])
  
  return(score_ci(pi_hat = pi_hat, n = n, alpha = alpha))
}
