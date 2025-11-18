# ehr_phenotyping
This repository contains the scripts for performing power analysis for genetic association testing that uses  from the manuscript "Measuring the accuracy of EHR-based phenotyping in the *All of Us* Research Program to optimize statistical power for genetic association testing".

* `pathogenic_counts.csv` contains example pathogenic variant carrier counts and sample sizes for a high stringency (HS) case set, two low stringency (LS) case sets, and a control set.  These can be adjusted to either replicate the analysis from the original manuscript or to apply this procedure to additional projects.
*  `power_functions.R` contains the main functions used calculating power under outcome misclassification.
*  `1_calculate_power.R` outputs predicted power with the specified case sets over a grid of target variant effect sizes and carrier frequencies.  Other parameters such as disease prevalence and significance level of association tests can be adjusted here.
*  `2_power_curves.R` produces plots of statistical power and can be used to reproduce Figures 3 and 4 from the manuscript.
