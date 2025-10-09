## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ItemRest)

## ----data_prep----------------------------------------------------------------
# Ensure the 'psych' package is available
if (requireNamespace("psych", quietly = TRUE)) {
  data(bfi, package = "psych")
  
  # Select the personality items (first 25 columns)
  analysis_data <- bfi[, 1:25]
  
  # Omit rows with missing values for this example
  analysis_data <- na.omit(analysis_data)
  
  # View the first few rows of the prepared data
  head(analysis_data)
}

## ----run_analysis-------------------------------------------------------------
# Run the analysis
if (exists("analysis_data")) {
  results <- itemrest(
    data = analysis_data,
    n_factors = 5,
    cor_method = "pearson"
  )
}

## ----print_optimal------------------------------------------------------------
if (exists("results")) {
  # Print the default optimal report
  print(results, report = "optimal")
}

## ----print_all----------------------------------------------------------------
if (exists("results")) {
  # Print the report for all tested strategies
  print(results, report = "all")
}

