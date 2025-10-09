
# ItemRest <img src="https://img.shields.io/badge/status-active-brightgreen" align="right" />

**ItemRest** is an R package developed to automate item removal
strategies in Exploratory Factor Analysis (EFA).  
It helps researchers identify low-quality items using statistical
criteria and simulate the impact of different removal combinations on
the factor structure and internal consistency of the scale.

## 🔧 Features

- Automatically identifies cross-loading and low-loading items based on
  customizable thresholds.
- Tests all possible combinations of flagged items.
- Reports explained variance, Cronbach’s alpha, and factor loading range
  for each solution.
- Highlights optimal strategies with no cross-loading items.
- Built on top of `psych`, `GPArotation`, `EFAtools`, and `qgraph`
  packages.

## 📦 Installation

Install the package directly from GitHub using the `devtools` package:

``` r
# First install devtools if not already installed:
install.packages("devtools")

# Then install ItemRest from GitHub:
devtools::install_github("ahmetcaliskan1987/ItemRest")
```

## Example

# Load the library

library(ItemRest) \# We will use the ‘bfi’ dataset from the ‘psych’
package for a realistic example. \# This requires the ‘psych’ package to
be installed. \# The ‘bfi’ dataset contains responses to 25 personality
items.

# 1. Prepare the data: Select the personality items (first 25 columns)

# and remove rows with missing values for this example.

data(bfi, package = “psych”) example_data \<- bfi\[, 1:25\] example_data
\<- na.omit(example_data)

# 2. Run the item removal analysis.

# Based on theory, the Big Five model has 5 factors.

# Let’s run the analysis with n_factors = 5.

results \<- itemrest( data = example_data, n_factors = 5, cor_method =
“pearson” \# Data is not ordinal, so pearson is appropriate )

# 3. Print the report for optimal strategies.

# This will show the final table after the analysis is complete.

print(results, report = “optimal”)

## 📄 License

This package is distributed under the MIT License.  
See the [LICENSE](LICENSE) file for more details.
