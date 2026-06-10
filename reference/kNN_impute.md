# Fast and Precise k-Nearest Neighbours Imputation for Mixed Data

Imputes missing values (NAs) in a given dataset using kNN algorithm. The
function supports Euclidean and Manhattan distances, designed for fast
numeric imputation and a Gower distance metric for mixed data. The
function operates in two modes: "Fast", which is meant to be used for
rapid, memory-efficient imputation and "Precise mode", which maximizes
the accuracy of the imputation.

## Usage

``` r
kNN_impute(
  data,
  k = 5,
  metric = c("gower", "euclidean", "manhattan"),
  mode = c("precise", "fast"),
  num_fun = c("median", "mean"),
  threads = NULL,
  maxColNa = 0.8,
  maxRowNa = 1
)
```

## Arguments

- data:

  A data.frame or matrix, which contains a missing values (NAs) to be
  imputed.

- k:

  An integer specifying the number of neighbours. Default: 5.

- metric:

  A string specifying the distance metric to use. Allowed values are
  `"gower"` (default), `"euclidean"`, or `"manhattan"`.

- mode:

  A string specifying the search algorithm. `"fast"` calculates global
  neighbors `"precise"` (default) recalculates neighbors for each
  missing cell. Uses weights in Gower distance metric.

- num_fun:

  A string specifying the function used to impute missing cells. Allowed
  values are `"median"` (default) or `"mean"`.

- threads:

  An integer specifying the number of threads used by the function.
  Defaults to all available CPU cores minus one.

- maxColNa:

  A numeric value between 0.0 and 1.0. Columns with fraction of missing
  values greater than this number will be removed. Default value is 0.8.

- maxRowNa:

  A numeric value between 0.0 and 1.0. Rows with fraction of missing
  values greater than this number will be removed. Default value is 1.0.

## Value

A data.frame with missing values replaced by imputed estimates.

## Examples

``` r
if (FALSE) { # \dontrun{
 df <- data.frame(
   PatientId = as.integer(1:8),
   HeartRate = c(72, 85, 70, 88, NA, 60, 78, 86),
   Cholesterol = c(190, 240, 185, 250, 200, NA, 195, 255),
   RiskLevel = as.factor(c("Low", "High", "Low", "High", "Moderate", "Low", "Low", NA))
 )
 
 #Fast imputation using Euclidean distance
 kNNFast <- kNN_impute(df, k = 2, metric = "euclidean", mode = "fast")
 # kNNFast$Cholesterol[6]
 # [1] 192.5
 
 #Highly accurate imputation using Gower distance
 kNNPrecise <- kNN_impute(df, k = 2, metric = "gower", mode = "precise")
 # kNNPrecise$RiskLevel[8]
 # [1] High

} # }
```
