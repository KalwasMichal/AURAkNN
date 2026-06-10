# AURAkNN: Accelerated Ultra-rapid R Architecture for k-NN imputation

**AURAkNN** is a high-performance R package designed for fast, precise,
and memory-efficient k-Nearest Neighbors (kNN) imputation of missing
values in mixed-type datasets.

## Installation

### Windows

In order to install this package you need to have Rtools installed.
Download int from:
[CRAN](https://cran.r-project.org/bin/windows/Rtools/). Then run in R:

``` r

install.packages("remotes")
remotes::install_github("KalwasMichal/AURAkNN")
library(AURAkNN)
```

### Linux (Ubuntu/Debian)

In case r-base-dev is not installed on your system type in terminal:

``` bash
sudo apt -get install r-base-dev
```

Then run in R:

``` r

install.packages("remotes")
remotes::install_github("KalwasMichal/AURAkNN")
library(AURAkNN)
```

### macOS

macOS required Xcode Command Line Tools for C compilation, make they are
installed on your system:

``` bash
xcode-select --install
```

Then run in R:

``` r

install.packages("remotes")
remotes::install_github("KalwasMichal/AURAkNN")
library(AURAkNN)
```

## Features

- **High-Performance C Engine**: Distance and imputation computations
  are implemented in C, maximizing efficiency.
- **OpenMP Parallelization**: Neighbors search and distance calculations
  are parallelized using OpenMP.
- **Numeric and Categorical Data**: Handling both numeric and
  categorical columns, with one-hot encoding and z-score normalization
  for Eclidean/Manhattan metric and native handling for Gower.
- **Flexible Distance Metrics**: Supports `"euclidean"` and
  `"manhattan"` metrics, intended for purely numeric data, and `"gower"`
  for mixed data.
- **Two Imputation Modes**:
  - `precise`: Finds neighbors for each missing cell using only known
    values. Combined with Gower distance, uses a correlation-based
    weights ((Spearman, Cramér’s V, R²)
  - `fast`: Calculates global neighbors for each row, allowing faster
    and memory-effiecient imputation.
- **Pre-filtering And Fallbacks**: Columns/rows exceeding configurable
  NA treshold are dropped before imputation. A fallback ensures no
  missing values left in the imputeted data set.

## Examples

### Fast imputation using Euclidean distance

``` r

library(AURAkNN)

df <- data.frame(
  PatientId   = as.integer(1:8),
  HeartRate   = c(72, 85, 70, 88, NA, 60, 78, 86),
  Cholesterol = c(190, 240, 185, 250, 200, NA, 195, 255),
  RiskLevel   = as.factor(c("Low", "High", "Low", "High", "Moderate", "Low", "Low", NA))
)

kNNFast <- kNN_impute(df, k = 2, metric = "euclidean", mode = "fast")
kNNFast$Cholesterol[6]
# [1] 192.5
```

### Highly accurate imputation using Gower distance

``` r

kNNPrecise <- kNN_impute(df, k = 2, metric = "gower", mode = "precise")
kNNPrecise$RiskLevel[8]
# [1] High
```

## Results

Performance was evaluated on the
[`ggplot2::diamonds`](https://ggplot2.tidyverse.org/reference/diamonds.html)
dataset — a random sample of **10 000 diamonds**. All 10 columns were
used, with **20% missing values** injected randomly into each column.
All methods were run with `k = 10`,Gower distance, 4 threads and
seed=2115.

| Column    | Type    | Description                                 |
|-----------|---------|---------------------------------------------|
| `carat`   | numeric | Diamond weight                              |
| `depth`   | numeric | Total depth percentage                      |
| `table`   | numeric | Width of top facet relative to widest point |
| `price`   | numeric | Price in USD                                |
| `x`       | numeric | Length (mm)                                 |
| `y`       | numeric | Width (mm)                                  |
| `z`       | numeric | Depth (mm)                                  |
| `cut`     | factor  | Cut quality (5 levels)                      |
| `color`   | factor  | Diamond color (7 levels)                    |
| `clarity` | factor  | Clarity grade (8 levels)                    |

### Numeric accuracy — RMSE

| Column  | AURAkNN precise/gower | AURAkNN fast/gower | VIM::kNN  |
|---------|-----------------------|--------------------|-----------|
| `carat` | 0.0531                | 0.1762             | 0.1694    |
| `depth` | 1.2647                | 1.3409             | 1.3509    |
| `table` | 1.8100                | 1.9733             | 1.9929    |
| `price` | 1261.7722             | 2106.1396          | 1350.4876 |
| `x`     | 0.1020                | 0.4229             | 0.2196    |
| `y`     | 0.1030                | 0.4157             | 0.2060    |
| `z`     | 0.1476                | 0.2556             | 0.1608    |

### Categorical accuracy

| Column    | AURAkNN precise/gower | AURAkNN fast/gower | VIM::kNN |
|-----------|-----------------------|--------------------|----------|
| `cut`     | 0.5755                | 0.4940             | 0.4405   |
| `color`   | 0.2385                | 0.2235             | 0.2380   |
| `clarity` | 0.2720                | 0.2740             | 0.2760   |

### Execution time

| Package               | Elapsed (s) |
|-----------------------|-------------|
| AURAkNN precise/gower | 2.667       |
| AURAkNN fast/gower    | 1.463       |
| VIM::kNN              | 248.168     |
