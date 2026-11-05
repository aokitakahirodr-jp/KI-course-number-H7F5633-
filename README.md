# KI course H7F5633 – Assignments

This repository contains my solutions and documentation for the course
assignments of KI course H7F5633.

**Author:** Takahiro Aoki

## Repository structure

```
assignment1/
├── assignment1.md      # Documentation of Assignment 1 (Tasks 1–6)
├── Task7_8.R           # R code for Tasks 7–8
├── assignment1.html    # Knitted report
├── data/               # Input data (see below)
└── figures/            # Saved plots (png, pdf, svg)
```

## How to reproduce

1. Open `Task7_8.R` in RStudio.
2. Install the required packages:

```r
install.packages(c("tidyverse", "svglite", "devtools"))
devtools::install_github("hirscheylab/tidybiology")
```

3. Knit the document.

## Data

- `magic_guys.csv` and `microarray_data.tab` were provided in the shared
  course folder (Course week 1).
- `CO2` is a built-in R dataset.
- `chromosome` and `proteins` are from the tidybiology package.

## Software
