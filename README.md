# CondIntCal (replication material)

This repository contains `R` code to reproduce the results in the pre-print  

> S. Allen, J. Burnello, and J. Ziegel (2025). 
> Assessing the conditional calibration of interval forecasts using decompositions of the interval score.
> ArXiv preprint [arXiv:2508.18034](https://arxiv.org/abs/2508.18034).

The methods developed in this paper are available in the more general `CondIntCal` package at https://github.com/sallen12/CondIntCal.


## Conditional calibration of interval forecasts

The above paper concerns calibration of forecasts in the form of central prediction intervals. 

The calibration of interval forecasts is typically assessed by checking whether the outcome falls within the interval with the nominal coverage probability. However, this is typically assessed unconditionally, therefore neglecting conditional biases that may be present in the forecast. This work introduces a conditional decomposition of the interval score that provides a measure of conditional calibration of interval forecasts. 

This repository provides the functionality to implement these interval score decompositions in practice. Functions are additionally available to assess the unconditional coverage and average length of the forecasts, and to display miscalibration-discrimination plots that highlight trade-offs between the conditional calibration and information content in the interval forecasts.

## Data

The data used in this study is a subset of the data available from https://github.com/yromano/cqr. The interval forecasts generated from the different prediction methods in the case study can be found in the `Case_Study/data` folder.

## Code

The simulation study results can be obtained by sourcing the `main.R` file in the `Sim_Study` folder. The case study results can be obtained by sourcing the `main.R` file in the `Case_Study` folder. The variable `var` in this latter file can be changed to get results for each of the three datasets considered. Both of these files use generic functions to evaluate interval forecasts available in the `utility_funcs.R` file.
