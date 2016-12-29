[![Coverage Status](https://coveralls.io/repos/vincentarelbundock/QuantileRegression.jl/badge.svg?branch=master)](https://coveralls.io/r/vincentarelbundock/QuantileRegression.jl?branch=master)
[![Build Status](https://travis-ci.org/vincentarelbundock/QuantileRegression.jl.svg?branch=master)](https://travis-ci.org/vincentarelbundock/QuantileRegression.jl)

# This package is no longer supported.

I will not accept PRs or act on requests for modifications. I'm leaving the repository up in case someone wants to fork and maintain.

# Quantile regression in the Julia language

A very simple (and mostly untested) implementation of quantile regression.

* https://github.com/vincentarelbundock/QuantileRegression.jl
* Author: Vincent Arel-Bundock (with changes by John Myles White and Patrick Kofod Mogensen)
* Contact: varel@umich.edu
* License: BSD-3
* Original code from the python statsmodels project
    - https://github.com/statsmodels/statsmodels

# Example

The file ``examples/qreg_example.jl`` shows how to use the functions provided here. It replicates part of the analysis in:

* Koenker, Roger and Kevin F. Hallock. "Quantile Regression". Journal of Economic Perspectives, Volume 15, Number 4, Fall 2001, Pages 143–156

We are interested in the relationship between income and expenditures on food for a sample of working class Belgian households in 1857 (the Engel data), so we estimate a least absolute deviation model.

    using QuantileRegression

    # Load data
    url = "http://vincentarelbundock.github.io/Rdatasets/csv/quantreg/engel.csv"
    Data = readtable("engel.csv")

    # Fit least absolute deviation model (quantile  = .5)
    > ResultQR = qreg(foodexp~income, Data, .5)
    > β = coeftable(ResultQR)
    > show(β)

    2x5 DataFrame:
                 Estimate   Std.Error   t value
    (Intercept)  81.4823    14.6345      5.56783
    income        0.560181   0.0131756  42.5164


The results look pretty close to Stata 12's ``qreg``:

    . insheet using engel.csv
    . qreg foodexp income, vce(iid, kernel(epan2))
    Median regression                                    Number of obs =       235
      Raw sum of deviations 46278.06 (about 582.54126)
      Min sum of deviations 17559.93                     Pseudo R2     =    0.6206

    ------------------------------------------------------------------------------
         foodexp |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
          income |   .5601805   .0131763    42.51   0.000     .5342206    .5861403
           _cons |   81.48233   14.63518     5.57   0.000     52.64815    110.3165
    ------------------------------------------------------------------------------

We can also compute and plot (using Julia's Winston) results for various quantiles. Full code to produce the figure is in the examples folder.

![](./examples/qreg_example_plot.png)
