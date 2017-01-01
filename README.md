[![Coverage Status](https://coveralls.io/repos/pkofod/QuantileRegression.jl/badge.svg?branch=master)](https://coveralls.io/r/pkofod/QuantileRegression.jl?branch=master)
[![Build Status](https://travis-ci.org/pkofod/QuantileRegression.jl.svg?branch=master)](https://travis-ci.org/pkofod/QuantileRegression.jl)

# Quantile regression in the Julia language

A very simple (and mostly untested) implementation of quantile regression.

* Install using `Pkg.clone(https://github.com/pkofod/QuantileRegression.jl)`
* Main author: Patrick Kofod Mogensen
* Contact: Use the [https://github.com/pkofod/QuantileRegression.jl/issues](issues) page
* License: BSD-3

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

# History
This package was originally created as a port of the reweighed least squares code
from the python project [https://github.com/statsmodels/statsmodels](statsmodels)
by Vincent Arel-Bundock. All contributions can be seen via the [https://github.com/pkofod/QuantileRegression.jl/graphs/contributors](contributors) page.
