# QReg

A very simple (and untested) implementation of quantile regression.

* https://github.com/vincentarelbundock/QReg
* Author: Vincent Arel-Bundock
* Contact: varel@umich.edu
* Date: 2013-07-13
* License: BSD-3
* Original code from the python statsmodels project
    - https://github.com/statsmodels/statsmodels

# Example

The file ``examples/qreg_example.jl`` shows how to use the functions provided here. It replicates part of the analysis in:

* Koenker, Roger and Kevin F. Hallock. "Quantile Regression". Journal of Economic Perspectives, Volume 15, Number 4, Fall 2001, Pages 143–156

We are interested in the relationship between income and expenditures on food for a sample of working class Belgian households in 1857 (the Engel data), and plot regression estimates at different quantiles.

![](https://github.com/vincentarelbundock/QReg/raw/master/examples/qreg_example_plot.png)
