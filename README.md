
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesMassBal

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/skoermer/BayesMassBal.svg?branch=master)](https://travis-ci.org/skoermer/BayesMassBal)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/skoermer/BayesMassBal?branch=master&svg=true)](https://ci.appveyor.com/project/skoermer/BayesMassBal)
[![Codecov test
coverage](https://codecov.io/gh/skoermer/BayesMassBal/branch/master/graph/badge.svg)](https://codecov.io/gh/skoermer/BayesMassBal?branch=master)
[![Travis build
status](https://travis-ci.com/skoermer/BayesMassBal.svg?branch=master)](https://travis-ci.com/skoermer/BayesMassBal)
<!-- badges: end -->

The goal of BayesMassBal is to allow users to easily conduct Bayesian
data reconciliation for a linearly constrained chemical or particulate
process at steady state.

Samples taken from a chemical process are always observed with noise.
Using data reconciliation, or mass balance methods, it is possible to
use the principle of [conservation of
mass](https://en.wikipedia.org/wiki/Conservation_of_mass) to filter the
noise. This technique is common in chemical engineering and mineral
processing engineering applications.

Typically, a mass balance produces point estimates of true mass flow
rates. However, using Bayesian methods one can obtain a more granular
view of process uncertainty. The `BayesMassBal` package provides
functions allowing the user to easily specify conservation of mass
constraints, organize collected data, conduct a Bayesian mass balance
using various error structures, and select the best model for their data
using [Bayes Factors](https://en.wikipedia.org/wiki/Bayes_factor).

The Bayesian mass balance uses Markov chain Monte Carlo methods to
obtain random samples from the distributions of constrained mass flow
rates. These samples can be used to generate plots, or for other
applications where sampling from such a distribution is useful.

## Installation

You can install the released version of BayesMassBal from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("BayesMassBal")
```

## Using `BayesMassBal`

After loading the package

``` r
library(BayesMassBal)
```

Functions are available to aid in Bayesian data reconciliation.

  - The `importObservations()` function can be used to import mass flow
    rate data from a `*.csv` file into `R` and organize it for use with
    the `BayesMassBal` package.
  - Toy data sets can be simulated using the `twonodeSim()` function for
    educational purposes, or for comparing the performance of Bayesian
    data reconciliation methods to other methods.
  - Using the `constrainProcess()` function, one can specify linear
    constraints in `R` or import them from a `*.csv` file.
  - The **B**ayesian **M**ass **B**alance function, `BMB()`, then can be
    used to generate samples from target distributions and approximate
    the log marginal likelihood for a specified model.
  - A summary table can be viewed in the console and saved using
    `summary.BayesMassBal()`.
  - The output from `BMB()` is a `"BayesMassBal"` object, which can be
    fed to `plot.BayesMassBal()` to easily plot the results.
  - A `"BayesMassBal"` object can also be used with the `BayesMassBal`
    function `mainEff()` to inspect how the main effect of a random
    variable and uncertainty in process performance are related.

An overview of a suggested workflow, including importing data into `R`,
specifying model constraints, using the `BMB` function, and making a
main effects plot, is available as a vignette:
`vignette("Two_Node_Process")`.
