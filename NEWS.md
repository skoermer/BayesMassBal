# BayesMassBal 1.0.0

## Major Changes

* Addition of output from using the output of the summary function on a "BayesMassBal" object
* ssEst function allows for fitting of a Bayesian lag-1 autoregressive model to a set of mass flow observations so that average behavior time to steady state can be inferred
* In BMB where independent variance is specified, inference is conducted on variance instead of precision
* Change in default prior specification for BMB function, ability to specify Jeffreys priors

# BayesMassBal 0.2.1

## Major Changes

* Markov chain Monte Carlo diagnostics are now available directly from the BMB function

## Bug Fixes

* Removed MASS from the Imports section in DESCRIPTION, as no functions from the MASS package are used.

# BayesMassBal 0.2.0

* Initial release
