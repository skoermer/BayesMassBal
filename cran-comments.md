## Test environments
* local windows install, R 4.0.2
* ubuntu 16.04 xenial (on travis-ci), R 3.6.3
* ubuntu 16.04 xenial (on travis-ci), R 4.0.0
* ubuntu 16.04 xenial (on travis-ci), R devel
* OSX (on travis-ci), R 3.6.3
* OSX (on travis-ci), R 4.0.2
* Windows Server (on appveyor), R 4.0.2
* Ran devtools::check_rhub() with defaults
* Ran devtools::check_win_devel() with defaults

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Downstream dependencies

There are currently no downstream dependencies for this package.


## Resubmission
This is a resubmission. In this version I have:

* Resets par() to the initial values in the last code chunk in the Two_Node_Process vignette.
* Resets par() to the initial values upon exiting the plot.BayesMassBal function
* Adds a test to check par() has been reset when using the plot.BayesMassBal function

* Changes commented out code in examples to code within \donttest{} and \dontrun{}

* Replaces the "MCMC"" acronym with "Markov chain Monte Carlo" in the Description field within DESCRIPTION
* The first use of the acronym "MCMC" or "HPDI" in any individual document is proceeded by the meaning of the acronym.

* The reference cited for approximating log-marginal likelihood is included in the Description field within DESCRIPTION using the author name, publication year, and <doi:...>.
