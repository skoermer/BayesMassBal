## Test environments

* Local Ubuntu 20.04.4 Focal, R 4.2.0
* Ubuntu Linux 20.04.1 LTS, R-devel, GCC
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* macOS 10.13.6 High Sierra, R-release, CRAN's setup
* Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
* Ran devtools::check_win_devel() with defaults
* Ran devtools::check_rhub() with defaults

## R CMD check results (local)

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### Notes

*  Known R-hub issue #503 https://github.com/r-hub/rhub/issues/503 : 

checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'

## Downstream dependencies

There are currently no downstream dependencies for this package.
