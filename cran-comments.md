## Resubmission
This is a resubmission. In this version I have: 

* Used single quotes for software names ('RevBayes') in package DESCRIPTION
* Added \value fields to:
  * densiTreeWithBranchData.Rd
  * plotAncStatesMAP.Rd
  * plotAncStatesPie.Rd
  * processAncStates.Rd
* Changed all calls from cat() to either stop() or close(bar)
* Ensured all tests plot to a pdf in a temp directory rather than package directory
* Put all examples that download data within \donttest{} statements


## Release summary
This is a new submission to CRAN.

## Test environments
Tested via GitHub actions on the following platforms (and R versions):
* macOS-latest (release R)
* windows-latest (release R)
* ubuntu-18.04 (release R)
* ubuntu-18.04 (devel R)
* ubuntu-20.04 (release R)
* ubuntu-20.04 (devel R)

Tested locally on OS X v. 10.15.7, R v. 4.1.0

Tested on Windows, Ubuntu, and Fedora with devtools::check_rhub()

Tested on Window with WindBuilder using devtools::check_win_release()

## R CMD check results 
There were no ERRORs or WARNINGs. 

There was one NOTE:
* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Carrie Tribble <ctribble09@gmail.com>’

  New submission

This is a new submission to CRAN.

## Downstream dependencies
There are currently no downstream dependencies for this package. 
