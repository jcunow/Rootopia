# Standard testthat entry point. `R CMD check` and devtools::test() run this.
library(testthat)
library(Rootopia)        # <- package name; change if your DESCRIPTION differs

test_check("Rootopia")
