library(testthat)
library(devtools)
load_all('../brms')
test_check("brms", path='../brms/tests')
