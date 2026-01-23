if(!require('survRM2')) install.packages('survRM2')

library('testthat')

# Create test data
set.seed(1234)
time0 <- rgeom(200, 0.06)
time1 <- rgeom(200, 0.053)
time <- c(time0, time1)
status <- rbinom(400, 1, 0.3)
arm <- c(rep(0, 200), rep(1, 200))
tau <- 40

withr::defer({
  detach(package:survRM2)
}, teardown_env())
