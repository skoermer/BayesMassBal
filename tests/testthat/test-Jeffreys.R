test_that("Jeffreys priors selection works", {

  C <- matrix(c(1,-1,0,-1,0,0,1,-1,0,-1), byrow = TRUE, ncol = 5, nrow = 2)
  X <- constrainProcess(C = C)

  y <- twonodeSim()$simulation

  expect_silent(BMB(X = X, y = y,BTE= c(1,10,1), priors = "Jeffreys", lml = FALSE, verb= 0))
  expect_silent(BMB(X = X, y = y,BTE= c(1,10,1),cov.structure = "location", priors = "Jeffreys", lml = FALSE, verb= 0))
  expect_silent(BMB(X = X, y = y,BTE= c(1,10,1),cov.structure = "component", priors = "Jeffreys", lml = FALSE, verb= 0))

})
