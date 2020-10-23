context("Testing par is reset when using plot() on a S3 BayesMassBal object")

test_that("par is reset within plot.BayesMassBal",{

  par_initial <-  par(no.readonly =TRUE)

  y <- importObservations(file = system.file("extdata", "twonode_example.csv",
                                               package = "BayesMassBal"),
                            header = TRUE, csv.params = list(sep = ";"))

  C <- matrix(c(1,-1,0,-1,0,0,1,-1,0,-1), byrow = TRUE, ncol = 5, nrow = 2)
  X <- constrainProcess(C = C)

  BMB_example <- BMB(X = X, y = y, cov.structure = "indep",
                     BTE = c(10,50,1), lml = FALSE, verb=0)

  plot(BMB_example,sample.params = list(beta = list(CuFeS2 = 1:3, gangue = 1:3)),layout = "trace",hdi.params = c(1,0.95))

  par_final <- par(no.readonly =TRUE)

  expect_that(identical(par_initial, par_final), equals(TRUE))
})
