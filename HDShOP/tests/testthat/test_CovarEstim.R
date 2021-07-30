
n<-3e2 # number of realizations
p<-.5*n # number of assets

x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)

Mtrx_trad <- CovarEstim(x, type="trad")

TM <- matrix(0, p, p)
diag(TM) <- 1

Mtrx_bgp <- CovarEstim(x, type="BGP14", TM=TM)
Mtrx_bgp_2 <- CovarEstim(x, type="BGP14", SCM=TM, TM=TM)
Mtrx_lw <- CovarEstim(x, type="LW20")



#### tests on matrices
test_that("CovarEstim works on matrices", {

  expect_is(Mtrx_trad, class="matrix")
  expect_is(Mtrx_bgp, class="matrix")
  expect_is(Mtrx_bgp_2, class="matrix")
  expect_is(Mtrx_lw, class="matrix")
})


#### reset the input data

x <- as.data.frame(x)

Mtrx_trad <- CovarEstim(x, type="trad")

TM <- matrix(0, p, p)
diag(TM) <- 1
TM <- as.data.frame(TM)

Mtrx_bgp <- CovarEstim(x, type="BGP14", TM=TM)
Mtrx_bgp_2 <- CovarEstim(x, type="BGP14", SCM=TM, TM=TM)
Mtrx_lw <- CovarEstim(x, type="LW20")



#### tests on data frames
test_that("CovarEstim works on data frames", {

  expect_is(Mtrx_trad, class="matrix")
  expect_is(Mtrx_bgp, class="matrix")
  expect_is(Mtrx_bgp_2, class="matrix")
  expect_is(Mtrx_lw, class="matrix")
})




