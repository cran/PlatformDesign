# testing functions

test_that("test cor.mat()", {
  res <- cor.mat(K = 2, M = 0, n = 101, n0 = 143)
  expect_equal(res$cor1, 0.4139344, tolerance=1e-6)
})

test_that("test fwer_critical()", {
  corMat1 <- cor.mat(K=2, M = 2, n=107, n0=198, n0t = 43)$cormat
  res <- fwer_critical(ntrt=4, fwer=0.025, corMat=corMat1)
  expect_equal(res$critical_val, 2.475233, tolerance=1e-6)

})

test_that("test one_stage_multiarm()", {
  res <- one_stage_multiarm(K = 2, fwer = 0.025, marginal.power = 0.8, delta = 0.4)
  expect_equal(res$n1, 101)
  expect_equal(res$n0_1, 143)
  expect_equal(res$N1, 345)
  expect_equal(res$z_alpha1, 2.220604, tolerance=1e-6)
  expect_equal(res$Power1, 0.9222971, tolerance=1e-6)

  res2 <-one_stage_multiarm(K = 1, fwer = 0.025, marginal.power = 0.8, delta = 0.4)
  expect_equal(res2$z_alpha1, 1.959964, tolerance=1e-6)

})

test_that("test one_design()", {
  res <- one_design(n2 = 107, n0_2 = 198, nt = 30, K = 2, M = 2, fwer = 0.025, marginal.power = 0.8, delta = 0.4)
  expect_equal(res$design_Karm$n1, 101)
  expect_equal(res$design_Karm$n0_1, 143)
  expect_equal(res$designs$critical_value2, 2.475233, tolerance=1e-6)
  expect_equal(res$designs$marginal.power2, 0.80011, tolerance=1e-6)
  expect_equal(res$designs$disjunctive.power2, 0.9853799, tolerance=1e-4)
})

test_that("test one_design2()", {
  K=2
  M=2
  nt <- 30
  A1 <- sqrt(K)
  n0t <- ceiling(A1 * nt)

  ntrt <- K+M
  multi <- one_stage_multiarm(K=2, fwer=0.025, marginal.power=0.8, delta=0.4)
  n1 <- multi$n1
  n0_1 <- multi$n0_1
  N1 <- multi$N1
  corMat1 <- multi$corMat1
  z_beta1 <- multi$z_beta1
  Power1 <- multi$Power1
  z_alpha1 <- multi$z_alpha1

  res <- one_design2(K = 2, M = 2,
    n1 = n1, n0_1 = n0_1, n2 = 107, n0_2 = 198,
    nt = nt, n0t = n0t, A1,
    fwer = 0.025,
    z_alpha1 = z_alpha1, z_beta1 = z_beta1)

expect_equal(res$N2, 669)
  expect_equal(res$A2, 2.012987, tolerance=1e-6)
  expect_equal(res$critical_value2, 2.475233, tolerance=1e-6)
  expect_equal(res$marginal.power2, 0.80011, tolerance=1e-6)
  expect_equal(res$disjunctive.power2, 0.9853799, tolerance=1e-5)

})
