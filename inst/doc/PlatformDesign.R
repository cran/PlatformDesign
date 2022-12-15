## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning = FALSE---------------------------------------------------
library(PlatformDesign)
library(ggplot2)

## -----------------------------------------------------------------------------
K <- 2
FWER_1 <- 0.025
beta1 <- 0.2
z_beta1 <- qnorm(1-beta1) #z_(1-beta1)
A1 <- sqrt(K)

## -----------------------------------------------------------------------------
multi <- one_stage_multiarm(K = 2, fwer = 0.025, marginal.power = 0.8, delta = 0.4)
corMat1 <- multi$corMat1
corMat1

## -----------------------------------------------------------------------------
multi$z_alpha1

## -----------------------------------------------------------------------------
multi$n1
multi$n0_1
multi$N1

## -----------------------------------------------------------------------------
multi$Power1

## -----------------------------------------------------------------------------
multi

## -----------------------------------------------------------------------------
nt <- 30
nt

n_0t <- ceiling(nt*A1)
n_0t

## -----------------------------------------------------------------------------
FWER_2 <- FWER_1
FWER_2
omega2_min <- 1-beta1
omega2_min
Omega2_min <- multi$Power1
Omega2_min

## -----------------------------------------------------------------------------
pair3 <- admiss(n1=101, n0_1=143, nt=30, ntrt=4, S=690)

ggplot(data=pair3, aes(x=Var1, y=Var2)) +
  geom_point() +
  geom_abline(intercept = 647, slope=-4, color="red") +
  geom_hline(yintercept=43, color="red")+
  geom_vline(xintercept=30, color="red") +
  xlim(0, 500)+
  ylim(0,1000)+
  xlab("n2")+
  ylab("n02")

## ---- eval=F------------------------------------------------------------------
#  design <- platform_design(nt = 30, K = 2, M = 2, fwer = 0.025, marginal.power = 0.8, delta = 0.4)
#  design$designs

## ---- eval=F------------------------------------------------------------------
#  design2 <- platform_design(nt=30, K=2, M=2, pwer=0.025, marginal.power=0.8, delta=0.4,seed=123)
#  design2$designs

## ----eval=F-------------------------------------------------------------------
#  start_time <- Sys.time()
#  test <- platform_design2(nt = 30, K = 2, M = 2, fwer = 0.025, marginal.power = 0.8, delta = 0.4)
#  end_time <- Sys.time()
#  end_time - start_time
#  # Time difference of 41.85487 secs

## ----eval=F-------------------------------------------------------------------
#  start_time <- Sys.time()
#  test2 <- platform_design(nt = 30, K = 2, M = 2, fwer = 0.025, marginal.power = 0.8, delta = 0.4)
#  end_time <- Sys.time()
#  end_time - start_time
#  # Time difference of 8.188013 mins

## -----------------------------------------------------------------------------
start_time <- Sys.time()
platform_design2(nt = 50, K = 2, M = 2, fwer = 0.025, marginal.power = 0.8, delta = 0.4)
end_time <- Sys.time()

# Time difference
end_time - start_time

