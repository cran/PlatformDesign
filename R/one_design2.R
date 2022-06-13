##' Provide design parameters based on n, n0, n2, n0_2
##'
##'Given n, n0, n2 and n0_2, provide the design parameters for one single design
##' @title Calculate design parameters based on user-defined n, n0, n2 and n0_2 for adding M new experimental arms
##' @param K the number of experimental arms in the original clinical trial.
##' @param  M a positive integer, the number of newly added experimental arms. When M = 0, the function returns the design of a K experimental arm trial.
##' @param n1 a positive integer, which is the sample size in each of the K experimental arms in the original clinical trial.
##' @param n0_1 a positive integer, which is the sample size of the controls for each of the K experimental arms in the original clinical trial.
##' @param n2 a positive integer, which is the sample size in each of the K+M experimental arms in the newly adapted trial. n2 > nt.
##' @param n0_2 a positive integer, which is the sample size of the controls for each of the K+M experimental arms in the newly adapted trial.
##' @param nt a positive integer, which is the number of patients already enrolled on each of the original K experimental arms when the new arms are added. nt should be smaller than n1 and n2.
##' @param n0t the number of controls already enrolled in the control arm when new arms are added
##' @param A1 sqrt(K)
##' @param  fwer the family-wise error rate, kept constant for original and the adapted design.
##' @param Power1 the disjunctive power of the original K-experimental arm design
##' @param z_alpha1 the critical value of the original K-experimental arm design
##' @param z_beta1 qnorm(1-beta1)
##' @param marginal.power the marginal power fore each control-experimental arm comparison in the original design.
##' @param delta standardized effect size
##' @param seed an integer for random number generation for numerically evaluating integration, default = 123
##' @noRd
##'
##' @import mvtnorm


one_design2 <- function(K, M, n1, n0_1, n2, n0_2, nt, n0t, A1, fwer, Power1, z_alpha1, z_beta1,
                        marginal.power, delta, seed=123) {

    ntrt <- K + M
    A2 <- (n0_2 - n0t)/(n2 - nt)
    nc <- n0t + n0_2  # sample size for control arm in new design
    N2 <- nc + ntrt * n2  # total sample size

    c2 <- cor.mat(K = K, M = M, n = n2, n0 = n0_2, n0t = n0t)
    corMat2 <- c2$cormat
    z_alpha2 <- fwer_critical(ntrt, fwer, corMat2)$critical_val

    a <- sqrt((1/n1 + 1/n0_1)/(1/n2 + 1/n0_2))
    z_beta2 <- a * (z_alpha1 + z_beta1) - z_alpha2
    marginal.power2 = pnorm(z_beta2)  # empirical marginal power, slightly differs from vzbeta, should be due to rounding?
    set.seed(seed)
    Power2 <- (1 - pmvnorm(lower = rep(-Inf, ntrt), upper = rep(-z_beta2, ntrt), mean = rep(0, ntrt), corr = corMat2))[1]

    res <- list(n1 = n1, n0_1 = n0_1, n2 = n2, n0_2 = n0_2, nt = nt, n0t = n0t, nc = nc, N2 = N2, A1 = A1, A2 = A2,
        cor1 = c2$cor1, cor2 = c2$cor2, critical_value1 = z_alpha1, critical_value2 = z_alpha2, marginal.power1 = marginal.power,
        marginal.power2 = marginal.power2, disjunctive.power1 = Power1, disjunctive.power2 = Power2, effect_size = delta)


    return(res)

}

