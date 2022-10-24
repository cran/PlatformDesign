##' Provide design parameters based on n1, n0_1, n2, n0_2, while controlling for FWER
##'
##'Given n1, n0_1, n2 and n0_2, provide the design parameters for one single design, while controlling for FWER
##' @title Calculate design parameters based on user-defined n, n0, n2 and n0_2 for adding M new experimental arms
##' @param K the number of experimental arms in the original clinical trial.
##' @param M a positive integer, the number of newly added experimental arms. When M = 0, the function returns the design of a K experimental arm trial.
##' @param n1 a positive integer, which is the sample size in each of the K experimental arms in the original clinical trial.
##' @param n0_1 a positive integer, which is the sample size of the controls for each of the K experimental arms in the original clinical trial.
##' @param n2 a positive integer, which is the updated sample size in each of the K+M experimental arms. n2 > nt.
##' @param n0_2 a positive integer, which is the updated sample size of the controls for each of the K+M experimental arms.
##' @param nt a positive integer, which is the number of patients already enrolled on each of the original K experimental arms when the new arms are added. nt should be smaller than n1 and n2.
##' @param n0t the number of controls already enrolled in the control arm when new arms are added
##' @param fwer the family-wise error rate, kept constant for original and the adapted design.
##' @param z_alpha1 the critical value of the original K-experimental arm design
##' @param z_beta1 qnorm(1-beta1)
##' @param seed an integer for random number generation for numerically evaluating integration, default = 123
##' @noRd
##'
##' @import mvtnorm
##' @examples
##' one_design2(K=2, M=2, n1=101, n0_1=143, n2=107, n0_2=198, nt=30, n0t=43, fwer=0.025,
##'  z_alpha1=2.220604, z_beta1=0.8416212,seed=123)
##' #$n2
##' #[1] 107
##' #
##' #$n0_2
##' #[1] 198
##' #
##' #$nt
##' #[1] 30
##' #
##' #$n0t
##' #[1] 43
##' #
##' #$nc
##' #[1] 241
##' #
##' #$N2
##' #[1] 669
##' #
##' #$A2
##' #[1] 2.012987
##' #
##' #$cor1
##' #[1] 0.3508197
##' #
##' #$cor2
##' #[1] 0.2746316
##' #
##' #$critical_value2
##' #[1] 2.475233
##' #
##' #$marginal.power2
##' #[1] 0.8001101
##' #
##' #$disjunctive.power2
##' #[1] 0.9853799
##' #
##' one_design2(K=1, M=3, n1=99, n0_1=99, n2=105, n0_2=204, nt=30, n0t=30, fwer=0.025, z_alpha1=1.959964, z_beta1=0.8416212,seed=123)

one_design2 <- function(K, M, n1, n0_1, n2, n0_2, nt, n0t, fwer, z_alpha1, z_beta1,
                        seed=123) {

    ntrt <- K + M
    A2 <- (n0_2 - n0t)/(n2 - nt)
    nc <- n0t + n0_2  # sample size for control arm in new design
    N2 <- nc + ntrt * n2  # total sample size

    c2 <- cor.mat(K = K, M = M, n = n2, n0 = n0_2, n0t = n0t)
    corMat2 <- c2$cormat

    if (is.null(fwer)==0){
        z_alpha2 <- fwer_critical(ntrt, fwer, corMat2)$critical_val
        }
     else{
         z_alpha2 <- z_alpha1
         }

    a <- sqrt((1/n1 + 1/n0_1)/(1/n2 + 1/n0_2))
    z_beta2 <- a * (z_alpha1 + z_beta1) - z_alpha2
    marginal.power2 = pnorm(z_beta2)  # empirical marginal power, slightly differs from vzbeta, might due to rounding
    set.seed(seed)
    Power2 <- (1 - pmvnorm(lower = rep(-Inf, ntrt), upper = rep(-z_beta2, ntrt), mean = rep(0, ntrt), corr = corMat2))[1]

    res <- list(n2 = n2, n0_2 = n0_2,
                nt = nt, n0t = n0t,
                nc = nc, N2 = N2,
                A2 = A2,
                cor1 = c2$cor1, cor2 = c2$cor2,
                critical_value2 = z_alpha2,
                marginal.power2 = marginal.power2,
                disjunctive.power2 = Power2)


    return(res)

}

