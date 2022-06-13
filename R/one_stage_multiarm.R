##' This design is to find the required sample sizes and the associated critical value to control the overall
##' type I error rate (FWER) and achieve the user-specified marginal (i.e., experimental-wise or pairwise)
##' power. Calculate required sample sizes for each of the experimental arm (n1), the control arm (n0_1),
##' the total sample size (N1) and the critical value z_alpha1 in a K-arm trial setting (K experimental arms
##' and 1 common control arm).
##'
##' Given the number of experimental arm (K), the family-wise type I error rate, the
##' marginal power for each experimental-control comparison and the standardized effect
##' size to calculate the sample size and other design parameters for the K-experimental
##' arm trial with one control arm.
##'
##' @title Calculate the sample size and other design parameters for an one-stage K-arm trial using the root-K
##'        rule for the allocation ratio
##' @param K the number of experimental arms
##' @param fwer the family-wise type I error rate
##' @param marginal.power the marginal power for each experimental-control comparison
##' @param delta the standardized clinical effect size expected to be detected
##'        in the trial
##' @param seed an integer for random number generation for numerically evaluating integration, default = 123
##'
##' @return \emph{n1} the sample size of each of the K experimental arms
##' @return \emph{n0_1} the sample size of the control arm
##' @return \emph{N1} the total sample size of a K-arm trial
##' @return \emph{z_alpha1} the critical value for the comparison between any of the K-experimental arm
##'              in the first period and its corresponding control
##' @return \emph{z_beta1} the value of the quantile function of the standard normal distribution with
##'              probability = marginal power of the K-arm trial
##' @return \emph{Power1} the disjunctive power of the K-arm trial defined as the
##'              probability of rejecting at least one of the K experimental arms under the
##'              alternative hypothesis
##' @return \emph{corMat1} the correlation matrix of the Z-test statistics
##'
##' @import mvtnorm
##' @export
##'
##' @author \verb{    }Xiaomeng Yuan, Haitao Pan
##'
##' @references \verb{    }Pan, H., Yuan, X. and Ye, J. (2022). An optimal two-period multi-arm confirmatory platform design
##' with adding new arms. Manuscript submitted for publication.
##' @references \verb{    }Dunnett, C. W. (1955). A multiple comparison procedure for comparing
##' several treatments with a control. Journal of the American Statistical
##' Association, 50(272), 1096-1121.
##'
##' @examples
##' one_stage_multiarm(K = 2, fwer = 0.025, marginal.power = 0.8, delta = 0.4)
##' #$n1
##' #[1] 101
##'
##' #$n0_1
##' #[1] 143
##'
##' #$N1
##' #[1] 345
##'
##' #$z_alpha1
##' #[1] 2.220604
##'
##' #$z_beta1
##' #[1] 0.8416212
##'
##' #$Power1
##' #[1] 0.9222971
##'
##' #$corMat1
##' #[,1]      [,2]
##' #[1,] 1.0000000 0.4142136
##' #[2,] 0.4142136 1.0000000
##'

one_stage_multiarm <- function(K, fwer, marginal.power, delta, seed=123) {

    z_beta1 <- qnorm(marginal.power)  #z_(1-beta1)
    A1 <- sqrt(K)

    rho0 <- 1/(A1 + 1)

    cor.mat2 <- function(K, cor1) {

        cormat <- matrix(0, nrow = K, ncol = K)

        for (i in 1:K) {
            for (j in 1:K) {
                if (i != j) {
                  cormat[i, j] <- cor1
                } else cormat[i, j] <- 1
            }
        }

        cormat
    }

    corMat1 <- cor.mat2(K = K, cor1 = rho0)

    fwer.old <- fwer_critical(K, fwer, corMat = corMat1, seed)
    z_alpha1 <- fwer.old$critical_val  # z_(1-alpha1)

    if (K == 1) {
        Power1 = marginal.power
    } else {

        set.seed(seed)
        Power1 <- 1 - pmvnorm(lower = rep(-Inf, K), upper = rep(-z_beta1, K), mean = rep(0, K), corr = corMat1)
    }

    n1 <- ceiling((z_alpha1 + z_beta1)^2/delta^2 * (1 + sqrt(K)/K))
    n0_1 <- ceiling(n1 * sqrt(K))
    N1 <- K * n1 + n0_1

    res <- list(n1 = n1, n0_1 = n0_1, N1 = N1, z_alpha1 = z_alpha1, z_beta1 = z_beta1, Power1 = Power1[1], corMat1 = corMat1)

    return(res)
}


