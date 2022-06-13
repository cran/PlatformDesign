##' Provide other design parameters for a two-period K+M trial, given n2 and n0_2, nt, K, M, fwer and marginal
##' power (of the the first period). This function serves for the purpose of spot-testing for any
##' pre-specified {n, n0_2} pair. Please use \emph{platform_design()} for finding optimal values of n and n0_2.
##'
##' Given n2 and n0_2, nt, K, M, fwer and marginal power (of the first period), provide other design
##' parameters for a two-period K+M trial.
##'
##' @title Calculate other design parameters of a two-period multi-arm platform design given updated sample sizes
##' @param n2 a positive integer, which is the sample size in each of the K+M experimental arms in the second period, n2 > nt
##' @param n0_2 a positive integer, which is the sample size of the concurrent control for each of the K+M experimental arms in the in the second period
##' @param nt a positive integer, the number of patients already enrolled on each of the K experimental arms in the first period when the new arms are added
##' @param K a positive integer, the number of experimental arms in the  first period in a two-period K+M trial
##' @param M a positive integer, the number of delayed (newly added) experimental arms added in the second period
##' @param fwer the family-wise error rate (FWER) to be controlled, default to be the same
##'        throughout the trial
##' @param marginal.power the marginal power to achieve in the first period in a two-period
##'        K+M trial
##' @param delta the standardized clinical effect size expected to be detected in
##'        the trial
##' @param seed an integer for random number generation for numerically evaluating integration, default = 123
##'
##' @return \bold{designs} contains the calculated design parameters for period 1 and 2 including:
##' @return \verb{    }\emph{n1} and \emph{n0_1}, the sample sizes of each of the K experimental arms and the
##'              control arm, respectively, in the first period
##' @return \verb{    }\emph{n2} and \emph{n0_2}, the updated sample sizes of each of the K + M experimental arms and
##'              its corresponding concurrent control, respectively, after adding M experimental arms in the second period
##' @return \verb{    }\emph{nt} and \emph{n0t}, the number of patients already enrolled on each of the K
##'              experimental arms and the control arm, respectively, in the first period when the new
##'              arms are added
##' @return \verb{    }\emph{nc}, the updated total sample size of the control arm after adding M
##'          experimental arms in the second period, i.e., the sum of concurrent (n0_2) and nonconcurrent
##'          (n0t) controls
##' @return \verb{    }\emph{N2}, the total sample size of the two-period K+M experimental arm (and 1 control
##'               arm) platform trial
##' @return \verb{    }\emph{A1},  the allocation ratio (control to experimental arm) before the M new
##'               experimental arms are added and after the initial K experimental arms end
##' @return \verb{    }\emph{A2},  the allocation ratio after the M new experimental arms are added and before
##'               the initial K experimental arms end
##' @return \verb{    }\emph{cor1}, the correlation of Z statistics between any two of the K initially opened
##'               experimental arms (or between any two of the M delayed arms)
##' @return \verb{    }\emph{cor2}, the correlation of Z statistics between any pair of one initially opened
##'               and one delayed experimental arm
##' @return \verb{    }\emph{critical_value1}, the critical value for the comparison between any of the K
##'              experimental arms in the first period and the corresponding control
##' @return \verb{    }\emph{critical_value2}, the critical value for the comparison between any of the K + M
##'              experimental arms in the second period and the corresponding control
##' @return \verb{    }\emph{marginal.power1} and \emph{marginal.power2}, the marginal power for the first
##'              and second period, respectively
##' @return \verb{    }\emph{disjunctive.power1} and \emph{disjunctive.power2}, the disjunctive power for the
##'               first and second period, respectively
##' @return \verb{    }\emph{effect_size}, the standardized clinical effect size
##'               expected to be detected in the trial
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
##' one_design(n2 = 107, n0_2 = 198, nt = 30, K = 2, M = 2, fwer = 0.025,
##'            marginal.power = 0.8, delta = 0.4)
##'
##' #$n1
##' #[1] 101
##'
##' #$n0_1
##' #[1] 143
##'
##' #$n2
##' #[1] 107
##'
##' #$n0_2
##' #[1] 198
##'
##' #$nt
##' #[1] 30
##'
##' #$n0t
##' #[1] 43
##'
##' #$nc
##' #[1] 241
##'
##' #$N2
##' #[1] 669
##'
##' #$A1
##' #[1] 1.414214
##'
##' #$A2
##' #[1] 2.012987
##'
##' #$cor1
##' #[1] 0.3508197
##'
##' #$cor2
##' #[1] 0.2746316
##'
##' #$critical_value1
##' #[1] 2.220604
##'
##' #$critical_value2
##' #[1] 2.475233
##'
##' #$marginal.power1
##' #[1] 0.8
##'
##' #$marginal.power2
##' #[1] 0.80011
##'
##' #$disjunctive.power1
##' #[1] 0.9222971
##'
##' #$disjunctive.power2
##' #[1] 0.9853799
##'
##' #$effect_size
##' #[1] 0.4
##'

one_design <- function(n2, n0_2, nt, K, M, fwer, marginal.power, delta, seed=123) {
    if (sum(c(n2, n0_2, nt, K, M)%%1) != 0) {
        stop("n2, n0_2, nt, K and M should all be integers.")
    }
    if (sum(c(n2, n0_2, nt, K, M) > 0) != length(c(n2, n0_2, nt, K, M))) {
        stop("n1, n0_1, n2, n0_2, nt, K and M should all > 0.")
    }
    if (n2 < nt)
        stop("n2 has to be greater than nt, use admiss() to choose from admissible set")

    ntrt <- K + M
    A1 <- sqrt(K)
    n0t <- ceiling(A1 * nt)
    A2 <- (n0_2 - n0t)/(n2 - nt)
    nc <- n0t + n0_2  # sample size of the control arm in new design
    N2 <- nc + ntrt * n2  # total sample size in new design

    if (n0_2 <= n0t)
        stop("n0_2 <= n0t, use admiss() to choose from admissible set")

    # original design parameters
    multi <- one_stage_multiarm(K, fwer, marginal.power, delta, seed=seed)
    n1 <- multi$n1
    if (nt >= n1)
        stop("nt has to be smaller than n1")
    n0_1 <- multi$n0_1
    N1 <- multi$N1
    corMat1 <- multi$corMat1

    if (K == 1) {
        z_alpha1 <- qnorm(1 - fwer)
    } else {
        z_alpha1 <- multi$z_alpha1
    }

    z_beta1 <- multi$z_beta1
    Power1 <- multi$Power1

    # update design parameters
    c2 <- cor.mat(K = K, M = M, n = n2, n0 = n0_2, n0t = n0t)
    corMat2 <- c2$cormat
    z_alpha2 <- fwer_critical(ntrt = ntrt, fwer, corMat2, seed)$critical_val
    a <- sqrt((1/n1 + 1/n0_1)/(1/n2 + 1/n0_2))
    z_beta2 <- a * (z_alpha1 + z_beta1) - z_alpha2
    marginal.power2 = pnorm(z_beta2)

    set.seed(seed)
    Power2 <- 1 - pmvnorm(lower = rep(-Inf, ntrt), upper = rep(-z_beta2, ntrt), mean = rep(0, ntrt), corr = corMat2)[1]

    res <- list(n1 = n1, n0_1 = n0_1, n2 = n2, n0_2 = n0_2, nt = nt, n0t = n0t, nc = nc, N2 = N2, A1 = A1, A2 = A2,
        cor1 = c2$cor1, cor2 = c2$cor2, critical_value1 = z_alpha1, critical_value2 = z_alpha2, marginal.power1 = marginal.power,
        marginal.power2 = marginal.power2, disjunctive.power1 = Power1, disjunctive.power2 = Power2, effect_size = delta)


    return(res)

}


