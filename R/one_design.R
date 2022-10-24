##' Provide other design parameters for a two-period K+M-experimental arm trial, given n2 and
##' n0_2, nt, K, M, fwer(or pwer) and marginal power (of the K-experimental arm trial). This
##' function serves for the purpose of spot-testing for any pre-specified {n,
##' n0_2} pair. Please use \emph{platform_design()} for finding optimal values
##' of n and n0_2, controlling for FWER (or PWER).
##'
##' Given n2 and n0_2, nt, K, M, fwer (or pwer) and marginal power (of the K-experimental arm trial), provide
##' other design parameters for a two-period K+M-experimental arm trial.
##'
##' @title Calculate other design parameters of a two-period K+M-experimental arm platform design given
##'  sample sizes
##' @param n2 a positive integer, which is the sample size in each experimental arm
##'  in the K+M-experimental arm trial
##' @param n0_2 a positive integer, which is the sample size of the concurrent control for each
##'  experimental arms in the K+M-experimental arm trial
##' @param nt a positive integer, the number of patients already enrolled on each of the K initial
##'  experimental arms when the M new arms are added
##' @param K a positive integer, the number of experimental arms in the first period in a
##'  two-period K+M-experimental arm trial
##' @param M a positive integer, the number of delayed (newly added) experimental arms added in
##'  the beginning of the second period of the K+M-experimental arm trial
##' @param fwer the family-wise error rate (FWER) to be controlled, default to be the same
##'        throughout the trial
##' @param pwer the pair-wise error rate (PWER) to be controlled, default to be the same
##'        throughout the trial
##' @param marginal.power the marginal power to achieve in the K-experimental arm (and K+M-experimental arm) trial
##' @param delta the standardized clinical effect size expected to be detected in
##'        the trial
##' @param seed an integer for random number generation for numerically evaluating integration,
##'  default = 123
##'
##' @return \bold{design_Karm} contains the  design parameters for the K-experimental arm
##'  trial including:
##' @return \verb{    }\emph{K}, the number of experimental arms in the K-experimental arm trial
##' @return \verb{    }\emph{n1}, the sample size for each of the K experimental
##'  arms in the k-experimental arm trial
##' @return \verb{    }\emph{n0_1}, the sample size of the common control arm in the K-experimental arm trial
##' @return \verb{    }\emph{N1} the total sample size of a K-experimental arm trial
##' @return \verb{    }\emph{z_alpha1}, the critical value for the comparison between any of the K
##'  experimental arms and the control in the K-experimental arm trial
##' @return \verb{    }\emph{FWER1}, the family-wise error rate for the K-experimental arm trial
##' @return \verb{    }\emph{z_beta1}, the quantile of the marginal power, i.e., qnorm(marginal
##'  power) for the K-experimental arm trial
##' @return \verb{    }\emph{Power1}, the disjunctive power for the K-experimental arm trial
##' @return \verb{    }\emph{cor0}, the correlation of Z-test statistics between any two of the K
##'  experimental arms
##' @return \verb{    }\emph{delta}, the standardized  effect size  expected to be detected in the
##'  trial
##' @return \bold{designs} contains the recommended optimal design parameters for the K+M-experimental arm
##'  trial including:
##' @return \verb{    }\emph{n2} and \emph{n0_2}, the sample sizes of each of the K+M-experimental arm
##'  experimental arms and its corresponding concurrent control, respectively, in the K+M-experimental arm trial
##' @return \verb{    }\emph{nt} and \emph{n0t}, the number of patients already enrolled on each
##'  of the K initial experimental arms and the common control arm, respectively, at the time the
##'  M new arms are added
##' @return \verb{    }\emph{nc}, the total sample size of the control arm for the K+M-experimental arm trial, i.e.
##'  , the sum of concurrent (n0_2) and nonconcurrent (n0t) controls
##' @return \verb{    }\emph{N2}, the total sample size of the two-period K+M-experimental arm trial
##' @return \verb{    }\emph{A1}, the allocation ratio (control to experimental arm) before the M
##'  new experimental arms are added and after the initial K experimental arms end
##' @return \verb{    }\emph{A2},  the allocation ratio after the M new experimental arms are
##'  added and before the initial K experimental arms end
##' @return \verb{    }\emph{cor1}, the correlation of Z-test statistics between any two of the K
##'  initially opened experimental arms (or between any two of the M newly added arms)
##' @return \verb{    }\emph{cor2}, the correlation of Z-test statistics between any pair of one
##'  initially opened and one newly added experimental arm
##' @return \verb{    }\emph{critical_value2}, the critical value for the comparison between each
##'              experimental arms and the corresponding control in the K+M-experimental arm trial
##' @return \verb{    }\emph{mariginal.power2}, the marginal power for the K+M-experimental arm trial
##' @return \verb{    }\emph{disjunctive.power2}, the disjunctive power for the K+M-experimental arm trial
##' @return \verb{    }\emph{FWER2}, the family-wise error rate for the K+M-experimental arm trial.
##' @return \verb{    }\emph{delta}, the standardized clinical effect size
##'               expected to be detected in the trial
##' @return \verb{    }\emph{save}, the number of patients saved in the K+M-experimental arm trial compared to
##'  conducting one K-experimental arm and one M-experimental arm trial, separately.

##'
##' @import mvtnorm
##' @export
##'
##' @author \verb{    }Xiaomeng Yuan, Haitao Pan
##'
##' @references \verb{    }Pan, H., Yuan, X. and Ye, J. (2022). An optimal two-period multiarm
##'  platform design with new experimental arms added during the trial. Manuscript submitted for
##'  publication.
##' @references \verb{    }Dunnett, C. W. (1955). A multiple comparison procedure for comparing
##' several treatments with a control. Journal of the American Statistical
##' Association, 50(272), 1096-1121.
##'
##' @examples
##' # control fwer
##' one_design(n2 = 107, n0_2 = 198, nt = 30, K = 2, M = 2, fwer = 0.025,
##'            marginal.power = 0.8, delta = 0.4)
##' #$design_Karm
##' #  K  n1 n0_1  N1 z_alpha1 FWER1   z_beta1    Power1      cor0 delta
##' #1 2 101  143 345 2.220604 0.025 0.8416212 0.9222971 0.4142136   0.4
##' #
##' #$designs
##' #   n2 n0_2 nt n0t  nc  N2       A1       A2      cor1      cor2 critical_value2
##' #1 107  198 30  43 241 669 1.414214 2.012987 0.3508197 0.2746316        2.475233
##' #
##' #  marginal.power2  disjunctive.power2 FWER2 delta save
##' #1    0.80011                0.9853799 0.025   0.4   21
##' #
##' # control pwer
##' one_design(n2 = 76, n0_2 = 140, nt = 30, K = 2, M = 2, pwer = 0.025,
##'              marginal.power = 0.8, delta = 0.4)
##' #$design_Karm
##' #  K n1 n0_1  N1 z_alpha1      FWER1   z_beta1    Power1      cor0 delta
##' #1 2 84  119 287 1.959964 0.04647892 0.8416212 0.9222971 0.4142136   0.4
##' #
##' #$designs
##' #  n2 n0_2 nt n0t  nc  N2       A1       A2      cor1      cor2 critical_value2
##' #1 76  140 30  43 183 487 1.414214 2.108696 0.3518519 0.2437831        1.959964
##' #  marginal.power2  disjunctive.power2      FWER2 delta save
##' #1       0.8001424           0.9867451 0.08807302   0.4   87



one_design <- function(n2, n0_2, nt, K, M, fwer=NULL, pwer=NULL, marginal.power, delta, seed=123) {
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

    # original design parameters, control for fwer or pwer
    if(is.null(fwer)==0){multi <- one_stage_multiarm(K=K, fwer=fwer, marginal.power=marginal.power,
                                                     delta=delta, seed=seed)}

    else {multi <- one_stage_multiarm(K=K, pwer=pwer, marginal.power=marginal.power,
                                      delta=delta, seed=seed)}

    n1 <- multi$n1
    if (nt >= n1)
        stop("nt has to be smaller than n1")
    n0_1 <- multi$n0_1
    N1 <- multi$N1
    corMat1 <- multi$corMat1
    z_alpha1 <- multi$z_alpha1
    z_beta1 <- multi$z_beta1
    Power1 <- multi$Power1

    design_Karm <- data.frame(K=multi$K,
                              n1=multi$n1, n0_1=multi$n0_1,N1=multi$N1,
                              z_alpha1=multi$z_alpha1,FWER1=multi$FWER1,
                              z_beta1=multi$z_beta1, Power1=multi$Power1,
                              cor0=1/(A1 + 1),
                              delta=multi$delta)

    # update design parameters
    c2 <- cor.mat(K = K, M = M, n = n2, n0 = n0_2, n0t = n0t)
    corMat2 <- c2$cormat

    if(is.null(fwer)==0){
            z_alpha2 <- fwer_critical(ntrt = ntrt, fwer, corMat2, seed)$critical_val}
     else{
            z_alpha2 <- z_alpha1
         }

    a <- sqrt((1/n1 + 1/n0_1)/(1/n2 + 1/n0_2))
    z_beta2 <- a * (z_alpha1 + z_beta1) - z_alpha2
    marginal.power2 = pnorm(z_beta2)

    set.seed(seed)
    Power2 <- 1 - pmvnorm(lower = rep(-Inf, ntrt), upper = rep(-z_beta2, ntrt), mean = rep(0, ntrt), corr = corMat2)[1]

    ## if control pwer, calculate fwer for optimal designs K+M-arm trial
    if(is.null(pwer)==0){
       FWER2 <- 1 - pmvnorm(lower = rep(-Inf, ntrt), upper = rep(z_alpha2, ntrt),
                                    mean = rep(0, ntrt), corr = corMat2)[1]
    }
    ## if control fwer, report fwer for optimal designs K+M-arm trial
    if(is.null(fwer)==0){
        FWER2 <- fwer
    }

    # Upper limit of total sample size S
    if (K == M) {
        S <- 2 * N1
    } else {
        if(is.null(fwer)==0){
            multi_m <- one_stage_multiarm(K=M, fwer=fwer, marginal.power=marginal.power,
                                          delta=delta, seed=seed)
        }
        else{multi_m <- one_stage_multiarm(K=M, pwer=pwer, marginal.power=marginal.power,
                                           delta=delta, seed=seed)}

        S <- multi$N1 +multi_m$N1
    }

    designs <- data.frame(
        n2 = n2, n0_2 = n0_2,
        nt = nt, n0t = n0t,
        nc = nc, N2 = N2,
        A1 = A1, A2 = A2,
        cor1 = c2$cor1, cor2 = c2$cor2,
        critical_value2 = z_alpha2,
        marginal.power2 = marginal.power2,
        disjunctive.power2 = Power2,
        FWER2=FWER2,
        delta = delta,
        save=S-N2)

    res <- list(design_Karm=design_Karm, designs=designs)

    return(res)

}


