##' This function can be used to design a K-experimental arm trial (with K experimental arm plus a common
##' control arm) given a pre-planned family-wise error rate (or pair-wise error rate) and with a
##' user-specified marginal power. It calculates required sample sizes for each of the
##' experimental arm (n1), the control arm (n0_1), the total sample size (N1), and the critical
##' value (z_alpha1) for each experimental arm-control comparison in the trial.
##'
##' Given the number of experimental arms (K), the family-wise type I error rate (or the pair-wise
##' type-I error-rate), the marginal power for each experimental-control comparison and the
##' standardized effect size, to calculate the sample sizes and other design parameters for the
##' K-experimental arm trial (with K-experimental arm in addition to one control arm).
##'
##' @title Calculate the sample sizes and other design parameters for an one-stage K-experimental arm trial
##'  using the root-K rule for the allocation ratio, controlling for FWER or PWER
##'
##' @param K the number of experimental arms
##' @param fwer the family-wise type I error rate, default to be null, users need to choose
##'  between controlling for fwer or pwer and input a value for this argument if choosing fwer
##' @param pwer the pair-wise type I error rate, default to be null, users need to input a value
##'  for this argument if controlling for pwer
##' @param marginal.power the marginal power for each experimental-control comparison
##' @param delta the standardized effect size expected to be detected
##'        in the trial
##' @param seed an integer used in random number generation for numerically evaluating integration,
##'  default = 123
##'
##' @return \emph{K} the number of experimental arms in the K-experimental arm trial (with K experimental arm
##'  plus a common control arm), e.g., for a 3-arm trial with 3 experimental arm and 1 control arm,
##'  K=3.
##' @return \emph{n1} the sample size for each of the K experimental arms
##' @return \emph{n0_1} the sample size of the common control arm
##' @return \emph{N1} the total sample size of a K-experimental arm trial
##' @return \emph{z_alpha1} the critical value for the comparison between any of the
##'  K-experimental arm and its corresponding control
##' @return \emph{FWER1} the family-wise type-I error rate
##' @return \emph{z_beta1} the quantile of the marginal power, i.e., qnorm(marginal power)
##' @return \emph{Power1} the disjunctive power of the K-experimental arm trial defined as the
##'              probability of rejecting at least one of the K experimental arms under the
##'              alternative hypothesis
##' @return \emph{corMat1} the correlation matrix of the Z-test statistics
##' @return \emph{delta} the standardized effect size expected to be detected in the K-experimental arm trial
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
##'  several treatments with a control. Journal of the American Statistical
##'  Association, 50(272), 1096-1121.
##'
##' @examples
##' # controlling for FWER
##' one_stage_multiarm(K = 2, fwer = 0.025, marginal.power = 0.8, delta = 0.4)
##' #$K
##' #[1] 2
##' #
##' #$n1
##' #[1] 101
##' #
##' #$n0_1
##' #[1] 143
##' #
##' #$N1
##' #[1] 345
##' #
##' #$z_alpha1
##' #[1] 2.220604
##' #
##' #$FWER1
##' #[1] 0.025
##' #
##' #$z_beta1
##' #[1] 0.8416212
##' #
##' #$Power1
##' #[1] 0.9222971
##' #
##' #$corMat1
##' #[,1]      [,2]
##' #[1,] 1.0000000 0.4142136
##' #[2,] 0.4142136 1.0000000
##' #
##' #$delta
##' #[1] 0.4
##' #
##' # controlling for pwer
##' one_stage_multiarm(K = 2, pwer = 0.025, marginal.power = 0.8, delta = 0.4)
##' #$K
##' #[1] 2
##' #
##' #$n1
##' #[1] 84
##' #
##' #$n0_1
##' #[1] 119
##' #
##' #$N1
##' #[1] 287
##' #
##' #$z_alpha1
##' #[1] 1.959964
##' #
##' #$FWER1
##' #[1] 0.04647892
##' #
##' #$z_beta1
##' #[1] 0.8416212
##' #
##' #$Power1
##' #[1] 0.9222971
##' #
##' #$corMat1
##' #[,1]      [,2]
##' #[1,] 1.0000000 0.4142136
##' #[2,] 0.4142136 1.0000000
##' #
##' #$delta
##' #[1] 0.4

one_stage_multiarm <- function(K, fwer=NULL, pwer=NULL, marginal.power, delta, seed=123) {

    if (sum(is.null(fwer)+is.null(pwer))!=1)
        stop("Users need to decide if controlling for fwer or pwer.
              Please input a value for either fwer or pwer.")

    # control fwer
    if(is.null(fwer)==0){
        if (sum(fwer <= 0, fwer >= 1) == 1)
            stop("0 < fwer < 1 not true")}

    # control pwer
    if(is.null(pwer)==0){
        if (sum(pwer <= 0, pwer >= 1) == 1)
            stop("0 < pwer < 1 not true")}

    if (sum(marginal.power <= 0, marginal.power >= 1) == 1)
        stop("0 < mariginal.power < 1 not true")

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

    # control for fwer or pwer
    if(is.null(fwer)==0){
            fwer.old <- fwer_critical(K, fwer, corMat = corMat1, seed)
            z_alpha1 <- fwer.old$critical_val  # z_(1-alpha1)
            FWER1=fwer
        }
     else{  z_alpha1 <- qnorm(1-pwer)
            if (K==1){FWER1=pwer}
             else{
                 FWER1 <- 1 - pmvnorm(lower = rep(-Inf, K), upper = rep(z_alpha1, K),
                                      mean = rep(0, K), corr = corMat1)[1]}
          }


    if (K == 1) {
        Power1 = marginal.power
    } else {
        set.seed(seed)
        Power1 <- 1 - pmvnorm(lower = rep(-Inf, K), upper = rep(-z_beta1, K),
                              mean = rep(0, K), corr = corMat1)[1]
    }

    n1 <- ceiling((z_alpha1 + z_beta1)^2/delta^2 * (1 + sqrt(K)/K))
    n0_1 <- ceiling(n1 * sqrt(K))
    N1 <- K * n1 + n0_1

    res <- list(K = K,
                n1 = n1, n0_1 = n0_1, N1 = N1,
                z_alpha1 = z_alpha1, FWER1=FWER1,
                z_beta1 = z_beta1, Power1 = Power1,
                corMat1 = corMat1, delta=delta)

    return(res)
}


