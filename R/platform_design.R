##' Find optimal design(s) for a two-period K+M experimental arm platform trial given a
##' user-specified family-wise error rate (or pair-wise error rate) and marginal power. The
##' K+M-experimental arm trial has K experimental arms and one control arm during the first period, and later M
##' experimental arms are added on the start of the second period. The one common control arm is
##' shared among all experimental arms across the trial. The function calculates required sample
##' sizes for each of the experimental arm (n2), the concurrent control (n0_2), the total sample
##' size (N2), the allocation ratios (A1 & A2), and the critical value (z_alpha1) for each
##' experimental arm-control comparison in the trial. The number of patients saved in a K+M-experimental arm
##' trial compared to conducting one K-experimental arm and one M-experimental arm trial separately is also provided. Users
##' can choose to control for either FWER or PWER in the trial.
##'
##' Providing an optimized design in terms of minimizing the total sample size for adding M
##' additional experimental arms in the middle of a clinical trial which originally
##' have K experimental arms and 1 control arm, given user-defined FWER (or PWER) and marginal
##' power. The optimal design for the K+M-experimental arm trial exists only if flag.dpmp = 0. It means that
##' the  optimal design can be found to keep both marginal and disjunctive power levels no less
##' than those in the corresponding K-experimental arm trial. If flag.dpmp = 1 and flag.mp = 1, it means the
##' optimal design to maintain both mariginal and disjunctive power levels can not be found, but
##' the a design with the disjunctive power no less than its counterpart in the K-experimental arm trial is
##' returned in \bold{designs}.
##'
##' @title Design an optimal two-period multiarm platform trial with new experimental arms added
##'  during the trial, controlling for FWER or PWER
##'
##' @param nt the number of patients already enrolled on each of the K initial experimental arms
##'  at the time the M new arms are added.
##' @param K the number of experimental arms in the first period in a two-period K+M-experimental arm trial
##' @param M the number of new experimental arms added at the start of the second period
##' @param fwer the family-wise type I error rate, default to be null, users need to choose
##'  between controlling for fwer or pwer and input a value for this argument if fwer is chosen
##' @param pwer the pair-wise type I error rate, default to be null, users need to choose between
##'        controlling for fwer or pwer and input a value for this argument if pwer is chosen
##' @param marginal.power the marginal power for each experimental-control comparison in the
##'        K-experimental arm trial. This is also the marginal power the algorithm aims to achieve
##'        in the K+M-experimental arm when min.marginal.power=marginal.power (default option).
##' @param min.marginal.power the marginal power the function aims to achieve in the K+M-experimental
##'        arm trial, default to be the same as the marginal power of the K-experimental arm trial.
##'        It will be the marginal power of the K+M-experimental arm if optimal design exists.
##'        Don't change the default unless you need to achieve a marginal power level different than
##'        that of the K-experimental arm trial.
##' @param delta the standardized effect size expected to be detected in the trial
##' @param seed an integer used in random number generation for numerically evaluating
##'  integration, default = 123
##'
##' @return The function returns a list, including \bold{design_Karm}, \bold{designs},
##'  \bold{flag.dp}, \bold{flag.mp}, and \bold{flag.dpmp}.
##'
##' @return \bold{design_Karm} contains the design parameters for the K-experimental arm trial including:
##' @return \verb{    }\emph{K}, the number of experimental arms
##' @return \verb{    }\emph{n1}, the sample size for each of the K experimental
##'  arms
##' @return \verb{    }\emph{n0_1}, the sample size of the common control arm
##' @return \verb{    }\emph{N1} the total sample size of a K-experimental arm trial
##' @return \verb{    }\emph{z_alpha1}, the critical value for the comparison between any of the K
##'  experimental arms and the control
##' @return \verb{    }\emph{FWER1}, the family-wise error rate
##' @return \verb{    }\emph{z_beta1}, the quantile of the marginal power, i.e., qnorm(marginal
##'  power)
##' @return \verb{    }\emph{Power1}, the disjunctive power
##' @return \verb{    }\emph{cor0}, the correlation of Z-test statistics between any two of the K
##'  experimental arms
##' @return \verb{    }\emph{delta}, the standardized  effect size expected to be detected in the
##'  K-experimental arm trial
##'
##' @return \bold{designs} contains the recommended optimal design parameters for the K+M-experimental arm
##'  trial including:
##' @return \verb{    }\emph{n2} and \emph{n0_2}, the sample sizes of each of the K+M experimental
##'         arms and its corresponding concurrent control, respectively
##' @return \verb{    }\emph{nt} and \emph{n0t}, the number of patients already enrolled on each
##'  of the K initial experimental arms and the control arm, respectively, at the time the M new
##'  arms are added
##' @return \verb{    }\emph{nc}, the total sample size of the control arm for the k+M trial, i.e.
##'  , the sum of the concurrent (n0_2) and nonconcurrent (n0t) controls
##' @return \verb{    }\emph{N2}, the total sample size of the two-period K+M-experimental arm trial
##' @return \verb{    }\emph{A1}, the allocation ratio (control to experimental arm) before the
##'  M new experimental arms are added and after the initial K experimental arms end
##' @return \verb{    }\emph{A2}, the allocation ratio (control to experimental arm) after the M
##'  new experimental arms are added and before the initial K experimental arms end
##' @return \verb{    }\emph{cor1}, the correlation of Z-test statistics between any two of the K
##'  initial experimental arms (or between any two of the M new arms)
##' @return \verb{    }\emph{cor2}, the correlation of Z-test statistics between any pair of one
##'  initially opened and one newly added experimental arm
##' @return \verb{    }\emph{critical_value2}, the critical value for the comparison between each
##'              experimental arm and the concurrent control in the K+M-experimental arm trial
##' @return \verb{    }\emph{mariginal.power2}, the marginal power for the K+M-experimental arm trial
##' @return \verb{    }\emph{disjunctive.power2}, the disjunctive power for the K+M-experimental arm trial
##' @return \verb{    }\emph{FWER2}, the family-wise type-I error rate for the K+M-experimental arm trial
##' @return \verb{    }\emph{delta}, the standardized effect size expected to be detected in the
##'  K+M-experimental arm trial
##' @return \verb{    }\emph{save}, the number of patients saved in the K+M-experimental arm trial compared to
##'  conducting one K-experimental arm and one M-experimental arm trial separately.
##'
##'@return  \bold{flag.dp}, \bold{flag.mp}, and \bold{flag.dpmp} indicate if the lower limit of
##'         disjunctive power, marginal power, or both of them has(have) met, respectively
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
##'  several treatments with a control. Journal of the American Statistical Association, 50(272),
##'  1096-1121.
##'
##' @examples
##' \donttest{platform_design(nt = 30, K = 2, M = 2, fwer = 0.025, marginal.power = 0.8,
##'  delta = 0.4)}
##' #flag.dpmp == 0, lower limits of marginal and disjunctive power are both met
##' #
##' #$design_Karm
##' #   K  n1 n0_1  N1  z_alpha1 FWER1   z_beta1      Power1      cor0   delta
##' # 1 2 101  143 345  2.220604 0.025 0.8416212   0.9222971 0.4142136   0.4
##' #
##' #$designs
##' #       n2 n0_2 nt n0t  nc  N2
##' #15669 107  198 30  43 241 669
##' #15994 106  202 30  43 245 669
##' #16315 105  206 30  43 249 669
##' #16632 104  210 30  43 253 669
##' #
##' #        A1       A2       cor1      cor2          critical_value2
##' #15669 1.414214 2.012987 0.3508197 0.2746316       2.475233
##' #15994 1.414214 2.092105 0.3441558 0.2708949       2.475790
##' #16315 1.414214 2.173333 0.3376206 0.2671464       2.476330
##' #16632 1.414214 2.256757 0.3312102 0.2633910       2.476854
##' #
##' #      marginal.power2 disjunctive.power2
##' #15669  0.8001100      0.9853799
##' #15994  0.8003363      0.9857541
##' #16315  0.8003878      0.9860900
##' #16632  0.8002699      0.9863903
##' #
##' #         FWER2    delta     save
##' #15669    0.025      0.4       21
##' #15994    0.025      0.4       21
##' #16315    0.025      0.4       21
##' #16632    0.025      0.4       21
##' #
##' #$flag.dp
##' #[1] 0
##' #
##' #$flag.mp
##' #[1] 0
##' #
##' #$flag.dpmp
##' #[1] 0

platform_design <- function(nt, K, M, fwer=NULL, pwer=NULL, marginal.power,
                            min.marginal.power = marginal.power, delta,
                            seed=123) {

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

    if (sum(c(nt, K, M)%%1) != 0)
        stop("nt, K and M should all be integers.")
    if (sum(c(nt, K, M) > 0) != length(c(nt, K, M)))
        stop("nt, K, and M should all > 0.")

    if (sum(marginal.power <= 0, marginal.power >= 1) == 1)
        stop("0 < mariginal.power < 1 not true")
    if (sum(min.marginal.power <= 0, min.marginal.power >= 1) == 1)
        stop("0 < min.mariginal.power < 1 not true")

    flag.dp <- flag.mp <- flag.dpmp <- 0

    ntrt <- K + M
    A1 <- sqrt(K)
    n0t <- ceiling(A1 * nt)

    # original design parameters, control for fwer or pwer
    if(is.null(fwer)==0){multi <- one_stage_multiarm(K=K, fwer=fwer, marginal.power=marginal.power,
                                                     delta=delta, seed=seed)}

     else {multi <- one_stage_multiarm(K=K, pwer=pwer, marginal.power=marginal.power,
                                       delta=delta, seed=seed)}

    n1 <- multi$n1
    n0_1 <- multi$n0_1
    if (nt >= n1)
        stop("nt has to be smaller than n1")
    if (n0t >= n0_1)
        stop("n0t has to be smaller than n0_1")

    N1 <- multi$N1
    corMat1 <- multi$corMat1
    z_beta1 <- multi$z_beta1
    Power1 <- multi$Power1
    min.disjunctive.power <- multi$Power1
    z_alpha1 <- multi$z_alpha1
    FWER1 <- multi$FWER1

    design_Karm <- data.frame(K=multi$K,
                              n1=multi$n1, n0_1=multi$n0_1,N1=multi$N1,
                              z_alpha1=multi$z_alpha1,FWER1=multi$FWER1,
                              z_beta1=multi$z_beta1, Power1=multi$Power1,
                              cor0=1/(A1 + 1),
                              delta=multi$delta)

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

    # admissible set
    pair <- admiss(n1, n0_1, nt, ntrt, S)

    # update design parameters
    Vone_design2 <- Vectorize(one_design2, vectorize.args = c("n2", "n0_2"))

    dat <- Vone_design2(K = K, M = M,
                        n1 = n1, n0_1 = n0_1,
                        n2 = pair$Var1, n0_2 = pair$Var2,
                        nt = nt, n0t = n0t,
                        fwer = fwer,
                        z_alpha1 = z_alpha1, z_beta1 = z_beta1,
                        seed=seed)

    dat.u <- unlist(dat)
    dat.d <- matrix(dat.u, nrow = 12, byrow = F)

    dat.f <- data.frame(n2 = dat.d[1, ], n0_2 = dat.d[2, ],
                        nt = dat.d[3, ], n0t = dat.d[4,],
                        nc = dat.d[5, ], N2 = dat.d[6, ],
                        A1 = A1, A2 = dat.d[7, ],
                        cor1 = dat.d[8, ], cor2 = dat.d[9,],
                        critical_value2 = dat.d[10, ],
                        marginal.power2 = dat.d[11,],
                        disjunctive.power2 = dat.d[12, ])

    # Select recommended designs
    dat.dp <- dat.f[dat.f$disjunctive.power2 >= min.disjunctive.power, ]
    if (nrow(dat.dp) == 0) {
        flag.dp <- 1
    }

    dat.mp <- dat.f[dat.f$marginal.power2 >= min.marginal.power, ]
    if (nrow(dat.mp) == 0) {
        flag.mp <- 1
    }

    ## designs meets both limits
    dat.dpmp <- dat.f[dat.f$disjunctive.power2 >= min.disjunctive.power & dat.f$marginal.power2 >= min.marginal.power,
        ]
    if (nrow(dat.dpmp) == 0) {
        flag.dpmp <- 1
    }

    if (flag.dpmp == 0) {
        dats <- dat.dpmp
        message("flag.dpmp == 0, lower limits of marginal and disjunctive power are both met")
    } else if (flag.dpmp == 1 & flag.dp == 0) {
        dats <- dat.dp
        warning("flag.mp == 1, lower limit of the marginal power can not be met, recommended design(s) selected from designs with disjunctive power >= Power1")
    } else if (flag.dpmp == 1 & flag.mp == 0) {
        dats <- dat.mp
        warning(" flag.dp == 1, lower limit of the disjunctive power can not be met, recommended design(s) selected from designs with marginal power >= min.marginal.power")
    } else {
        dats <- NULL
        warning("The lower limit of neither marginal nor disjuctive power is met, please redefine nt or min.marginal.power.")
    }

    ## Select designs with minimum sample size
    designs <- dats[dats$N2 == min(dats$N2), ]

    ## if control pwer, calculate fwer for optimal designs K+M trial
    if(is.null(pwer)==0){
        FWER2 <- rep(0, nrow(designs))
        for (i in 1:nrow(designs)){
          corMat2 <- cor.mat(K, M, n=designs$n2[i], n0=designs$n0_2[i], n0t=designs$n0t[i])$cormat

          FWER2[i] <- 1 - pmvnorm(lower = rep(-Inf, ntrt), upper = rep(z_alpha1, ntrt),
                                  mean = rep(0, ntrt), corr = corMat2)
        }

        designs$FWER2 <- FWER2

    }
    ## if control fwer, report fwer for optimal designs K+M trial
    if(is.null(fwer)==0){
        designs$FWER2 <- fwer
    }

    designs$delta <- delta
    designs$save <- S-designs$N2


    res <- list(design_Karm = design_Karm,
                designs = designs,
                flag.dp = flag.dp,
                flag.mp = flag.mp,
                flag.dpmp = flag.dpmp)

    return(res)

}


