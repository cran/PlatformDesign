##' Find optimal design(s), provide the design parameters for a two-period K+M experimental
##' arm platform trial
##'
##' Providing an optimized design in terms of minimizing the total sample size for adding M additional
##' experimental arms in the middle of a clinical trial which originally (in the first period) have K
##' experimental arms and 1 control arm, given user-defined FWER and marginal power.
##'
##' @title Design of an optimal two-period multi-arm platform trial with delayed arms
##' @param nt the number of patients already enrolled on each of the K experimental arms in
##'        the first period at the time the new arms are added
##' @param K the number of experimental arms in the  first period in a two-period K+M trial
##' @param M the number of delayed (newly added) experimental arms added in the second period
##' @param fwer the family-wise error rate (FWER) to be controlled, default to be the same throughout
##'        the trial
##' @param marginal.power the marginal power to achieve in the first period in a two-period K+M trial
##' @param min.marginal.power the user-defined lower limit of the marginal power in the K+M
##'        trial (with K+M experimental arms and 1 common control arm), default to be the marginal
##'        power in the first period
##' @param delta the standardized clinical effect size expected to be detected in the trial
##' @param seed an integer for random number generation for numerically evaluating integration, default = 123
##'
##' @return The function returns a list, including \bold{designs}, \bold{flag.dp}, \bold{flag.mp}, and \bold{flag.dpmp}.
##' @return \bold{designs} contains the recommended optimal design parameters for periods 1 and 2
##'          including:
##' @return \verb{    }\emph{n1} and \emph{n0_1}, the sample sizes of each of the K experimental arms and the
##'              concurrent control, respectively, in the first period
##' @return \verb{    }\emph{n2} and \emph{n0_2}, the updated sample sizes of each of the K + M experimental
##'         arms and its corresponding concurrent control, respectively, after adding M experimental arms in the second period
##' @return \verb{    }\emph{nt} and \emph{n0t}, the number of patients already enrolled on each of the K
##'              experimental arms and the control arm, respectively, in the first period at the time the M new
##'              arms are added
##' @return \verb{    }\emph{nc}, the updated total sample size of the control arm
##'          after adding M experimental arms in the second period, i.e., the sum of concurrent (n0_2)
##'          and nonconcurrent (n0t) controls
##' @return \verb{    }\emph{N2}, the total sample size of the two-period K+M experimental arm (and 1 control
##'               arm) platform trial
##' @return \verb{    }\emph{A1},  the allocation ratio (control to experimental arm) before the M new
##'               experimental arms are added and after the initial K experimental arms end
##' @return \verb{    }\emph{A2},  the allocation ratio after the M new experimental arms are added and before
##'               the initial K experimental arms end
##' @return \verb{    }\emph{cor1}, the correlation of Z-test statistics between any two of the K initially opened
##'               experimental arms (or between any two of the M newly added arms)
##' @return \verb{    }\emph{cor2}, the correlation of Z-test statistics between any pair of one initially
##'          opened and one newly added experimental arm
##' @return \verb{    }\emph{critical_value1}, the critical value for the comparison between any of the K
##'              experimental arms in the first period and the corresponding control
##' @return \verb{    }\emph{critical_value2}, the critical value for the comparison between any of the K + M
##'              experimental arms in the second period and the corresponding control
##' @return \verb{    }\emph{marginal.power1} and \emph{marginal.power2}, the marginal power for the first
##'              and second period, respectively
##' @return \verb{    }\emph{disjunctive.power1} and \emph{disjunctive.power2}, the disjunctive power for the
##'               first and second period, respectively
##' @return \verb{    }\emph{standardized_effect_size}, the standardized clinical effect size
##'               expected to be detected in the trial
##'
##'@return  \bold{flag.dp}, \bold{flag.mp}, and \bold{flag.dpmp} indicate if the lower limit of
##'         disjunctive power, marginal power, and both of them has(have) met, respectively.
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
##' \donttest{platform_design(nt = 30, K = 2, M = 2, fwer = 0.025, marginal.power = 0.8, delta = 0.4)}
##'
##' #flag.dpmp == 0, lower limits of marginal and disjunctive power are both met
##' #$designs
##' #      n1  n0_1  n2 n0_2 nt n0t  nc  N2
##' #15669 101  143 107  198 30  43 241 669
##' #15994 101  143 106  202 30  43 245 669
##' #16315 101  143 105  206 30  43 249 669
##' #16632 101  143 104  210 30  43 253 669
##'
##' #        A1       A2       cor1      cor2         critical_value1 critical_value2
##' #15669 1.414214 2.012987 0.3508197 0.2746316        2.220604        2.475233
##' #15994 1.414214 2.092105 0.3441558 0.2708949        2.220604        2.475790
##' #16315 1.414214 2.173333 0.3376206 0.2671464        2.220604        2.476330
##' #16632 1.414214 2.256757 0.3312102 0.2633910        2.220604        2.476854
##'
##' #          marginal.power1 marginal.power2 disjunctive.power1 disjunctive.power2
##' #15669             0.8       0.8001100          0.9222971          0.9853799
##' #15994             0.8       0.8003363          0.9222971          0.9857541
##' #16315             0.8       0.8003878          0.9222971          0.9860900
##' #16632             0.8       0.8002699          0.9222971          0.9863903
##'
##'  #          standardized_effect_size
##' #15669          0.4
##' #15994          0.4
##' #16315          0.4
##' #16632          0.4
##'
##' #$flag.dp
##' #[1] 0
##'
##' #$flag.mp
##' #[1] 0
##'
##' #$flag.dpmp
##' #[1] 0
##'

platform_design <- function(nt, K, M, fwer, marginal.power, min.marginal.power = marginal.power, delta,
                            seed=123) {

    if (sum(c(nt, K, M)%%1) != 0)
        stop("nt, K and M should all be integers.")
    if (sum(c(nt, K, M) > 0) != length(c(nt, K, M)))
        stop("nt, K, and M should all > 0.")
    if (sum(fwer < 0, fwer > 1) == 1)
        stop("0 <= fwer <= 1 not true")
    if (sum(marginal.power < 0, marginal.power > 1) == 1)
        stop("0 <= mariginal.power <= 1 not true")
    if (sum(min.marginal.power < 0, min.marginal.power > 1) == 1)
        stop("0 <= min.mariginal.power <= 1 not true")

    flag.dp <- flag.mp <- flag.dpmp <- 0

    ntrt <- K + M
    A1 <- sqrt(K)
    n0t <- ceiling(A1 * nt)

    # original design parameters
    multi <- one_stage_multiarm(K, fwer, marginal.power, delta, seed)
    n1 <- multi$n1
    n0_1 <- multi$n0_1
    if (nt > n1)
        stop("nt has to be smaller than n1")
    if (n0t > n0_1)
        stop("n0t has to be smaller than n0_1")

    N1 <- multi$N1
    corMat1 <- multi$corMat1
    z_beta1 <- multi$z_beta1
    min.disjunctive.power <- multi$Power1
    z_alpha1 <- multi$z_alpha1

    # Upper limit of total sample size S
    if (K == M) {
        S <- 2 * N1
    } else {
        S <- ceiling((z_alpha1 + z_beta1)^2/delta^2 * (2 + 2 * sqrt(K) + 2 * sqrt(M) + K + M))
    }

    # admissible set
    pair <- admiss(n1, n0_1, nt, ntrt, S)

    # update design parameters
    Vone_design2 <- Vectorize(one_design2, vectorize.args = c("n2", "n0_2"))
    dat <- Vone_design2(K = K, M = M, n1 = n1, n0_1 = n0_1, n2 = pair$Var1, n0_2 = pair$Var2, nt = nt, n0t = n0t,
        A1, fwer = fwer, Power1 = multi$Power1, z_alpha1 = z_alpha1, z_beta1 = z_beta1, marginal.power = marginal.power,
        delta = delta, seed=seed)

    dat.u <- unlist(dat)
    dat.d <- matrix(dat.u, nrow = 19, byrow = F)
    dat.f <- data.frame(n1 = dat.d[1, ], n0_1 = dat.d[2, ], n2 = dat.d[3, ], n0_2 = dat.d[4, ], nt = dat.d[5, ], n0t = dat.d[6,
        ], nc = dat.d[7, ], N2 = dat.d[8, ], A1 = dat.d[9, ], A2 = dat.d[10, ], cor1 = dat.d[11, ], cor2 = dat.d[12,
        ], critical_value1 = dat.d[13, ], critical_value2 = dat.d[14, ], marginal.power1 = dat.d[15, ], marginal.power2 = dat.d[16,
        ], disjunctive.power1 = dat.d[17, ], disjunctive.power2 = dat.d[18, ], standardized_effect_size = dat.d[19,
        ])

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

    res <- list(designs = designs, flag.dp = flag.dp, flag.mp = flag.mp, flag.dpmp = flag.dpmp)

    return(res)

}


