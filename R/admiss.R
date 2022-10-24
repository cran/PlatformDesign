##' Find the admissible set of the (n2, n0_2) pairs, given n1, n0_1, nt, ntrt and S.
##'
##' Given n1, n0_1, nt, ntrt and S, using three constraints to find the admissible set of
##' the (n2, n0_2) pairs. See the vignettes for details.
##' @title Find the admissible set in a two-period
##'        K+M-experimental arm platform design with delayed arms
##' @param n1 the sample size in each of the K experimental arms in the K-experimental arm trial
##' @param n0_1 the sample size of the common control arm in the K-experimental arm trial
##' @param nt the number of patients already enrolled on each of the K initial
##'  experimental arms when the new arms are added
##' @param ntrt the number of experimental arms in the K+M-experimental arm trial, i.e, K+M
##' @param S the upper limit of the total sample size for the K+M-experimental arm trial. It
##'        usually takes the value of the sum of the sample sizes of two separate clinical
##'        trials (one with K and another with M experimental arms, each having one control
##'        arm). The total sample size of K (or M)-arm trial can be calculated
##'        using function one_stage_multiarm().
##'
##' @return a dataframe which contains all candidate values of n2 and n0_2 in its first and second
##'  column, respectively
##'
##' @export
##' @examples
##' admiss(n1=101, n0_1=143, nt=30, ntrt=4, S=690)


admiss <- function(n1, n0_1, nt, ntrt, S) {
    A1 = n0_1/n1

    n2_lb <- nt + 1
    n0_2_lb <- ceiling(A1 * nt) + 1

    # ntrt*n2 + n0_2 + ceiling(A1*nt) < S
    n2_ub <- (S - ceiling(A1 * nt) - n0_2_lb)/ntrt
    n0_2_ub <- S - ceiling(A1 * nt) - ntrt * n2_lb

    r1 <- seq(n2_lb, n2_ub)
    r2 <- seq(n0_2_lb, n0_2_ub)
    candidate <- expand.grid(r1, r2)
    candidate$cs <- candidate$Var1 * ntrt + candidate$Var2 + ceiling(A1 * nt)

    candidate$index <- ifelse(candidate$cs <= S, 1, 0)
    pair <- candidate[candidate$index == 1, c(1,2)]

    return(pair)

}
