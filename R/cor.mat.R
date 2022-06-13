##' Calculate the correlation matrix of the z-statistics in the two-period multi-arm platform
##' design with delayed arms, given K, M, n, n0 and n0t.
##'
##' Given K, M, n, n0 and n0t, calculate the correlation matrix of the z-statistics in the
##' two-period K+M experimental arm trial with one common control arm.
##'
##' @title Calculate the correlation matrix of the z-statistics for in a two-period
##'        multi-arm platform design with delayed arms
##' @param K the number of experimental arms in the  first period in a two-period K+M trial
##' @param M the number of delayed additional experimental arms added in the second period,
##'        default = 0 for calculating the correlation matrix of the Z-test statistics for a
##'        K-experimental arm trial in the first period
##' @param n a positive integer, which is the sample size in each of the experimental arms
##' @param n0 a positive integer, which is the sample size of the control for each of the
##'        experimental arms
##' @param n0t the number of patients already enrolled in the control arm when delayed
##'        experimental arms are added, default to NULL for calculating correlation matrix of
##'        the k-experimental arm trial in the first period
##'
##' @return \emph{cormat}, the correlation matrix of Z-test statistics in the
##'         two-period K+M experimental arm trial with one common control arm, or that in the
##'         k-experimental arm trial in the first period when M = 0
##'
##' @export
##' @examples
##' cor.mat(K = 2, M = 0, n = 101, n0 = 143)
##'
##' #$cormat
##' #        [,1]      [,2]
##' #[1,] 1.0000000 0.4139344
##' #[2,] 0.4139344 1.0000000
##'
##' #$cor1
##' #[1] 0.4139344
##'
##' #$cor2
##' #NULL
##'
##'  cor.mat(K = 2, M = 2, n = 107, n0 = 198, n0t = 43)
##'
##' #$cormat
##' #      [,1]      [,2]      [,3]      [,4]
##' #[1,] 1.0000000 0.3508197 0.2746316 0.2746316
##' #[2,] 0.3508197 1.0000000 0.2746316 0.2746316
##' #[3,] 0.2746316 0.2746316 1.0000000 0.3508197
##' #[4,] 0.2746316 0.2746316 0.3508197 1.0000000
##'
##' #$cor1
##' #[1] 0.3508197
##'
##' #$cor2
##' #[1] 0.2746316
##'

cor.mat <- function(K, M = 0, n, n0, n0t = NULL) {
    ntrt <- K + M
    old <- seq(1, K, 1)

    if (M == 0) {
        new <- NULL
    } else {
        new <- seq(K + 1, K + M, 1)
    }

    cor1 <- 1/(n0/n + 1)
    if (is.null(n0t)) {
        cor2 <- NULL
    } else {
        cor2 <- (n0 - n0t)/(n0^2/n + n0)
    }

    cormat <- matrix(0, nrow = ntrt, ncol = ntrt)

    for (i in 1:ntrt) {
        for (j in 1:ntrt) {
            if (i != j) {
                if (sum(c(i, j) %in% old) == 2 | sum(c(i, j) %in% new) == 2) {
                  cormat[i, j] <- cor1
                } else {
                  cormat[i, j] <- cor2
                }
            } else cormat[i, j] <- 1
        }
    }

    res <- list(cormat = cormat, cor1 = cor1, cor2 = cor2)

    return(res)  # return the cor matrix
}



