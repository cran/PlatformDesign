##' Calculate the critical value and the marginal type-I error rate given the number of experimental
##' arms, the family-wise type I error rate and the correlation matrix of the z-statistics.
##'
##' Use the number of experimental arms, the family-wise type I error rate and the correlation
##' matrix of the Z-test statistics to calculate the marginal type I error rate and the critical value.
##' @title Calculate the critical value and the marginal type-I error rate
##' @param ntrt the number of experimental arms in the trial
##' @param fwer the family-wise error rate (FWER) to be controlled, default to be the same
##'        throughout the trial
##' @param corMat the correlation matrix of the Z-test statistics
##' @param seed an integer for random number generation for numerically evaluating integration, default = 123
##'
##' @return \verb{    }\emph{pairwise_alpha} the marginal type-I error rate for the comparison between
##'         any of the experimental arm and its corresponding control
##' @return \verb{    }\emph{critical_val}, the critical value for the comparison between any of
##'         the experimental arm and the corresponding controls
##' @return Other values returned are inputs.
##'
##' @import mvtnorm
##' @export
##'
##' @author \verb{    }Xiaomeng Yuan, Haitao Pan
##'
##' @references \verb{    }Dunnett, C. W. (1955). A multiple comparison procedure for comparing
##' several treatments with a control. Journal of the American Statistical
##' Association, 50(272), 1096-1121.
##'
##' @examples
##'
##'corMat1 <- cor.mat(K=2, M = 2, n=107, n0=198, n0t = 43)$cormat
##'fwer_critical(ntrt=4, fwer=0.025, corMat=corMat1)
##'
##' #$ntrt
##' #[1] 4
##'
##' #$fwer
##' #[1] 0.025
##'
##' #$corMat
##' #      [,1]      [,2]      [,3]      [,4]
##' #[1,] 1.0000000 0.3508197 0.2746316 0.2746316
##' #[2,] 0.3508197 1.0000000 0.2746316 0.2746316
##' #[3,] 0.2746316 0.2746316 1.0000000 0.3508197
##' #[4,] 0.2746316 0.2746316 0.3508197 1.0000000
##'
##' #$pairwise_alpha
##' #[1] 0.006657461
##'
##' #$critical_val
##' #[1] 2.475233
##'

fwer_critical <- function(ntrt, fwer, corMat, seed=123) {

    if (ntrt == 1) {
        alpha_pair <- fwer
        c <- qnorm(1 - fwer)

        res <- list(ntrt = ntrt, fwer = fwer, corMat = corMat, pairwise_alpha = alpha_pair, critical_val = c)

    } else {
        # function to solve for finding critical value for given fwer
        fwer_dunnett <- function(c, ntrt, fwer, corr) {
            int <- 1 - pmvnorm(lower = rep(-Inf, ntrt), upper = rep(c, ntrt), mean = rep(0, ntrt), corr = corr)

            return(as.double(int) - fwer)
        }

        find_critical <- function(ntrt, fwer) {
            critical = uniroot(fwer_dunnett, interval = c(0, qnorm(1 - fwer/ntrt) + 0.01), ntrt = ntrt, fwer = fwer,
                corr = corMat)
            return(critical)
        }

        set.seed(seed)
        c <- find_critical(ntrt, fwer)$root

        alpha_pair <- 1 - pnorm(c)

        res <- list(ntrt = ntrt, fwer = fwer, corMat = corMat, pairwise_alpha = alpha_pair, critical_val = c)
    }





    return(res)


}


