#' RFT p-values
#'
#' Calculates p-values of a statistical field using random field theory
#' 
#' @param D Image dimensions.
#' @param c Number of clusters.
#' @param k Spatial extent in resels
#' @param u Statistical threshold.
#' @param n Number of statistical field in conjunction.
#' @param resels Resel measurements of the search region.
#' @param df Degrees of freedom expressed as c(degrees of interest, degrees of error).
#' @param fieldType:
#' \itemize{
#' \item{"T"}{T-field} 
#' \item{"F"}{F-field} 
#' \item{"X"}{Chi-square field"} 
#' \item{"Z"}{Gaussian field}
#' }
#' 
#' @return The probability of obtaining the specified cluster
#' \itemize{
#' {"Pcor"}{"corrected p-value"}
#' {"Pu"}{"uncorrected p-value"}
#' {"Ec"}{"expected number of clusters"}
#' {"ek"}{"expected number of resels per cluster"}
#' }
#' 
#' @details 
#' 
#' This function calculates p-values of a thresholded statistical field at various levels:
#' 
#' set-level
#' rft.pval(D, c, k, u, n, resels, df, fieldType)
#' 
#' cluster-level
#' rft.pval(D, 1, k, u, n, resels, df, fieldType)
#' 
#' peak-level
#' rft.pval(D, 1, 0, u, n, resels, df, fieldType)
#' 
#' Where set-level refers to obtaining the set of clusters, cluster-level refers to a specific 
#' cluster, and peak-level refers to the maximum (or peak) of a single cluster.
#' 
#' @references 
#' Friston K.J., (1994) Assessing the Significance of Focal Activations Using Their Spatial Extent.
#' Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power.
#' Worlsey K.J., (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' 
#' @Author Zachary P. Christensen
#' 
#' @seealso rftResults, resels
#' 
#' @note: function currently in beta phase
#' @examples
#' 
#' # generate some data as if we just fitted a linear regression
#' outimg1 <- makeImage(c(10, 10, 10), rt(1000))
#' maskimg <- getMask(outimg1)
#' 
#' # create clusters using arbitrary threshold
#' clusters <- image2ClusterImages(outimg1, minClusterSize=1, minThresh = 2, maxThresh = Inf)
#' fwhm <- estSmooth(outimg1, maskimg)
#' resels <- resels(mask, fwhm$fwhm)
#' peak <- max(clusters[[1]])
#' peakP <- rftPval(3, 1, 0, 2, 1, resels, c(1, 1), fieldType="T")
#' 
#' 
#' @export rftPval
rftPval <- function(D, c, k, u, n, resels, df = c(idf, rdf), fieldType) {
    if (missing(fieldType)) {
        stop("Must specify fieldType")
    } else if (missing(df)) {
        stop("Must specify df")
    } else if (missing(resels)) {
        stop("Must specify resels")
    } else if (missing(u) && missing(k) && missing(c)) {
        stop("Must atleast specify one of u, k, or c")
    }
    G <- sqrt(pi) / gamma((1:(D + 1) / 2))
    ec <- euler(u, df, fieldType)
    ec <- pmax(ec[1:(D + 1)], .Machine$double.eps)
    P <- toeplitz(as.numeric(ec * G))
    P[lower.tri(P)] <- 0
    if (n != round(n)) {
        n <- round(n)
        warning("rounding exponent `n' to", n)
    }
    phi <- diag(nrow = nrow(P))
    pot <- P
    while (n > 0) {
        if (n %% 2)
            phi <- phi %*% pot
        n <- n %/% 2
        pot <- pot %*% pot
    }
    P <- phi
    P <- P[1,]
    EM <- (resels[1:(D + 1)] / G) * P # maxima in all dimensions
    Ec <- sum(EM) # number of overall expected maxima/clusters
    EN <- P[1] * resels[D + 1] # number of resels in entire image
    ek <- EN / EM[D + 1] # expected number of resels per cluster

    rfB <- (gamma(D / 2 + 1) / ek) ^ (2 / D)
    Punc <- exp( - rfB * (k ^ (2 / D))) # cumulative cluster-size distribution from which uncorrected P values are calculated

    Pcor <- 1 - ppois(c - 1, lambda = (Ec + .Machine$double.eps) * Punc)
    z <- list(Pcor = Pcor, Punc = Punc, Ec = Ec, ek = ek)
    z
}