


#' 
#'
#'
#' quickly normalize dataset of warped images to more easily visualize individual differences from group (not for statistical analysis)
#'
#' mask <-
#' paths <- paste()
#' ilist <- imageFileNames2ImageList(paths)
#' ucsd 105
#' @export scaleSubjectWarp

scaleSubjectWarp <- function(ilist, mask) {
    imat <- imageListToMatrix(ilist, mask)
    mimat <- colMeans(imat)
    r <- imat - matrix(rep(mimat, nrow(imat)), nrow = nrow(imat), ncol = ncol(imat), byrow = TRUE)
    matrixToImages(r, mask)
}
