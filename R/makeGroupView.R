#'
#' @param plist list of paths to image files
#' @param mask
#' @param statdir
#' @param verbose
#'
#' @return
#' 
#'
#'
#'
#' @description
#' 
#' Some image modalities do not lend themselves to visualizing individual 
#' subjects. For example, the jacobian image from a VBM analysis does not 
#' account for the innevitable deformation that will occur between all subjects
#' to a template. Typically this is not an issue because most statistical 
#' methods will account for consistant deformation patterns throughout groups.
#' However, this may wash out interesting findings when visualizing single
#' subjects. 
#' 
#' Using an ANCOVA and contrasting each subject with the entire group we can 
#' easily account for these differences and create presentable subject images.
#' (Note: this type of analysis would not typically be appropriate for testing
#' a hypothesis.) 
#' 
#' The process is as follows:
#' 
#' The provided list of file paths (\code{plist}) are used to extract images 
#' which are converted into an image matrix (subjects x voxels). This image 
#' matrix is then treated as a response variable while an identity matrix of
#' the subjects is used as a design matrix (this allows each subject to be
#' treated as an individual group). Each subject is contrasted with the entire
#' group in order to account for normal differences betweeen images. 
#'
#' If \code{statdir = TRUE} then the directory from which each image was
#' obtained becomes the home for the adjusted subject image and quick view 
#' image. Alternatively, one may provide a list of desired paths/names to 
#' prefix the each image in the corresponding \code{plist}. 
#'
#'
#'
#' @export makeGroupView

makeGroupView <- function(plist, mask, statdir = NULL, verbose = TRUE) {
    if (statdir == "TRUE")
        statdir <- lapply(plist, function(x)(paste(dirname(x), "/", sep = "")))
    ilist <- imageFileNames2ImageList(plist)
    if (missing(mask))
        mask <- getMask(ilist[[1]])
    y <- imageListToMatrix(ilist, maks)
    n <- ncol(y)
    x <- diag(n)
    z <- .lm.fit(x, y)
    mrss <- sqrt(z$residuals / (n - 1))
    statList <- list()
    if (verbose)
        progress <- txtProgressBar(min = 0, max = n, style = 3)
    for (i in 1:n) {
        conmat <- matrix(rep(-1, n), nrow = 1)
        conmat[1, i] <- n - 1

        se <- t(as.matrix( * (mrss * conmat[i,] %*% XX %*% conmat[i,])))
        se[se == 0] <- 1 # instead of letting 0 / 0 = NULL this leads to 0 / 1 = 0
        statmat <- (conmat[i,] %*% z$coefficients) / se

        statList <- lappend(statList, makeImage(mask, statmat))

        if (!is.null(statdir)) {
            brain <- renderSurfaceFunction(surfimg = list(bm), alphasurf = 0.1, smoothsval = 1.5)
            make3ViewPNG(diag(4), diag(4), diag(4), paste(newpath, "QuickView.png", sep = ""))
        }
        make3ViewPNG( diag(4), diag(4), diag(4), tempfile(fileext = '.png') )
        if (verbose)
            setTxtProgressBar(progress, i)
    }
    if (verbose)
        close(progress)

}
