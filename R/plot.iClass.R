# to do:
# figure out how to deal with multiple responses in formula plot

#' @export
#' @docType methods
#' @details \strong{plot} Render plots describing iData contents.
#' @rdname iClass-plot
setGeneric("plot", function(groups, vars, iData, mask = NULL) {
  # argument checks
  if (missing(x))
    stop("Must specific iData object.")
  nogroup <- ifelse(missing(groups), TRUE, FALSE)
  novar <- ifelse(missing(vars), TRUE, FALSE)
  
  if (nogroup && novars)
    plot(0)
  else {
    if (novar)
      demog <- getDemog(x, groups)
    else
      demog <- getDemog(x, groups, vars)
    
    # actual function
    if (!nogroup) {
      if (novar)
        demog <- data.frame(matrix(0, nrow = nrow(iData@index), 0))
      if (missing(mask) && !nogroup)
        mask <- antsImageClone(iData@iList[[groups]]@mask)
      
      # create mean variable representing mask
      for (i in seq_len(length(groups))) {
        gdim <- dim(iData@iList[[groups[i]]])
        clustvec <- mask[iData@iList[[groups[i]]]@mask != 0]
        clustseq <- seq_len(gdim[2])[clustvec != 0]
        gvec <- rep(NA, nrow(iData@index))
        for (j in seq_len(gdim[1]))
          gvec[iData@index[groups[i]][j]] <-  mean(iData@iList[[groups[i]]]@iMatrix[j, clustseq])
        demog <- cbind(demog, gvec)
        colnames(demog)[ncol(demog)] <- groups[i]
      }
    }
    plot(demog)
  }
}, valueClass = c(groups = "character", 
                  vars = "character",
                  iData = "iData",
                  mask = "antsImage"))


#' @export
#' @docType methods
#' @details \strong{plot} Create plot of variables against specific cluster
#'  within contrast.
#' @rdname iClass-plot
plot.iModel <- function(iModel, contrast, cluster) {
  if (length(contrast) > 1)
    stop("Contrast must be of length 1.")
  if (is.null(iModel@C[[contrast]]$clusterImage)) {
    warning("No results to plot.")
    return(0)
  } else {
    # create mask for specific cluster
    clustimg <- antsImageClone(iModel@C[[contrast]]$clusterImage)
    if (missing(cluster))
      clustimg[clustimg != 0] <- 1
    else {
      clustimg[clustimg != cluster] <- 0
      clustimg[clustimg != 0] <- 1
    }
    plot(iModel@iData, iModel@y, mask = clustimg)
  }
}



setGeneric("plot", function(formula, data, mask) {
  z <- iFormula(formula, data)

  plot(formula, demog)
  }, valueClass = c(formula = "formula",
                    data = "iData",
                    mask = "antsImage"))



