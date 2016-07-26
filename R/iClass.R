# to do:
# fix iData$split() function

#' iData/iGroup reference class
#' 
#' @details 
#' iData/iGroup refers to a set of reference classes used for representing 
#' image data. The base compositional unit of is the \code{\link{iGroup}} class
#' that allows a set of images to be represented by an image matrix, name, 
#' binary mask (of class \code{\link{antsImage}}), and information describing 
#' the images' modality (i.e. fMRI, PET, VBM, etc.). 
#' 
#' @usage
#' iGroup(name, iMatrix, mask, modality, ...)
#' 
#' iData(iList, demog, index, ...)
#' 
#' @note
#' Although the use of reference classes and pointers allows more efficient
#' memory usage, it may have unfamiliar side effects to those used to typical
#' R syntax and unfamiliar with other languages such as Java and C++ (see
#' examples below and \link{refClass-class} for more details).
#' 
#' @example 
#' 
#' # create iGroup object
#' ilist <- getANTsRData("population")
#' mask <- getMask(ilist[[1]])
#' imat <- imageListToMatrix(ilist, mask)
#' iGroup1 <- iGroup(name = "pop1", iMatrix = imat, mask = mask, modality = "T1")
#' 
#' ilist <- lappend(ilist, ilist[[1]])
#' imat <- imageListToMatrix(ilist, mask)
#' iGroup2 <- iGroup(name = "pop2", iMatrix = imat, mask = mask, modality = "T1")
#' 
#' # save iGroup object
#' tmpfile <- tempfile(fileext = ".h5")
#' iGroupWrite(iGroup1, tmpfile)
#' 
#' # load saved iGroup object
#' print(iGroupRead(tmpfile, "pop1"))
#' 
#' 
#' demog <- data.frame(id = c("A", "B", "C", NA),
#'   age = c(11, 7, 18, 22), sex = c("M", "M", "F", "F"))
#'   
#' bool1 <- c(TRUE, TRUE, TRUE, FALSE)
#' bool2 <- c(TRUE, TRUE, TRUE, TRUE)
#' 
#' # create iData object that holds demographics info
#' mydata <- iData(demog = demog)
#' 
#' # add iGroup object to iData
#' mydata$addGroup(iGroup1, bool1)
#' mydata$addGroup(iGroup2, bool2)
#' 
#' # ensure only active voxels are included in mask
#' mydata$checkMask("pop1", "pop2")
#' 
#' # save iData object
#' iDataSave(mydata, tmpfile)
#' # which wraps the following alternative reference class approach
#' # mydata$save(tmpfile)
#' 
#' # load saved iData object
#' print(iDataLoad(tmpfile))
#' # which wraps the following alternative reference class approach
#' # object <- iData()
#' # object$load(tmpfile)
#' 
#' # split iData object into k-folds or train/test groups
#' split_list <- mydata$split(0.3)
#' 
#' # retreive demographic information specific to an iGroup
#' mydata$getDemog("pop1", "age", "sex")
#'
#' @export iGroup/iData-class


#' iGroup reference class
#' 
#' @param name
#' @param iMatrix
#' @param mask
#' @param modality
#' 
#' @slot name
#' @slot iMatrix
#' @slot mask
#' @slot modality
#' @slot K
#' 
#' @methods
#' initialize
#' checkMask
#' iGroupRead
#' iGroupWrite
#' print
#'
#'
#' @export iGroup
iGroup <- setClass("iGroup",
                   slots = list(
                     name = "character",
                     iMatrix = "matrix",
                     mask = "antsImage",
                     modality = "character",
                     K = "list")
                   )

iGroupMake <- function(name = "iGroup1", iMatrix, mask, modality = "fMRI", rowslist, HParam, RT) {
  if (!missing(rowslist) | !missing(HParam) | !missing(RT)) {
    if (missing(rowslist) | missing(HParam) | missing(RT))
      stop("rowslist, HParam, and RT must all be provided if any one is.")
    if ((length(rowslist) != length(HParam)) | (length(rowslist) != length(RT)))
      stop("Length of rowslist, HParam, and RT must be the same.")
    K <- list()
    for (i in seq_len(length(rowslist)))
      K[[i]] <- list(row = rowslist[[i]], HParam = HParam[i], RT = RT[i])
    K <- rftFilter(K)
  } else 
    K <- 1
  return(iGroup(name = name, iMatrix = iMatrix, mask = antsImageClone(mask),
                modality = modality, K = K))
}


#' @param filename 
#' @param name
#' @describeIn iGroup
iGroupRead <- function(filename, name) {
  if (!usePkg("h5"))
    stop( "Please install package h5 in order to use this function." )
  if (!file.exists(filename))
    stop("file does not exist")
  file <- h5file(filename)
  dnames <- list.datasets(file["iData/iList"])
  if (!any(basename(dirname(dnames)) == name))
    stop("name is note a saved iGroup within specified filename.")
  tmp <- file.path("iData/iList", name)
  name <- name
  iMatrix <- file[file.path(tmp, "iMatrix")][]
  modality <- file[file.path(tmp, "modality")][]
  tmpmask <- file.path(tmp, "mask")
  mymask <- as.antsImage(file[tmpmask][])
  k = antsSetSpacing(mymask, h5attr(file[tmpmask], "spacing"))
  k = antsSetOrigin(mymask, h5attr(file[tmpmask], "origin"))
  k = antsSetDirection(mymask, h5attr(file[tmpmask], "direction"))
  mask <- mymask
  h5close(file)
  return(iGroup(name, iMatrix, mask, modality))
}

#' @describeIn iGroup
setMethod("print", "iGroup", function(x) {
  cat("iGroup object:\n")
  cat(" ", x@name, "\n")
  cat(" images =", nrow(x@iMatrix), "\n")
  cat(" voxels =", ncol(x@iMatrix), "\n")
  cat(" modality =", x@modality, "\n")
  cat("___\n")
})

#' @describeIn iGroup
iGroupWrite <- function(object, filename) {
  if (class(object) != "iGroup")
    stop("object must be of class iGroup")
  if (!usePkg("h5"))
    stop("Please install package h5 in order to use this function.")
  if (file.exists(filename))
    stop("Stopping because filename exists.")
  file <- h5file(filename)
  tmp <- file.path("iData/iList", object@name)
  file[file.path(tmp, "iMatrix")] <- object@iMatrix
  file[file.path(tmp, "modality")] <- object@modality
  file[file.path(tmp, "mask")] <- as.array(object@mask)
  h5attr(file[file.path(tmp, "mask")] , "spacing") <- antsGetSpacing(object@mask)
  h5attr(file[file.path(tmp, "mask")] , "direction") <- antsGetDirection(object@mask)
  h5attr(file[file.path(tmp, "mask")] , "origin") <- antsGetOrigin(object@mask)
  h5close(file)
  return(TRUE)
}

## iData functions----

#' iData reference class
#' 
#' @param iList
#' @param demog
#' @param index
#' @param ...
#' 
#' @field iList
#' @field demog
#' @field index
#' 
#' @method 
#' initialize
#' addGroup
#' checkMask
#' getDemog
#' load
#' save
#' show
#' split
#' 
#' 
#' 
#' @export iData
iData <- 
  setRefClass("iData",
              fields = list(
                iList = "list",
                demog = "data.frame",
                index = "data.frame"),
              methods = list(
                initialize = function(iList, demog, index, ...) {
                  ll <- list(...)
                  if (missing(iList)) {
                    iList <<- list("iGroup1" = iGroup())
                  } else {
                    lnames <- rep("tmp", length(iList))
                    for (i in 1:length(iList))
                      lnames[i] <- iList[[i]]$name
                    iList <<- iList
                    names(iList) <<- lnames
                  }
                  
                  demog <<- if (missing(demog)) data.frame() else as.data.frame(demog)
                  
                  if (missing(index)) {
                    if (!missing(iList) && !missing(demog)) {
                      tmpbool <- rep.int(NA, nrow(demog))
                      tmpindex <- data.frame(tmpbool)
                      tmpindex[seq_len(nrow(iList[[1]]$iMatrix)), 1] <- seq_len(nrow(iList[[i]]$iMatrix))
                      if (length(iList) > 1) {
                        for (i in 2:length(iList)) {
                          tmpindex <- cbind(tmpindex, data.frame(tmpbool))
                          tmpindex[seq_len(nrow(iList[[i]]$iMatrix)), i] <- seq_len(nrow(iList[[i]]$iMatrix))
                        }
                      }
                      colnames(tmpindex) <- names(iList)
                    } else
                      tmpindex <- data.frame()
                  } else {
                    if (class(index) == "logical") {
                      tmpbool <- index
                      tmpbool[tmpbool == TRUE] <- seq_len(sum(tmpbool))
                      tmpindex <- data.frame(tmpbool)
                    } else 
                      tmpindex <- index
                    if (!missing(iList) && !missing(demog)) {
                      if (ncol(tmpindex) != length(iList))
                        stop("Not all iGroups are represented in the index.")
                      if (nrow(tmpindex) != nrow(demog))
                        stop("All rows of demog must be represented by rows in index.")
                      for (i in seq_len(ncol(tmpindex))) {
                        if (sum(tmpindex[, i]) != sum(seq_len(nrow(iList[[i]]$iMatrix))))
                          stop(paste("Number of non NA rows in index does not equal number of rows in iGroup ", names(iList)[i], ".", sep = ""))
                      }
                      colnames(tmpindex) <- names(iList)
                    } else
                      tmpindex <- index
                  }
                  index <<- tmpindex
                },
                addGroup = function(iGroup, bool) {
                  'add iGroup to iList field'
                  if (class(iGroup) != "iGroup")
                    stop("iGroup must be of class iGroup.")
                  if (length(bool) != nrow(demog))
                    stop("The number of rows in demog must be equal to the length of bool.")
                  if (sum(bool) != nrow(iGroup@iMatrix))
                    stop("Number of TRUE elements in bool is not equal to number of images in iGroup")
                  if (any(names(iList) == iGroup@name))
                    stop(paste("iGroup of name ", iGroup@name, " already exists in iList", sep = ""))
                  bool[TRUE] <- seq_len(sum(bool))
                  if (any(dim(index) == 0)) {
                    index <<- data.frame(bool)
                    colnames(index) <<- iGroup@name
                  } else {
                    index <<- cbind(index, bool)
                    colnames(index)[ncol(index)] <<- iGroup@name
                  }
                  if (length(iList) != 0)
                    iList <<- list(iGroup)
                  else
                    iList <<- lappend(iList, iGroup)
                  names(iList)[length(iList)] <<- iGroup@name
                  colnames(index) <<- names(iList)
                },
                checkMask = function(...) {
                  'ensures only active voxels are present in mask of iGroup of the given name'
                  lname <- c(...)
                  if (class(lname) != "character")
                    stop("Argument must be a character specifying an existing iGroup within iList.")
                  for (i in seq_len(length(lname))) {
                    if (!any(names(iList) == lname[i]))
                      stop("Argument must be a character specifying an existing iGroup within iList.")
                    nimg <- nrow(iList[[lname[i]]]@iMatrix)
                    mask_vec <- iList[[lname[i]]]@iMatrix
                    mask_vec[mask_vec != 0] <- 1
                    mask_vec <- colSums(mask_vec)
                    mask_vec[mask_vec != nimg] <- 0
                    mask_vec[mask_vec == nimg] <- 1
                    iMatrix <- iList[[lname[i]]]@iMatrix[, as.logical(mask_vec)]
                    mask <- makeImage(iList[[lname[i]]]@mask, mask_vec)
                    iList[[lname[i]]] <<- iGroup(name = iList[[lname[i]]]@lname[i], iMatrix = iMatrix, mask = mask, modality = iList[[lname[i]]]@modality)
                  }
                  return(TRUE)
                },
                dropGroup = function(...) {
                  'drop iGroup objects with provided names from iList'
                  tmp <- c(...)
                  iList <<- iList[[-tmp]]
                },
                load = function(filename, verbose = TRUE) {
                  'loads previously saved iData from h5 file'
                  if (!usePkg("h5"))
                    stop( "Please install package h5 in order to use this function." )
                  if (!file.exists(filename))
                    stop("file does not exist")
                  file <- h5file(filename)
                  
                  # load each iGroup
                  if (verbose)
                    cat("Loading iList...")
                  inames <- unique(basename(dirname(list.datasets(file["iData/iList"]))))
                  tmplist <- list()
                  for (i in 1:length(inames)) {
                    if (verbose)
                      cat( "\n ", inames[i])
                    tmplist[[i]] <- iGroupRead(filename, inames[i])
                  }
                  names(tmplist) <- inames
                  iList <<- tmplist
                  if (verbose)
                    cat(". \n")
                  
                  # load index
                  if (verbose)
                    cat("Loading index. \n")
                  index <<- data.frame(file["iData/index"][])
                  colnames(index) <<- h5attr(file["iData/index"], "colnames")
                  
                  # load demog
                  if (verbose)
                    cat("Loading demog. \n")
                  tmp <- file["iData/demog/matrix"][]
                  dnam <- h5attr(tmp, "colnames")
                  colnames(tmp) <- dnam
                  for (i in 1:length(dnam)) {
                    colnum <- paste("iData/demog/col", i, sep = "")
                    tmpattr <- list.attributes(file[colnum])
                    for (j in 1:length(tmpattr))
                      attr(tmp[dnam[i]], tmpattr[j]) <- h5attr(file[colnum], tmpattr[j])
                  }
                  demog <<- tmp
                  h5close(file)
                },
                getDemog = function(groups, vars) {
                  'get variables indexed according to iGroup indicated by groups'
                  if (!any(groups == groupss(iList)))
                    stop("groups does not match and iGroups in iList.")
                  if (length(groups) == 1)
                    rindex <- rowgroupss(index[groups])[!is.na(index[groups])]
                  else
                    rindex <- as.logical(rowSums(index[groups]))
                  if (missing(vars))
                    return(demog[rindex, ])
                  else
                    return(demog[rindex, vars])
                },
                getGroup = function(...) {
                  'get iGroups '
                  groups <- c(...)
                  return(iList[[groups]])
                },
                getGroupIndex = function(groups) {
                  'create index for retreiving complimentary iGroup subjects'
                  out <- index[groups]
                  for (i in seq_len(nrow(index))) {
                    if (any(is.na(out[i, ])))
                      out[i, ] <- out[i, ] * -1
                  }
                  return(out)
                },
                save = function(filename, verbose = TRUE) {
                  'saves iData object to h5 file'
                  if (!usePkg("h5"))
                    stop("Please install package h5 in order to use this function.")
                  if (file.exists(filename))
                    stop("Stopping because filename exists.")
                  file <- h5file(filename)
                  
                  # save each iGroup
                  if (verbose)
                    cat("Saving iGroup... \n")
                  for (i in 1:length(iList)) {
                    if (verbose)
                      cat(" ", names(iList)[i], "\n")
                    iGroupWrite(iList[[i]], filename)
                  }
                  
                  # save index
                  if (verbose)
                    cat("Saving index. \n")
                  file["iData/index"] <- data.matrix(index)
                  h5attr(file["iData/index"], "colnames") <- colnames(index)
                  
                  # save demog
                  if (verbose)
                    cat("Saving demog. \n")
                  file["iData/demog/matrix"] <- data.matrix(demog)
                  dnam <- colnames(demog)
                  h5attr(file["iData/demog/matrix"], "colnames") <- dnam
                  if (verbose)
                    cat("Saving demographic information. \n")
                  for (i in 1:length(dnam)) {
                    dattrfile <- file[paste("iData/demog/col", i, sep = "")]
                    dattr <- names(attributes(demog[dnam[i]]))
                    if (!is.null(dattr)) {
                      for (j in 1:length(dattr))
                        h5attr(dattrfile, dattr[j]) <- attr(demog[dnam[i]], dattr[j])
                    }
                  }
                  h5close(file)
                  return(TRUE)
                },
                show = function() {
                  cat("iData object \n")
                  cat("iList contains: \n")
                  for (i in seq_len(length(iList)))
                    print(iList[[i]])
                  cat("demog contains: \n")
                  cat(" ", ncol(demog), "variables \n")
                  cat(" ", nrow(demog), "rows \n")
                },
                split = function(nsplit) {
                  'splits iData into two groups with specified ratios if nsplit < 1 or into n folds if nsplit > 1'
                  ndemog <- nrow(demog)
                  if (nsplit > 1) {
                    mylist <- list()
                    split_ids <- sample(ndemog) %% round(nsplit) + 1
                    for (j in sort(unique(split_ids))) {
                      splitrow <- seq_len(ndemog)[split_ids == j]
                      ilist <- list()
                      for (i in seq_len(iList)) {
                        ilist[[i]] <- iGroup(name = iList[[i]]$name,
                                             iMatrix = iList[[i]]$iMatrix[splitrow, ],
                                             mask = antsImageClone(iList[[i]]$mask),
                                             modality = iList[[i]]$modality)
                      }
                      names(ilist) <- names(iList)
                      mylist[[j]] <- iData(iList = ilist, demog = demog[splitrow, ], index = index[splitrow, ])
                      names(mylist)[length(mylist)] <- paste("split", j, sep = "")
                    }
                    return(mylist)
                  } else if (nsplit < 1) {
                    splitrow <- sample(1:ndemog, nsplit * ndemog)
                    
                    ilist1 <- list()
                    for (i in 1:length(iList)) {
                      ilist1[[i]] <- iGroup(name = iList[[i]]$name,
                                            iMatrix = iList[[i]]$iMatrix[splitrow, ],
                                            mask = antsImageClone(iList[[i]]$mask),
                                            modality = iList[[i]]$modality)
                      names(ilist1)[length(ilist1)] <- names(iList)[i]
                    }
                    
                    ilist2 <- list()
                    for (i in 1:length(iList)) {
                      ilist2[[i]] <- iGroup(name = iList[[i]]$name,
                                            iMatrix = iList[[i]]$iMatrix[-splitrow, ],
                                            mask = antsImageClone(iList[[i]]$mask),
                                            modality = iList[[i]]$modality)
                      
                      names(ilist2)[length(ilist2)] <- names(iList)[i]
                    }
                    mylist[[1]] <- iData(iList = ilist1, demog = demog[splitrow, ], index = index[splitrow, ])
                    mylist[[2]] <- iData(iList = ilist2, demog = demog[-splitrow, ], index = index[-splitrow, ])
                    names(mylist) <- c("split1", "split2")
                    return(mylist)
                  } else
                    cat("no split occured. \n")
                })
              )

#' @details loads previously saved iData from h5 file
#' @describeIn iData
iDataRead <- function(filename, verbose = TRUE) {
  object <- iData()
  object = object$load(filename = filename, verbose = verbose)
  return(object)
}

#' @details saves iData object to h5 file
#' @describeIn iData
iDataWrite <- function(object, filename, verbose = TRUE) {
  if (class(object) != "iData")
    stop("object must be of class iData.")
  object$save(filename = filename, verbose = verbose)
}


#' iData Model Formulae 
#' 
#' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param iData an object of class \code{\link{iData}} containing data represented in the provided formula
#' @param ...
#' 
#' 
#' 
#' @export iFormula
iFormula <- 
  function(formula, iData, control = ...) {
    lhs <- form1[[2]]
    groups <- all.vars(lhs)
    for (i in seq_len(length(groups))) {
      if (!any(names(iData$iList) == groups[i]))
        stop(paste(y, " is not an iGroup found within iData", sep = ""))
    }
    groupIndex <- iData$getGroupIndex(groups)
    
    rhs <- delete.response(formula)
    vars <- all.vars(rhs)

    out <- list(formula = formula, y = groups, groupIndex = groupIndex, RHS = RHS)
    class(out) <- "iFormula"
    return(out)
  }

model.frame.iFormula <- function(object, iData) {
  model.frame(object$RHS, data = iData$getDemog(object$y))
}

## iModel----
#' @importFrom stats coef
#' @S3method coef iModel
coef.iModel <- function(object) {
  object$beta
}

#' @importFrom stats fitted
#' @S3method fitted iModel
fitted.iModel <- function(object) {
  object$X %*% object$beta
}

#' @importFrom stats model.matrix
#' @S3method model.matrix iModel
model.matrix.iModel <- function(object) {
  object$X
}

#' @importFrom stats residuals
#' @S3method residuals iModel
residuals.iModel <- function(object) {
  object$res
}

iModel <-
  setRefClass("iModel",
              fields = list(
                iData = "iData",
                y = "character",
                frame = "data.frame",
                X = "list",
                XV = "matrix",
                beta = "matrix",
                betaCov = "matrix",
                res = "matrix",
                mrss = "matrix",
                K = "ANY",
                xVi = "list",
                XX = "matrix",
                dims = "list",
                trRV = "numeric",
                trRVRV = "numeric",
                rdf = "numeric",
                fwhm = "numeric",
                resels = "numeric",
                rpvImage = "antsImage",
                nvox = "numeric",
                call = "call"),
              methods = list(
                setBeta_Res = function(x) {
                  'set the model coefficients'
                  KWY <- rftFilter(K, X$W %*% Y)
                  beta <<- X$XX %*% KWY
                  res <<- .res(X$KWX, KWY)
                  mrss <<- colSums(res^2) / trRV
                  # if (missing(x))
                  #   beta <<- XX %*% crossprod(X, iData$imgList[[y]]$imageMatrix)
                  # else
                  #   beta <<- x
                },
                setCall = function(call) {
                  'set call field'
                  call <<- call
                },
                setV2R = function (sample, verbose = NULL) {
                  'estimate voxels in resels space.'
                  smooth <- estSmooth(resid, modData$mask, dims$rdf, scaleResid = FALSE, sample, verbose)
                  fwhm <<- smooth$fwhm
                  rpvImage <<- smooth$rpvImage
                  resels <<- resels(modData$mask, smooth$fwhm)
                },
                setTrRV = function() {
                  'set trace of RV, RVRV, and estimate residual degrees of freedom.'
                  out <- .trRV(KWX, X$V)
                  trRV <<- out$trRV
                  trRVRV <<- out$trRVRV
                  rdf <<- out$rdf
                },
                setW = function(x) {
                  'set the weight matrix'
                  if (missing(weights)) {
                    iV <- sqrt(MASS::ginv(xVi$V))
                    weights <- iV * (abs(iV) > 1e-6)
                  } else {
                    if (class(weights) == "numeric" && length(weights) == nimg)
                      X$W <<- diag(weights)
                    else if (all(dim(weights) == nimg))
                      X$W <<- weights
                    else
                      stop("weights must be a matrix of nimg x nimg or a vector of length nimg")
                  }
                },
                setX = function(x) {
                  'set design matrix'
                  X$X <<- x
                  dims$nimg <<- nrow(X)
                  dims$npred <<- ncol(X)
                  dims$nvox <<- ncol(iData$imgList[[y]]$imageMatrix)
                },
                setXVi = function(method = NULL) {
                  if (is.null(xVi$V))
                    xVi$V <<- diag(ncol(X))
                  if (!is.null(method)) {
                    xVi <<- estNonSphericity()
                  }
                },
                setKWX = function() {
                  'set the pseudoinverse of X'
                  X$KWX <<- .setX(rftFilter(K, W %*% X$X))
                  X$XX <<- .pinvx(KWX)
                },
                setX_V = function() {
                  X$V <<- rftFilter(K, rftFilter(K, W %*% tcrossprod(xVi$V, W)))
                },
                setBetaCov = function() {
                  betaCov <<- XX %*% tcrossprod(XV, XX)
                },
                show = function() {
                  'method for printing rftModel'
                  cat("Random Field Theory model fitted by ", call[[1]])
                  cat("Call: \n")
                  print(call)
                  cat("\nCoefficients: \n")
                  
                  for (i in seq_len(ncol(X))) {
                    cat(colnames(X)[i], "\n")
                  }
                  
                  cat("\nDegrees of interest = ", dims$idf, "\n")
                  cat("Residual degrees of freedom = ", dims$rdf, "\n")
                  cat("Voxels = ", nvox, "\n")
                  cat("FWHM = ", fwhm, "\n")
                  cat("Resels = ", resels, "\n\n")
                })
             )

# setRes = function(x, scale = TRUE) {
#   'return residual forming matrix or set residuals'
#   if (missing(x))
#     x <- iData$imgList[[y]]$imageMatrix - X %*% B
#   mrss <<- colSums(x^2) / dims$rdf
#   if (scale)
#     x <- t(t(x) * (1 / mrss))
#   res <<- x
# },
