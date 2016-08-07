# to do:

#' iData/iGroup reference class
#' 
#' @details 
#' iData/iGroup refers to a set of S4 classes used for representing image data.
#' The base compositional unit is the \code{\link{iGroup}} class that allows a 
#' set of images to be represented by an image matrix, name, binary mask (of 
#' class \code{\link{antsImage}}), information describing the images' modality
#' (i.e. fMRI, PET, VBM, etc.), and optionally the image filter information. 
#' 
#' @usage
#' iGroup(iMatrix, name, mask, modality, ...)
#' 
#' iData(iList, demog, index, ...)
#' 
#' @example 
#' 
#' # create iGroup object
#' ilist <- getANTsRData("population")
#' mask <- getMask(ilist[[1]])
#' imat <- imageListToMatrix(ilist, mask)
#' iGroup1 <- iGroup(imat, "pop1", mask, modality = "T1")
#' 
#' # ensure only active voxels are included in mask
#' 
#' 
#' ilist <- lappend(ilist, ilist[[1]])
#' imat <- imageListToMatrix(ilist, mask)
#' iGroup2 <- iGroup(imat, "pop2", mask, modality = "T1")
#' 
#' # save iGroup object
#' tmpfile <- tempfile(fileext = ".h5")
#' iGroupWrite(iGroup1, tmpfile)
#' 
#' # load saved iGroup object
#' (iGroup1_reload <- iGroupRead(tmpfile))
#' 
#' demog <- data.frame(id = c("A", "B", "C", NA),
#'   age = c(11, 7, 18, 22), sex = c("M", "M", "F", "F"))
#'   
#' bool1 <- c(TRUE, TRUE, TRUE, FALSE)
#' bool2 <- c(TRUE, TRUE, TRUE, TRUE)
#' 
#' # create iData object that holds demographics info
#' mydata <- iData(iGroup1, bool1, demog)
#' 
#' # add iGroup object to iData
#' mydata <- add(mydata, iGroup2, bool1)
#' 
#' # save iData object
#' tmpdir <- "iData_test"
#' iDataWrite(mydata, tmpdir)
#' 
#' # load saved iData object
#' (mydata_reload <- iDataRead(tmpdir))
#' 
#' # split iData object into k-folds or train and test groups
#' mydatasplit <- iDataSplit(mydata, 0.3)
#' 
#' # retreive demographic information specific to an iGroup
#' getDemog(mydata, "pop1", c("age", "sex"))
#'
#' @export iGroup/iData-class


# iGroup----
#' iGroup-class
#' 
#' @param .Object
#' @param iMatrix image by voxel matrix
#' @param name name for the iGroup object (used for reference in \code{\link{iData}} and formulas)
#' @param mask mask with number of non-zero entries defining the matrix columns.
#' @param modality image modality that iGroup will represent
#' @param rowslist row of iMatrix constituting block/partitions
#' @param HParam cut-off period in seconds
#' @param RT observation interval in seconds
#' @param checkMask logical. ensure mask only represents active voxels
#' @param filename optional filename to save iGroup object (default = \code{tempfile(fileext = ".h5")})
#' 
#' @slot name iGroup name
#' @slot file h5file connection
#' @slot iMatrix h5file pointer to matrix
#' @slot location h5file location
#' @slot modality image modality represented
#' @slot mask mask with number of non-zero entries defining the matrix columns.
#' @slot K filter information
#' 
#' @details 
#' See \linke{iGroup/iData-class} for examples
#' 
#' 
#' @export iGroup
iGroup <- setClass("iGroup",
                   slots = list(
                     name = "character",
                     file = "H5File",
                     iMatrix = "DataSet",
                     location = "character",
                     modality = "character",
                     mask = "antsImage",
                     K = "ANY",
                     dim = "numeric")
                   )

setMethod("initialize", "iGroup", function(.Object, iMatrix = matrix(1, 1, 1), name, mask,
                                           modality, rowslist, HParam, RT, checkMask = TRUE, filename) {
  if (!usePkg("h5"))
    stop("Please install package h5 in order to use this function.")
  
  if (!missing(iMatrix) && !missing(mask) && checkMask) {
    nimg <- nrow(iMatrix)
    mask_vec <- iMatrix
    mask_vec[mask_vec != 0] <- 1
    mask_vec <- colSums(mask_vec)
    mask_vec[mask_vec != nimg] <- 0
    mask_vec[mask_vec == nimg] <- 1
    iMatrix <- iMatrix[, as.logical(mask_vec)]
    mask <- makeImage(mask, mask_vec)
  }
  
  # filename
  if (missing(filename)) {
    tmpfile <- tempfile(fileext = ".h5")
    file <- h5file(tmpfile)
    .Object@location <- tmpfile
  } else {
    if (file.exists(filename))
      stop("filename already exists.")
    file <- h5file(filename)
    .Object@location <- filename
  }
  .Object@file <- file
  
  # modality
  if (missing(modality)) {
    file["modality"] <- "fMRI"
    .Object@modality <- "fMRI"
  } else {
    file["modality"] <- modality
    .Object@modality <- modality
  }
  
  # name
  if (missing(name)) {
    file["name"] <- "unnamed"
     .Object@name <- "unnamed"
  } else {
    file["name"] <- name
    .Object@name <- name
  }
  
  # mask
  ## configure
  if (missing(mask))
    .Object@mask <- makeImage(c(1, 1, 1), 1)
  else
    .Object@mask <- antsImageClone(mask)
  ## write
  file["mask"] <- as.array(.Object@mask)
  h5attr(file["mask"] , "spacing") <- antsGetSpacing(.Object@mask)
  h5attr(file["mask"] , "direction") <- antsGetDirection(.Object@mask)
  h5attr(file["mask"] , "origin") <- antsGetOrigin(.Object@mask)
  
  # K
  ## configure
  if (!missing(rowslist) | !missing(HParam) | !missing(RT)) {
    if (missing(rowslist) | missing(HParam) | missing(RT))
      stop("rowslist, HParam, and RT must all be provided if any one is.")
    if ((length(rowslist) != length(HParam)) | (length(rowslist) != length(RT)))
      stop("Length of rowslist, HParam, and RT must be the same.")
    K <- list()
    nk <- length(rowslist)
    K <- data.frame(Filters = rep("F1", nk), HParam = rep(1, nk), RT = rep(1, nk))
    for (i in seq_len(nk)) {
      K$Filters[rowslist[[i]]] <- paste("F", i, sep = "")
      K$HParam[rowslist[[i]]] <- HParam[i]
      K$RT[rowslist[[i]]] <- RT[i]
    }
  } else 
    K <- 1
  .Object@K <- K
  ## write
  if (class(.Object@K) == "list") {
    nk <- length(.Object@K)
    file[file.path("K", "Filters")] <- .Object@K$Filters
    file[file.path("K", "HParam")] <- .Object@K$HParam
    file[file.path("K", "RT")] <- .Object@K$RT
  } else {
    file["K"] <- .Object@K
  }
  
  # iMatrix
  .Object@dim <- dim(iMatrix)
  if (missing(iMatrix)) {
    file["iMatrix"] <- matrix(1, 1, 1)
  } else {
    chunk <- iControl()$chunksize(.Object@dim[1])
    if (chunk < .Object@dim[2]) {
      file["iMatrix"] <- iMatrix[, seq_len(chunk)]
      file["iMatrix"] <- cbind(file["iMatrix"], (chunk + 1):.Object@dim[2])
    } else
      file["iMatrix"] <- iMatrix
  }
  .Object@iMatrix <- file["iMatrix"]
  return(.Object)
})

setMethod("show", "iGroup", function(object) {
  cat("iGroup object:\n")
  cat("     name =", object@name, "\n")
  cat("   images =", nrow(object), "\n")
  cat("   voxels =", ncol(object), "\n")
  cat(" location =", object@location, "\n")
  cat(" modality =", object@modality, "\n")
  cat("___\n")
})

#' @details retrieve dimensions of iGroup's iMatrix slot
#' @describeIn iGroup
setMethod("dim", "iGroup", function(x) {
  return(x@dim)
})

#' @details retrieve name of iGroup object
#' @describeIn iGroup
setMethod("names", "iGroup", function(x) {
  return(x@name)
})

#' @details set name of iGroup
#' @describeIn iGroup
setMethod("names<-", "iGroup", function(x, value) {
    x@name <- value
    x@file["name"][] <- value
    return(x)
})

#' @details subset iGroup objects
#' @describeIn iGroup
setMethod("[", "iGroup", function(x, i) {
  out <- iGroup(x@iMatrix[i, ], x@name, mask, modality = x@modality, checkMask = FALSE)
  if (class(x@K) == "numeric")
    out@K <- 1
  else
    out@K <- x@K[i, ]
  return(out)
})

#' @details read/load iGroup from h5 file
#' @describeIn iGroup
iGroupRead <- function(filename) {
  if (!usePkg("h5"))
    stop( "Please install package h5 in order to use this function." )
  if (!file.exists(filename))
    stop("file does not exist")
  file <- h5file(filename)
  out <- iGroup()
  
  out@location <- filename
  out@iMatrix <- file["iMatrix"]
  out@modality <- file["modality"][]
  out@name <- file["name"][]
  
  # mask
  mymask <- as.antsImage(file["mask"][])
  k = antsSetSpacing(mymask, h5attr(file["mask"], "spacing"))
  k = antsSetOrigin(mymask, h5attr(file["mask"], "origin"))
  k = antsSetDirection(mymask, h5attr(file["mask"], "direction"))
  out@mask <- mymask
  
  # K
  ## configure
  if (file["K"][] == 1)
    out@K <- 1
  else {
    knames <- unique(basename(dirname(list.datasets(file[tmpk]))))
    Filters <- file[file.path("K", "rows")][]
    HParam <- file[file.path("K", "Hparam")][]
    RT <- file[file.path("K", "RT")][]
    out@K <- data.frame(Filters, HParam, RT)
  }
  out@dim <- dim(out@iMatrix[])
    
  return(out)
}

#' @details write/save iGroup to h5 file
#' @describeIn iGroup
iGroupWrite <- function(x, filename) {
  oldfile <- x@location
  if (file.exists(filename))
    stop("filename already exists.")
  out <- iGroup(x@iMatrix[], x@name, x@mask, modality = x@modality, checkMask = FALSE, filename = filename)
  out@K <- x@K
  if (class(x@K) == "data.frame") {
    out@file[file.path("K", "Filters")] <- x@K$Filters
    out@file[file.path("K", "HParam")] <- x@K$HParam
    out@file[file.path("K", "RT")] <- x@K$RT
  }
  
  return(TRUE)
}

# iData----
#' 
#' @param .Object
#' @param x
#' @param bool
#' @param demog
#' 
#' @slot iList list of iGroup objects
#' @slot demog demographic information
#' @slot index index that coordinates iGroup rows with demographic rows
#' 
#' @details 
#' See \link{iGroup/iData-class} for examples
#' 
#' @export iData
iData <- setClass("iData",
                  slots = list(
                    iList = "list",
                    demog = "data.frame",
                    index = "data.frame")
                  )

setMethod("initialize", "iData", function(.Object, x, bool, demog) {
  if (missing(x)) {
    iList <- list(iGroup())
    names(iList) <- iList[[1]]@name
    .Object@iList <- iList
    .Object@index <- data.frame()
    .Object@demog <- data.frame()
  } else {
    # check x
    if (class(x) == "iGroup")
      x <- list(x)
    if (class(x) != "list")
      stop("x must be an iGroup object or list of iGroup objects.")
    lname <- c()
    for (i in seq_len(length(x))) {
      lname <- c(lname, x[[i]]@name)
      if (class(x[[i]]) != "iGroup")
        stop("All members of x must be of class iGroup.")
    }
    names(x) <- lname
    .Object@iList <- x
    
    # check bool
    if (missing(bool))
      .Object@index <- data.frame()
    else {
      if (class(bool) == "logical")
        bool <- list(bool)
      if (length(x) != length(bool))
        stop("Each iGroup object must have a corresponding bool vector listed.")
      
      if (missing(demog))
        .Object@index <- data.frame()
      else {
        for (i in seq_len(length(x))) {
          if (!missing(demog)) {
            if (length(bool[[i]]) != nrow(demog))
              stop(paste("The number of elements in bool[[", i, "]] does not equal the number of rows in demog.", sep = ""))
          .Object@demog <- demog
          }
          tmpsum <- sum(bool[[i]])
          if (nrow(x[[i]]) != tmpsum)
            stop(paste("The number of true elements in bool does not equal number of rows in iGroup object ", x[[i]]@name, ".", sep = ""))
          tmpbool <- as.numeric(bool[[i]])
          tmpbool[tmpbool == 1] <- seq_len(tmpsum)
          if (i == 1)
            index <- data.frame(tmpbool)
          else 
            index <- cbind(index, tmpbool)
        }
        colnames(index) <- lname
        .Object@index <- index
      }
    }
  }
  return(.Object)
})

setMethod("show", "iData", function(object) {
  cat("iData object \n")
  cat("iList contains: \n")
  for (i in seq_len(length(object@iList)))
    print(object@iList[[i]])
  cat("demog contains: \n")
  cat(" ", ncol(object@demog), "variables \n")
  cat(" ", nrow(object@demog), "rows \n")
})

#' @details retrieve names of iGroups in iList slot
#' @describeIn iData
setMethod("names", "iData", function(x) {
    out <- c()
    for (i in seq_len(length(x@iList)))
      out <- c(out, x@iList[[i]]@name)
    return(out)
})

#' @details replace names of iGroups within iList slot
#' @describeIn iData
setMethod("names<-", "iData", function(x, value) {
  if (length(value) != length(x@iList))
    stop("names must be the same length as the length of iList")
  for (i in seq_len(length(x@iList)))
    names(x@iList[[i]]) <- value
  names(x@iList) <- value
  return(x)
})

#' @details add iGroup to iList field
#' @describeIn iData
add <- function(x, iGroup, bool) {
  if (class(iGroup) != "iGroup")
    stop("iGroup must be of class iGroup.")
  if (length(bool) != nrow(x@demog))
    stop("The number of rows in demog must be equal to the length of bool.")
  if (sum(bool) != nrow(iGroup@iMatrix))
    stop("Number of TRUE elements in bool is not equal to number of images in iGroup")
  if (any(names(x@iList) == iGroup@name))
    stop(paste("iGroup of name ", iGroup@name, " already exists in iList", sep = ""))
  
  bool[bool == TRUE] <- seq_len(sum(bool))
  if (any(dim(x@index) == 0)) {
    index <- data.frame(bool)
    colnames(index) <- iGroup@name
  } else {
    index <- cbind(x@index, bool)
    colnames(index)[ncol(index)] <- iGroup@name
  }
  x@index <- index
  
  if (length(x@iList) != 0) {
    x@iList <- lappend(x@iList, iGroup)
    names(x@iList)[length(x@iList)] <- iGroup@name
  } else {
    x@iList <- list(iGroup)
    names(x) <- iGroup@name
  }
  return(x)
}

#' @details drop iGroup objects with provided names from iList
#' @describeIn iData
substract <- function(x, ...) {
  groups <- c(...)
  for (i in seq_len(length(groups))) {
    x@iList <- x@iList[[-which(names(x@iList) == groups[i])]]
    x@index <- x@index[, -which(names(x@index) == groups[i])]
  }
  return(x)
}

#' @details get variables indexed according to iGroup indicated by groups
#' @describeIn iData
getDemog <- function(x, groups, vars) {
  if (!any(groups == names(x)))
    stop("groups does not match any iGroups in iList slot.")
  if (length(groups) == 1) {
    rindex <- as.logical(x@index[groups])
  } else {
    for (i in seq_len(length(groups))) {
      if (i == 1)
        rindex <- x@index[groups[i]]
      else
        rindex <- rindex * x@index[groups[i]]
    }
    rindex <- as.logical(rindex[, 1])
  }
  if (missing(vars))
    return(x@demog[rindex, ])
  else
    return(x@demog[rindex, vars])
}

#' @details retrieve list of iGroups
#' @describeIn iData
getGroups <- function(x, ...) {
  if (class(x) != "iData")
    stop("x must be of class iData.")
  ll <- c(...)
  return(x@iList[ll])
}

#' @details create index for retreiving complimentary iGroup subjects
#' @describeIn iData
getGroupIndex <- function(x, ...) {
  out <- x@index[c(...)]
  for (i in seq_len(nrow(x@index))) {
    if (any(is.na(out[i, ])))
      out[i, ] <- out[i, ] * -1
  }
  return(out)
}

#' @details subset iData objects
#' @describeIn iData
setMethod("[", "iData", function(x, i) {
  out <- iData(x@iList[[1]][x@index[i, 1]], as.logical(x@index[i, 1]), x@demog[i, ])
  lg <- length(x@iList)
  if (lg > 1) {
    for (j in 2:lg)
      out <- add(out, x@iList[[j]][x@index[i, j]], as.logical(x@index[i, j]))
  }
  return(out)
})

#' @details loads previously saved iData from its set directory
#' @describeIn iData
iDataRead <- function(dirname, verbose = TRUE) {
  if (!usePkg("h5"))
    stop( "Please install package h5 in order to use this function." )
  if (!dir.exists(dirname))
    stop("dirname does not exist.")
  if (!dir.exists(file.path(dirname, "iList")))
    stop("dirname does not appear to be the directory of an iData object.")
  if (!file.exists(paste(dirname, "/demog.h5", sep = "")))
    stop("dirname does not appear to be the directory of an iData object.")

  # read each iGroup
  if (verbose)
    cat("Reading iList...")
  x <- list()
  
  inames <- c()
  ifiles <- list.files(file.path(dirname, "iList"))
  ifiles <- file.path(dirname, "iList", ifiles)
  for (i in seq_len(length(ifiles))) {
    x[[i]] <- iGroupRead(ifiles[i])
    inames <- c(inames, x[[i]]@name)
    if (verbose)
      cat( "\n ", inames[i])
  }
  names(x) <- inames
  if (verbose)
    cat(". \n")
  
  # read index
  if (verbose)
    cat("Reading index. \n")
  file <- h5file(file.path(dirname, "demog.h5"))
  index <- data.frame(file["index"][])
  colnames(index) <- h5attr(file["index"], "colnames")
  
  # read demog
  if (verbose)
    cat("Reading demog. \n")
  demog <- as.data.frame(file["demog/matrix"][])
  dnam <- h5attr(file["demog/matrix"], "colnames")
  colnames(demog) <- dnam
  for (i in seq_len(length(dnam))) {
    dattrfile <- file[paste("demog/col", i, sep = "")]
    dattr <- list.attributes(dattrfile)
    if (!is.null(dattr)) {
      for (j in seq_len(length(dattr)))
        attr(demog[dnam[i]], dattr[j]) <- h5attr(dattrfile, dattr[j])
    }
  }
  
  object <- iData()
  object@iList <- x
  object@demog <- demog
  object@index <- index
  
  return(object)
}

#' @details write/save iData object to its own directory
#' @describeIn iData
iDataWrite <- function(x, dirname, verbose = TRUE) {
  if (dir.exists(dirname))
    stop("dirname already exists")
  dir.create(dirname)
  
  # write iList
  if (verbose)
    cat("Writing iGroup... \n")
  dir.create(file.path(dirname, "iList"))
  for (i in seq_len(length(x@iList))) {
    cat(" ", x@iList[[i]]@name, "\n")
    tmpname <- file.path(dirname, "iList", paste("iGroup", i, ".h5", sep = ""))
    iGroupWrite(x@iList[[i]], tmpname)
  }
  
  file <- h5file(file.path(dirname, "demog.h5"))
  # write index
  if (verbose)
    cat("Writing index. \n")
  file["index"] <- data.matrix(x@index)
  h5attr(file["index"], "colnames") <- colnames(x@index)
  
  # write demog
  if (verbose)
    cat("Writing demog. \n")
  file["demog/matrix"] <- data.matrix(demog)
  dnam <- colnames(demog)
  h5attr(file["demog/matrix"], "colnames") <- dnam
  for (i in seq_len(length(dnam))) {
    dattrfile <- file[paste("demog/col", i, sep = "")]
    dattr <- names(attributes(demog[dnam[i]]))
    if (!is.null(dattr)) {
      for (j in seq_len(length(dattr)))
        h5attr(dattrfile, dattr[j]) <- attr(demog[dnam[i]], dattr[j])
    }
  }
  
  h5close(file)
  return(TRUE)
}

#' @details splits iData into two groups with specified ratios if nsplit < 1 or into n folds if nsplit > 1
#' @describeIn iData
iDataSplit <- function(x, nsplit) {
  ndemog <- nrow(x@demog)
  if (nsplit > 1) {
    out <- list()
    split_ids <- sample(ndemog) %% round(nsplit) + 1
    for (j in sort(unique(split_ids))) {
      splitrow <- seq_len(ndemog)[split_ids == j]
      out[[j]] <- x[splitrow]
      names(out)[length(out)] <- paste("split", j, sep = "")
    }
  } else if (nsplit < 1) {
    splitrow <- sample(1:ndemog, nsplit * ndemog)
    out <- list(split1 = x[splitrow], split2 = x[-splitrow])
  } else
    cat("no split occured. \n")
  return(out)
}