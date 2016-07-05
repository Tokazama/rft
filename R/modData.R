#' Reference class for modData
#' 
#' Creates an object that holds demographic and corresponding image information
#' in its own environment that allows models to be fitted without creating 
#' redundant copies of large objects. 
#' 
#' @param imageMatrix specifies
#' @param mask mask with number of non-zero entries defining the matrix columns
#' @param demog data.frame that identifies the variables of interest
#' 
#' 
#' @field imageMatrix corresponds to imageMatrix argument 
#' @field mask corresponds to mask argument
#' @field demog corresponds to demog argument
#' @field nvox number of voxels in image
#' @field nimg number of images represented by imageMatrix
#' @examples 
#' 
#' \dontrun {
#'   (mydata <- modData(imageMatrix = mat, mask = mask, demog = dframe))
#'   
#'   # ensure that only active voxels are included in mask and imageMatrix
#'   mydata$check()
#'   
#'   # add new images with respective demographics data
#'   mydata$addImages(newmat, newdemog)
#'   
#'   # add new variables to existing demographics data
#'   mydata$addVar(newvars)
#'   
#'   # save data as h5 file
#'   mydata$save(tempfile(fileext=".h5"))
#'   
#'   # load existing h5 file
#'   mydata$load(tempfile(fileext=".h5"))
#'   
#'   # split mydata into testing and training group
#'   (mydatasplit <- mydata$split(2/3))
#'   
#'   # split mydata into multiple groups
#'   (mydatasplit <- mydata$split(4))
#'   
#'   # bind data back into single group
#'   (mydata <- modData$bind(mydatasplit))
#'   
#'   # export fields as list
#'   listmydata <- mydata$getAll()
#'   
#' }
#'
#' @export modData
modData <- 
  setRefClass(Class = "modData",
              fields = list(
                imageMatrix = "ANY",
                mask = "antsImage",
                demog = "data.frame",
                nimg = "numeric",
                nvox = "numeric"))

modData$methods(
  addImages = function(newImageMatrix, newDemog) {
    'add new images to imageMatrix and corresponding demographics data to demog'
    imgdim <- dim(newImageMatrix)
    if (imgdim[2] != nvox)
      stop("newImageMatrix and imageMatrix do not have the same number of voxels")
    if (!missing(newDemog)) {
      demogdim <- dim(newDemog)
      if (imgdim[1] != demogdim[1])
        stop("rows of newImageMatrix and newDemog are not equal")
      else if (colnames(newDemog) != colnames(demog))
        stop("newDemog does not have the same variables as demog")
      else
        demog <<- rbind(demog, newDemog)
    } else {
      for (i in 1:imgdim[1])
        demog <<- rbind(demog, NA)
    }
    imageMatrix <<- rbind(imageMatrix, newImageMatrix)
  },
  addVar = function(newVars) {
    'add variables to existing demog data.frame'
    if (nrow(newVars) != nrow(demog))
      stop("nrow(newVars does not equal number of rows in demog")
    demog <<- cbind(demog, newVars)
  },
  bind = function(modDataSplit) {
    'binds modDataSplit object into a modData object'
    if (class(modDataSplit) != "modDataSplit")
      stop("modDataSplit must be of class modDataSplit")
    nimg <<- attr(modDataSplit, "imgTotal")
    nvox <<- modDataSplit$split1$nvox
    imat <- matrix(nrow = nimg, ncol = nvox)
    dframe <- as.data.frame(matrix(nrow = nimg, ncol = ncol(modDataSplit$split1$demog)))
    colnames(dmat) <- colnames(modDataSplit$split1$demog)
    for (i in 1:length(modDataSplit)) {
      splitnum <- paste("split", i, sep = "")
      splitrow <- attr(modDataSplit, splitnum)
      imat[splitrow, ] <- modDataSplit[splitnum]$imageMatrix
      dframe[splitrow, ] <- modDataSplit[splitnum]$demog
    }
    imageMatrix <<- imat
    demog <<- dframe
  },
  check = function() {
    'check imageMatrix and mask for non-active voxels'
    nimg <<- nrow(imageMatrix)
    nvox <<- ncol(imageMatrix)
    if (nrow(demog) != nimg)
      stop("nrow(demog) is not equal to nrow(imageMatrix)")
    mask_vec <- as.matrix(imat)
    mask_vec[mask_vec != 0] <- 1
    mask_vec <- colSums(mask_vec)
    mask_vec[mask_vec != nimg] <- 0
    mask_vec[mask_vec == nimg] <- 1
    tmp <- as.matrix(imageMatrix)[, as.logical(mask_vec)]
    
    nvox <<- ncol(tmp)
    imageMatrix <<- tmp
    mask <<- makeImage(mask, mask_vec)
  },
  copy = function() {
    'copy modData object'
    return(modData(imageMatrix = imageMatrix, mask = mask, demog = demog))
  },
  getAll = function() {
    'returns list of all fields'
    list(imageMatrix = imageMatrix, mask = mask, demog = demog, nimg = nimg, nvox = nvox)
  },
  show = function() {
    'concise print of modData'
    cat("modData: \n")
    
    cat("  imageMatrix: \n")
    cat("   ", nimg, "images \n")
    cat("   ", nvox, "voxels \n")
    
    cat("  demog: \n")
    cat("   ", ncol(demog), "variables ")
  },
  load = function(filename) {
    'loads previously saved modData from h5 file'
    if (!usePkg("h5"))
      stop( "Please install package h5 in order to use this function." )
    if (!file.exists(filename))
      stop("file does not exist")
    file <- h5file(filename)
    tmp <- file["modData/demog/rftDemogMatrix"][]
    dnam <- h5attr(tmp, "colnames")
    colnames(tmp) <- dnam
    for (i in 1:length(dnam)) {
      colnum <- paste("modData/demog/col", i, sep = "")
      tmpattr <- list.attributes(file[colnum])
      for (j in 1:length(tmpattr))
        attr(tmp[dnam[i]], tmpattr[j]) <- h5attr(file[colnum], tmpattr[j])
    }
    demog <<- tmp
    imageMatrix <<- as.antsMatrix(file["modData/imageMatrix"][], "double")
    temp <- file["modData/mask"]
    masktmp <- as.antsImage( temp[] )
    k=antsSetSpacing( mask, h5attr( temp, "spacing" ) )
    k=antsSetOrigin( mask, h5attr( temp, "origin" ) )
    k=antsSetDirection( mask, h5attr( temp, "direction" ) )
    mask <<- masktmp
    h5close(file)
  },
  save = function(filename) {
    'save modData to h5 file'
    if ( ! usePkg( "h5" ) )
      stop( "Please install package h5 in order to use this function." )
    if (file.exists(filename)) {
      if (overwrite) {
        warning("will overwrite existing filenam")
        file.remove(filename)
      } else
        stop( "stopping because file.exists( filename )" )
    }
    file <- h5file(filename)
    # save demographics data.frame and attributes
    file["modData/demog/rftDemogMatrix"] <- data.matrix(demog)
    dnam <- colnames(demog)
    h5attr(file["modData/demog/rftDemogMatrix"], "colnames") <- dnam
    for (i in 1:length(dnam)) {
      dattrfile <- file[paste("modData/demog/col", i, sep = "")]
      dattr <- names(attributes(demog[dnam[i]]))
      if (!is.null(dattr)) {
        for (j in 1:length(dattr))
          h5attr(dattrfile, dattr[j]) <- attr(demog[dnam[i]], dattr[j])
      }
    }
    # save imageMatrix
    file["modData/imageMatrix"] <- imageMatrix
    # save mask
    file["modData/mask"] <- as.array( imageMask )
    h5attr(file["modData/mask"], "spacing") <- antsGetSpacing( imageMask )
    h5attr(file["modData/mask"], "direction") <- antsGetDirection( imageMask )
    h5attr(file["modData/mask"], "origin") <- antsGetOrigin( imageMask )
    file["modData/filename"] <- filename
    h5close(file)
    return(TRUE)
  },
  split = function(nsplit) {
    'splits modData into two groups with specified ratios if nsplit < 1 or into n folds if nsplit > 1'
    if (nsplit > 1) {
      mylist <- list()
      split_ids <- sample(nimg) %% round(nsplit) + 1
      for (i in sort(unique(split_ids))) {
        splitrow <- (1:nimg)[split_ids == i]
        mylist[[paste("split", i, sep = "")]] <- modData(imageMatrix = imageMatrix[splitrow, ],
                                                         mask = mask, demog = demog[splitrow, ])
        attr(mylist, paste("split", i, sep = "")) <- splitrow
      }
      attr(mylist, "imgTotal") <- nimg
      class(mylist) <- "modDataSplit"
      return(mylist)
    } else if (nsplit < 1) {
      splitrow <- sample(1:nimg, nsplit * nimg)
      mylist <- list(split1 = modData(imageMatrix = imageMatrix[splitrow, ],
                                      mask = mask, demog = demog[splitrow, ]),
                     split2 = modData(imageMatrix = imageMatrix[-spitrow, ],
                                      mask = mask, demog = demog[-spitrow, ]))
      attr(mylist, "imgTotal") <- nimg
      attr(mylist, "split1") <- splitrow
      attr(mylist, "split2") <- 1:nimg[-splitrow]
      class(mylist) <- "modDataSplit"
      return(mylist)
    } else
      cat("no split occured. \n")
  })
