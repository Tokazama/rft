#' Reference class for imgData
#' 
#' Creates an object that holds demographic and corresponding image information
#' in its own environment that allows models to be fitted without creating 
#' redundant copies of large objects. 
#' 
#' @param imageMatrix specifies
#' @param mask
#' @param bool 
#' @param name name representing the image group
#' @param demog data.frame that identifies the variables of interest
#' 
#' 
#' @field imgList list of image groups. Each group contains the following: 
#' \itemize{
#'   \item{imageMatrix} {matrix where each row represents an image}
#'   \item{mask} {antsImage mask where 1s represents active voxels}
#'   \item{bool} {logical vector whose length is equal to the rows in demog. If 
#'   true then the corresponding row is represented in the imageMatrix.}
#' }
#' @field demog demographic information for image groups
#' @examples 
#' 
#'
#' @export imgData
imgData <- 
  setRefClass(Class = "imgData",
              fields = list(imgList = "list", demog = "data.frame"),
              methods = list(
                addGroup = function(imageMatrix, mask, bool, name) {
                  'add image group to imgList field'
                  if (sum(bool) != nrow(imageMatrix))
                    stop('number of true elements in bool does not equal rows in imageMatrix')
                  if (length(bool) != nrow(imgData$demog))
                    stop('length(bool) must equal rows in demog')
                  
                  ilist <- list(imageMatrix = imageMatrix, mask = mask, bool = bool)
                  imgList <<- lappend(imgList, ilist)
                  names(imgList)[length(imgList)] <- name
                  return(imgData)
                },
                check = function(name) {
                  'checks the mask in an image group to ensure only active voxels are represented'
                  if (any(names(imgList) == name)) {
                    mask_vec <- imgList[name]$imageMatrix
                    mask_vec[mask_vec != 0] <- 1
                    mask_vec <- colSums(mask_vec)
                    mask_vec[mask_vec != nimg] <- 0
                    mask_vec[mask_vec == nimg] <- 1
                    imgList[name]$imageMatrix <<- imgList[name]$imageMatrix[, as.logical(mask_vec)]
                    imgList[name]$mask <<- makeImage(imgList[name]$mask, mask_vec)
                  } else
                    stop(paste(name, "is not an image group within imgList"))
                },
                copy = function() {
                  'make copy of imgData object'
                  out <- list()
                  for (i in 1:length(imgList)) {
                    out[[i]]$mask <- antsImageClone(imgList[[i]]$mask)
                    out[[i]]$imageMatrix <- imgList[[i]]$imageMatrix
                    out[[i]]$bool <- imgList[[i]]$bool
                  }
                  out <- imgData(imgList = out, demog = demog)
                  names(out$imgList) <- names(imgList)
                  return(out)
                },
                load = function(filename) {
                  'loads previously saved imgData from h5 file'
                  if (!usePkg("h5"))
                    stop( "Please install package h5 in order to use this function." )
                  if (!file.exists(filename))
                    stop("file does not exist")
                  file <- h5file(filename)
                  tmp <- file["demog/rftDemogMatrix"][]
                  dnam <- h5attr(tmp, "colnames")
                  colnames(tmp) <- dnam
                  for (i in 1:length(dnam)) {
                    colnum <- paste("demog/col", i, sep = "")
                    tmpattr <- list.attributes(file[colnum])
                    for (j in 1:length(tmpattr))
                      attr(tmp[dnam[i]], tmpattr[j]) <- h5attr(file[colnum], tmpattr[j])
                  }
                  demog <<- tmp
                  
                  nilists <- length(list.groups(file["imgList"]))
                  ilists <- paste("imgList/ilist", nlist, sep = "")
                  for (i in 1:nilists) {
                    imgList[[i]]$imageMatrix <<- file[file.path(ilists[i], "imageMatrix")][]
                    imgList[[i]]$bool <<- file[file.path(ilists[i], "bool")][]
                    tmpmask <- file.path(ilists[i], "mask")
                    mymask <<- as.antsImage(file[tmpmask][])
                    k = antsSetSpacing(mymask, h5attr(file[tmpmask], "spacing"))
                    k = antsSetOriginh(mymask, h5attr(file[tmpmask], "origin"))
                    k = antsSetDirection(mymask, h5attr(file[tmpmask], "direction"))
                    imgLists[[i]]$mask <- antsImageClone(mymask)
                    names(imgList)[length(imgList)] <- file[file.path(ilists[i], "name")][]
                  }
                  h5close(file)
                },
                save = function(filename) {
                  'save imgData to h5 file'
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
                  file["demog/rftDemogMatrix"] <- data.matrix(demog)
                  dnam <- colnames(demog)
                  h5attr(file["demog/rftDemogMatrix"], "colnames") <- dnam
                  for (i in 1:length(dnam)) {
                    dattrfile <- file[paste("demog/col", i, sep = "")]
                    dattr <- names(attributes(demog[dnam[i]]))
                    if (!is.null(dattr)) {
                      for (j in 1:length(dattr))
                        h5attr(dattrfile, dattr[j]) <- attr(demog[dnam[i]], dattr[j])
                    }
                  }
                  
                  # save imgList
                  nlist <- length(imgList)
                  ilists <- paste("imgList/ilist", nlist, sep = "")
                  inames <- names(imgList)
                  for (i in 1:nlist) {
                    file[file.path(ilists[i], "name")] <- inames[i]
                    file[file.path(ilists[i], "imageMatrix")] <- imgList[[i]]$imageMatrix
                    file[file.path(ilists[i], "bool")] <- imgList[[i]]$bool
                    tmpmask <- file.path(ilists[i], "mask")
                    file[tmpmask] <- as.array(imgList[[i]]$mask)
                    h5attr(file[tmpmask], "spacing") <- antsGetSpacing(imgList[[i]]$mask)
                    h5attr(file[tmpmask], "direction") <- antsGetDirection(imgList[[i]]$mask)
                    h5attr(file[tmpmask], "origin") <- antsGetOrigin(imgList[[i]]$mask)
                  }

                  h5close(file)
                  return(TRUE)
                },
                show = function() {
                  'concise print of modData'
                  cat("imgData: \n\n")
                  
                  cat("imgList: \n")
                  for (i in 1:length(imgList))
                    cat(" ", paste(names(imgList)[i]), "\n")
                  cat("\n")
                  
                  cat("demog: \n")
                  ndemog <- ncol(demog)
                  cat(" ", ndemog, "variables \n")
                  if (ndemog <= 10) {
                    for (i in 1:ndemog)
                      cat(paste(names(demog)[i]), "\n")
                    cat("\n")
                  }
                },
                split = function(nsplit) {
                  'splits imgData into two groups with specified ratios if nsplit < 1 or into n folds if nsplit > 1'
                  ndemog <- nrow(demog)
                  if (nsplit > 1) {
                    mylist <- list()
                    split_ids <- sample(ndemog) %% round(nsplit) + 1
                    for (j in sort(unique(split_ids))) {
                      splitrow <- (1:ndemog)[split_ids == j]
                      ilist <- list()
                      for (i in 1:length(imgList)) {
                        ilist[[i]] <- list(imageMatrix = imgList[[i]]$imageMatrix[splitrow, ],
                                      mask = antsImageClone(imgList[[i]]$mask),
                                      bool = imgList[[i]]$bool[splitrow])
                        
                        names(ilist)[length(ilist)] <- names(imgList)[i]
                      }
                      mylist[[j]] <- imgData(imgList = ilist, demog = demog[splitrow, ])
                      names(mylist)[length(mylist)] <- paste("split", j, sep = "")
                    }
                    return(mylist)
                  } else if (nsplit < 1) {
                    splitrow <- sample(1:ndemog, nsplit * ndemog)
                    
                    ilist1 <- list()
                    for (i in 1:length(imgList)) {
                      ilist1[[i]] <- list(imageMatrix = imgList[[i]]$imageMatrix[splitrow, ],
                                         mask = antsImageClone(imgList[[i]]$mask),
                                         bool = imgList[[i]]$bool[splitrow])
                      
                      names(ilist1)[length(ilist1)] <- names(imgList)[i]
                    }
                    
                    ilist2 <- list()
                    for (i in 1:length(imgList)) {
                      ilist2[[i]] <- list(imageMatrix = imgList[[i]]$imageMatrix[-splitrow, ],
                                          mask = antsImageClone(imgList[[i]]$mask),
                                          bool = imgList[[i]]$bool[-splitrow])
                      
                      names(ilist2)[length(ilist2)] <- names(imgList)[i]
                    }
                    mylist[[1]] <- imgData(imgList = ilist1, demog = demog[splitrow, ])
                    mylist[[2]] <- imgData(imgList = ilist2, demog = demog[-splitrow, ])
                    names(mylist) <- c("split1", "split2")
                    return(mylist)
                  } else
                    cat("no split occured. \n")
                }
                varBind = function(newVars) {
                  'add variables to existing demog data.frame'
                  if (nrow(newVars) != nrow(demog))
                    stop("nrow(newVars does not equal number of rows in demog")
                  demog <<- cbind(demog, newVars)
                }
                )
              )

imgDataMake <- function(imageMatrix, mask, bool, name, demog) {
  if (sum(bool) != nrow(imageMatrix))
    stop('number of true elements in bool does not equal rows in imageMatrix')
  if (length(bool) != nrow(demog))
    stop('length(bool) must equal rows in demog')
  
  ilist <- list(imageMatrix = imageMatrix, mask = mask, bool = bool)
  imgData(imgList = list(name = ilist), demog = demog)
}
