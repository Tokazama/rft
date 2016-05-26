rftMakeContrast <- function (object, contrast, conType, statdir) {
  con_dims <- dim(contrast)
  if (is.null(con_dims))
    contrast <- as.matrix(contrast, nrow = 1)
  ncon <- con_dims[1]
  if (object$dims[2] != con_dims[2])
    stop("Each contrast length must be equal to the number of columns in the design matrix, including any intercepts. \n")
  if (missing(conType)){
    conType = "T"
  }
  
  if (is.null(rownames(contrast)))
    conNames <- paste("contrastImage", 1:ncon, sep = "")
  else
    conNames <- rownames(contrast)
  
  if (any(conType == "T"))
    XX <- MASS::ginv(crossprod(X))
  
  RV <- object$R %*% object$V
  df2 <- (sum(diag(RV))^2) / sum(diag(RV %*% RV))  # trRV^2 / trRVRV
  mrss <- object$rss / sum(diag(RV))
  ilist <- list()
  for (i in 1:ncon) {
    c <- matrix(contrast[i, ], ncol = 1)
    # t contrast-------------------------------------------------------------------
    if (conType[i] == "T") {
      se <- sqrt(mrss * (crossprod(c, XX) %*% c))
      statvec <- crossprod(c, object$coefficients) / se
      df1 <- 0
    }
    
    # F contrast-------------------------------------------------------------------
    if (conType == "F") {
      C0 <- diag(object$dims[2]) - c %*% MASS::ginv(c)
      X0 <- X %*% C0
      R0 <- diag(object$dims[1]) - X0 %*% MASS::ginv(X0)
      M <- R0 - object$R
      MV <- M %*% object$V
      traceMV <- sum(diag(MV))

      statvec <- ((t(object$coefficients) %*% t(X) %*% M %*% X %*% object$coefficients) / traceMV) / mrss
      
      df1 <- traceMV^2 / sum(diag(MV %*% MV))  # trMV^2 / trMVMV
    }
    ilist <- lappend(ilist, list(contrastImage = makeImage(mask, statvec),
                                 dof = c(df1, df2),
                                 conType = conType[i]))
    if (!is.null(statdir))
      antsImageWrite(ilist[[i]]$contrastImage, file = paste(statdir, conNames, ".nii.gz", sep = ""))
  }
  names(ilist) <- conNames
  ilist
}
