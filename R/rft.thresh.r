#' Produces a threshold value based on cluster or voxel level statistics
#' 
#' @param x statistical map of class antsImage 
#' @param pval-probability of false positive
#' @param nvox minimum desired cluster size (in voxels)
#' @param n number of images in conjunction
#' @param fwhm-full width at half maxima
#' @param mask-antsImage mask
#' @param df-degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:
#' \itemize{
#' \item{"T"}{T-field} 
#' \item{"F"}{F-field} 
#' \item{"X"}{Chi-square field"} 
#' \item{"Z"}{Gaussian field}
#' }
#' @param threshType:
#' \itemize{
#'	\item{"cRFT"}{computes a threshold per expected cluster level probability}
#'	\item{"pRFT"}{uses the mask and pval calculates the minimum statistical threshold}
#'	\item{"cFDR"}{uses an uncorrected threshold at the alpha level and then computes and FDR threshold based on cluster maxima}
#'	\item{"pFDR"}{computes the fdr threshold for the entire field of voxels}
#' }
#' @param pp primary threshold used in FDR methods
#' @param verbose enables verbose output
#' @return Outputs a statistical value to be used for threshold a SPM
#' @description
#' 
#'	A statistical threshold level is predicted using a p-value (pval) and 
#'	suprathreshold cluster level (ka). The input statistical parametric map (SPM) 
#'	is then thresholded and clusters are extracted. Random field theory (RFT) is 
#'	used to produce the cluster-level statistics.  
#'	
#'	It is important to note that there is an inverse relationships between the 
#'	'pval' and 'ka' input and the calculated threshold. Calculating the actual
#'	cluster-level statistics utilizes the 'ka' and the threshold value. Therefore,
#'	the 'pval' and 'ka' should be used according to the type of analysis (fMRI,
#'	PET, or VBM) and region of interest. This has been validated with a power analysis
#'	that utilizes the peviously discussed values herein and signal characteristics 
#'	(Friston et al., 1996). Therefore, if the ultimate goal of an analysis is to 
#'	produce cluster-level statistics then lower thresholds may ultimately lead to  
#'	higher statistical power. Alternatively, voxel-level statistics have increased
#'	power with higher thresholds.
#'
#'	Consult provided reference material for further details concerning how to produce
#'	apropriate thesholds for your analysis.
#'
#'	
#' @References
#' Friston K.J., (1994) Assessing the Significance of Focal Activations Using Their Spatial Extent
#' Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power
#' @Author Zachary P. Christensen
#' @note: function currently in beta phase
#' @examples
#'
#'
#'
#' @export rft.thresh
rft.thresh <-function(x, pval, nvox, n, fwhm, resels, df, fieldType, threshType, pp=.001,verbose=FALSE){
  D <-x@dimension
  vox2res <-1/prod(fwhm)
  k <-nvox*vox2res
  nvox <-length(as.vector(x[x !=0]))
  # find minimum threshold value to acheive desired using RFT
  if (threshType=="cRFT" | threshType=="pRFT"){
    u <-max(x)
    if (threshType=="cRFT"){
      alpha <-rft.pval(D, 1, k, u, n, resels, df, fieldType)$Pcor
    }else if(threshType=="pRFT"){
      alpha <-rft.pval(D, 1, 0, u, n, resels, df, fieldType)$Pcor
    }
    if (alpha > pval){
      stop("No voxels survive threshold given the parameters")
    }
    while(alpha < pval){
      u <-u-.01
      if (threshType=="cRFT"){
        alpha <-rft.pval(D, 1, k, u, n, resels, df, fieldType)$Pcor
      }else if(threshType=="pRFT"){
        alpha <-rft.pval(D, 1, 0, u, n, resels, df, fieldType)$Pcor
      }
    }
  }else if (threshType=="cFDR" | threshType=="pFDR"){
    # FDR based thresholds 
    # find initial threshold for given fieldType
    if (fieldType=="Z"){
      stat <-qnorm(1-pp)
    }else if(fieldType=="T"){
      stat <-qt(1-pp, df[2])
    }else if(fieldType=="F"){
      stat <-qf(1-pp, df[1], df[2])
    }else if(fieldType=="X"){
      stat <-qchisq(1-pp, df[1],df[2])
    }
    if (threshType=="cFDR"){
      fdrclust <-labelClusters(x,k,u,Inf)
      fdrlabs <-unique(clust[clust > 0])
      cmax <-c()
      
      if (verbose=="True"){
        progress <-txtProgressBar(min=0, max=n, style=3)
      }
      for (i in 1:length(fdrlabs)){
        cmax <-c(cmax,rft.pval(D, 1, 0, max(x[fdrclust]), n, resels, df, fieldType)$Ec)
        if (verbose=="TRUE"){
          setTxtProgressBar(progress,i)
        }
      }
      if (verbose=="TRUE"){
        close(progress)
      }
    }else if(threshType=="pFDR"){
      fdrclust <-antsImageClone(x)
      fdrclust[fdrclust < stat] <-0
      cmax <-as.numeric(fdrclust)
      if (verbose=="True"){
        progress <-txtProgressBar(min=0, max=n, style=3)
      }
      for (i in 1:length(cmax)){
        cmax[i] <-rft.pval(D, 1, 0, max(x[fdrclust]), n, resels, df, fieldType)$Ec
        if (verbose=="TRUE"){
          setTxtProgressBar(progress,i)
        }
      }
      if (verbose=="TRUE"){
        close(progress)
      }
    }
    cmax <-cmax/rft.pval(D, 1, 0, stat, n, resels, df, fieldType)$Ec
    cmax <-sort(cmax,decreasing=TRUE)
    # find Q
    fdrvec <-sort((pval*(1:length(cmax))/length(cmax)),decreasing=TRUE)
    i <- 1
    while (cmax[i] >= fdrvec[i]){
      i <-i+1
    }
    if (fieldType=="Z"){
      u <-qnorm(1-cmax[i])
    }else if(fieldType=="T"){
      u <-qt(1-cmax[i], df[1])
    }else if(fieldType=="F"){
      u <-qf(1-cmax[i], df[1], df[2])
    }else if(fieldType=="X"){
      u <-qchisq(1-cmax[i], df[1],df[2])
    }
  }
  u
  }
