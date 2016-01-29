#' Produces a threshold value based on cluster or voxel level statistics
#' 
#'
#' @param img statistical map of class antsImage 
#' @param pval-probability of false positive
#' @param ka-minimum desired cluster size
#' @param fwhm-full width at half maxima
#' @param mask-antsImage mask
#' @param df-degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:
#' \describe{
#' \item{"T"}{T-field} 
#' \item{"F"}{F-field} 
#' \item{"X"}{Chi-square field"} 
#' \item{"Z"}{Gaussian field}
#' }
#' @param threshType:
#' \describe{
#'	\item{"cluster"}{computes a threshold per expected cluster level probability}
#'	\item{"voxel"}{uses the mask and pval calculates the minimum statistical threshold}
#'	\item{"cfdr"}{uses an uncorrected threshold at the alpha level and then computes and FDR threshold based on cluster maxima}
#'	\item{"vfdr"}{computes the fdr threshold for the entire field of voxels}
#' }
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
#' @note: function currently in beta phase. Waiting for acceptance of peer-reviewed paper
#' @examples
#'
#' var1<-vardata[,10]
#' subs<-nrow(varmat)
#' voxels<-ncol(varmat)
#' regpval<-matrix(nrow=1,ncol=voxels)
#' regtstat<-matrix(nrow=1,ncol=voxels)
#' resmat<-matrix(OL,nrow=subs,ncol=voxels)
#' for (i in 1:voxels){
#' vox<-varmat[,i]
#' regfit<-lm(vox~var1)
#' ##Extract statistical values
#' resmat[,i]<-residuals(regfit)
#' regsum<-summary(regfit)
#' regtstat[,i]<-regsum$coefficients[3,3]
#' }
#' fwhm<-estPooled.smooth(res,rdf,mask)
#'
#'
#' @export rft.thresh
rft.thresh <-function(StatImg, pval, k, fwhm, resels, df, fieldType, threshType){
  if (missing(threshType)){
    stop("Must specify threshold type")
  }
  if (missing(pval)){
    stop("Must specify pval")
  }
  D <-mask@dimension
  if (threshType=="cluster" | threshType=="voxel"){
    if (missing(fieldType)){
      stop("Must specify fieldType if computing voxel or cluster level threshold")
    }
    alpha <-pval-1
    u <-15
    Mfwhm <-mean(fwhm)
    while(alpha < pval){
      stat <-u-.01
      if (threshType=="cluster"){
        alpha <-rft.pval(D, 1, k, u, resels, df, fieldType)
      }else if(threshType=="voxel"){
        alpha <-rft.pval(1, 1, maxpeak, resels, df, fieldType)
      }else{
        cat("Must specify appropriate threshType \n")
      }
    }
  thresh <-u
  }else if (threshType=="cfdr" | threshType=="vfdr"){
    if (threshType=="cfdr"){
      if (fieldType=="Z"){
        stat <- qnorm(1 - pval)
      }else if(fieldType=="T"){
        stat <- qt(1 - pval, df = df[2])
      }else if(fieldType=="F"){
        stat <- qf(1 - pval, df1 = df[1], df2 = df[2])
      }else if(fieldType=="X"){
        stat <-qchisq(1-pval, df[1],df[2])
      }
      statimg <-image2ClusterImages(StatImg, minClusterSize=1,minThresh=stat,maxThresh=Inf)
      cmax <-c()
      for (i in 1:length(clist)){
        cmax <-cbind(cmax,max(clist[[i]]))
        }
      }else if (threshType=="vfdr"){
        cmax <-as.array(StatImg)
      }
      if (fieldType=="Z"){
        p <-sort(1 - pnorm(cmax),decreasing=TRUE)
      }else if(fieldType=="T"){
        p <-sort(1 - pt(cmax, df = df[1]),decreasing=TRUE)
      }else if(fieldType=="F"){
        p <-sort(1 - pf(cmax, df1 = df[1], df2 = df[2]),decreasing=TRUE)
      }else if(fieldType=="X"){
        p <-sort(1-pchisq(cmax, df[1],df[2]),decreasing=TRUE)
      }
      pfdr <-sort((alpha*(1:length(p))/length(p)),decreasing=TRUE)
      i <- 1
      while (p[i] >= pfdr[i]){
        i <-i+1
      }
      if (fieldType=="Z"){
        thresh <- qnorm(1 - p[i])
      }else if(fieldType=="T"){
        thresh <- qt(1 - p[i], df = df[1])
      }else if(fieldType=="F"){
        thresh <- qf(1 - p[i], df1 = df[1], df2 = df[2])
      }else if(fieldType=="X"){
        thresh <-qchisq(1-p[i], df[1],df[2])
      }
  }
  thresh
  }
