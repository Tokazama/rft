#' RFT Statistical Results
#'
#' Returns RFT based statistical results for a single statistical image
#'
#' @param x statistical field image of class antsImage 
#' @param resels resel values for the mask
#' @param fwhm full width at half maxima
#' @param df degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:
#' \itemize{
#' \item{"T"}{T-field}
#' \item{"F"}{F-field}
#' \item{"X"}{Chi-square field"}
#' \item{"Z"}{Gaussian field}
#' }
#' @param RPVImg resels per voxel image
#' @param k minimum desired cluster size
#' @param thresh a numeric value to threshold the statistical field or a character of the following methods:
#' \itemize{
#'	\item{"cRFT"}{computes a threshold per expected cluster level probability }
#'	\item{"pRFT"}{uses the mask and pval calculates the minimum statistical threshold}
#'	\item{"cFDR"}{uses an uncorrected threshold at the alpha level and then computes and FDR threshold based on cluster maxima}
#'	\item{"pFDR"}{computes the fdr threshold for the entire field of voxels}
#' }
#' @param pval the p-value for estimating the threshold
#' @param pp the initial p-value for thresholding (only used for FDR methods)
#' @param n number of images in conjunction
#' @param tex output statistical table in latex format (uses package \code{xtable})
#' @param statdir directory where output is saved (if not specified images are not saved)
#' @param verbose enables verbose output
#'
#' @return Outputs a statistical value to be used for threshold a statistical field image
#' \itemize{
#' \item{"SetStats"}{"set-level statistics and number of clusters"}
#' \item{"ClusterStats"}{"cluster-level statistics and descriptors"}
#' \item{"PeakStats"}{"peak-level statistics and descriptors"}
#' \item{"LabeledClusters"}{"image of labeled clusters"}
#' \item{"u"}{"the threshold used"}
#' }
#' 
#' @description 
#' 
#' \code{rft.pval} is used to compute all family-wise error (FWE) corrected statistics while \code{p.adjust} is used to compute all false-discovery rate
#' based statistics. All statistics herein involve implementation of random field theory (RFT) to some degree.
#' 
#' Both cluster-level and peak-level statistics are described by the uncorrected p-value along with the FDR and FWE corrected p-values for each cluster. 
#' Peak-level statistics are described by the maximum statistical value in each cluster and the comparable Z statistic. The ClusterStats table also contains
#' coordinates for each cluster and the number of voxels therein. By default thresh="pRFT" and pval=.05. Alternatively, the user may use a specific numeric 
#' value for thresholding the statistical field. \code{rft.thresh} more fully describes using appropriate thresholds for statistical fields.
#' 
#' 
#' @References
#' Chumbley J., (2010) Topological FDR for neuroimaging
#' Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power
#' Worsley K.J., (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' 
#' @Author Zachary P. Christensen
#' @kewords rft.pval
#' @note: function currently in beta phase
#' @examples
#' 
#' # estimatation of a single images smoothness
#' outimg1 <-makeImage(c(10,10,10), rt(1000))
#' outimg2 <-makeImage(c(10,10,10), rt(1000))
#' mask <-getMask(outimg1)
#' imat <-imageListToMatrix(list(outimg1, outimg2), mask)
#' variable <-rnorm(2)
#' 
#' # fit the model
#' residuals <-matrix(nrow=2,ncol=ncol(imat))
#' tvals <-matrix(nrow=1,ncol=ncol(imat))
#' for (i in 1:ncol(imat)){
#'  fit <-lm(variable~imat[,i])
#'  residuals[,i] <-residuals(fit)
#'  tvals[,i] <-summary(fit)$coefficients[2,3]
#' }
#' 
#' df <-c(1,fit$df.residuals)
#' # scale residuals
#' mrss <-colMeans(res^2)/df
#' sr <-r/sqrt(mrss) #scaled residuals
#' fwhm <-estSmooth(sr,mask,df)
#' resels <-rft.resels(mask,fwhm$fwhm)
#' 
#' # specifying thresholds
#' timg <-makeImage(mask,tvals)
#' results <-rft.results(timg, resels, fwhm$fwhm, df, fieldType="T", thresh="pRFT", pval=.05) # threshold to create peak values with p-value of .05 (default)
#' results <-rft.results(timg, resels, fwhm$fwhm, df, fieldType="T", thresh="cRFT", pval=.05) # threshold to create clusters with p-value of .05
#' results <-rft.results(timg, resels, fwhm$fwhm, df, fieldType="T", thresh="pFDR", pval=.05, pp=.01) # initial threshold at p-value of .001 
#'                                                                                            followed by peak FDR threshold at p-value of .05
#' results <-rft.results(timg, resels,fwhm,df,fieldType="T",thresh="cFDR", pval=.05, pp=.01) # initial threshold at p-value of .001 
#'                                                                                            followed by cluster FDR threshold at p-value of .05
#' # correcting for non-isotropic
#' results <-rft.results(timg, resels, fwhm$fwhm, df, fieldType="T", fwhm$RPVImg)
#'
#'
#' @export rft.results
rft.results <-function(x, resels, fwhm, df, fieldType, RPVImg, k=1, thresh="pRFT", pval=.05, pp=.01, n=1, tex=F, statdir=F, verbose=FALSE){
  if (missing(x)){
    stop("Must specify x")
  }
  if (missing(resels)){
    stop("Must specify resels")
  }
  if (missing(fwhm)){
    stop("Must specify fwhm")
  }
  if (missing(df)){
    stop("Must specify degrees of freedom (df)")
  }
  if (missing(fieldType)){
    stop("Must specify fieldType")
  }

  if (class(thresh)=="character"){
	if (verbose=="TRUE"){
	  cat("Calculating threshold \n")
	}
  threshType=thresh
  u <-rft.thresh(x, pval, k, n, fwhm, resels, df, fieldType, threshType, pp,verbose=verbose)
  }else if(class(thresh)=="numeric"){
    u <-thresh
  }else{
    stop("thresh must be a numeric value or a list of specifications for calculating a threshold")
  }
  D <-x@dimension
  vox2res <-1/prod(fwhm)
  k <-k*vox2res
  nvox <-length(as.vector(x[x !=0]))
  
  mask <-getMask(x)
  # extract clusters at threshold
  clust <-labelClusters(x,k,u,Inf)
  if (!(missing(statdir))){
    antsImageWrite(clust,file=paste(statdir))
  }
  labs <-unique(clust[clust > 0])
  nclus <-length(labs) # number of clusters
  stats <-matrix(nrow=nclus, ncol=12)
  stats[1:nclus,10:12] <-getCentroids(clust)[1:nclus,1:3]

  # set-level stat
  Pset <-rft.pval(D, nclus, k, u, n, resels, df, fieldType)$Pcor
  stats <-as.table(stats)
  Ez <-rft.pval(D, 1, 0, u, n, resels, df, fieldType)$Ec
  EC <-c()
  for (i in 1:nclus){
    cname <-paste("Cluster",i,sep="")
    cat("Calculating statistics for",cname," \n",sep=" ")
    rownames(stats)[i] <-cname
    cmask <-antsImageClone(mask)
    cmask[clust !=labs[i]] <-0
    stats[i,4] <-sum(as.array(cmask)) # number of voxels in cluster

    if (missing(RPVImg)){
      # follows isotropic image assumptions
      K <-stats[i,4]*vox2res
    }else{
      # extract resels per voxel in cluster (for non-isotropic image)
      rkc <-RPVImg[cmask==1]
      lkc <-sum(rkc)/stats[i,4]
      iv <-matrix(rft.resels(cmask,c(1,1,1)),nrow=1)
      iv <-iv %*% matrix(c(1/2,2/3,2/3,1),ncol=1)
      K <-iv*lkc
    }
    # max cluster value
    stats[i,8] <-max(x[cmask==1])
	
	# uncorrected p-value (peak)
    if (fieldType=="Z"){
      stats[i,7] <-1-pnorm(stats[i,8])
    }else if(fieldType=="T"){
      stats[i,7] <-1-pt(stats[i,8], df[2])
    }else if(fieldType=="F"){
      stats[i,7] <-1-pf(stats[i,8], df[1], df[2])
    }else if(fieldType=="X"){
      stats[i,7] <-1-pchisq(stats[i,8], df[1],df[2])
    }
    # cluster/peak-level stats
    ppeak <-rft.pval(D, 1, 0, stats[i,8], n, resels, df, fieldType) 
    stats[i,5] <-ppeak$Pcor # fwe p-value (peak)

    pclust <-rft.pval(D, 1, K, u, n, resels, df, fieldType)
    stats[i,1] <-pclust$Pcor # fwe p-value (cluster)
    stats[i,3] <-pclust$Punc
    
	EC <-cbind(EC,rft.pval(D, 1, 0, stats[i,8], n, resels, df, fieldType)$Ec/Ez) 
    # comparable Z stat
    stats[i,9] <-qnorm(1-stats[i,7])
  }
  # FDR correction (see Chumbley J., (2010) Topological FDR for neuroimaging)
  stats[,2] <-p.adjust(stats[,3],"BH") # cluster FDR
  stats[,6] <-p.adjust(EC,"BH") # peak FDR
  colnames(stats) <-c("Pfwe", "Pfdr", "P", "Voxels",
                      "Pfwe", "Pfdr", "P", "MaxStat",
                      "Z", "xc", "yc", "zc")
  setstats <-round(as.table(matrix(c(Pset,nclus),nrow=1,ncol=2)),4)
  colnames(setstats) <-c("p-value","Clusters")
  rownames(setstats) <-c("")
  
  cluststats <-round(stats[1:nclus,c(1:4,10:12)],4)
  peakstats <-round(stats[1:nclus,5:9],4)
  results <-list(setstats, cluststats, peakstats, clust, u)
  names(results) <-c("SetStats", "ClusterStats", "PeakStats", "LabeledClusters", "u")
  if (tex=="T"){
    texTables <-list(xtable(results$SetStats),xtable(results$ClusterStats),xtable(results$PeakStats))
    names(texTables) <-c("SetStats", "ClusterStats", "PeakStats")
    results <-lappend(results,texTables)
    names(results)[length(results)] <-"texTables"
  }
  if (statdir=="T"){
    write.csv(stats, file=paste(statdir))
  }
  results
}
