#' Utilizes RFT to produce cluster and voxel level statistics
#' 
#' @param D image dimensions
#' @param thresh statistical threshold for SPM
#' @param ka minimum desired cluster size
#' @param mask antsImage mask
#' @param fwhm full width at half maxima
#' @param resel resel values for the mask
#' @param StatImg statistical field image of class antsImage 
#' @param df degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:
#' \describe{
#' \item{"T"}{T-field} 
#' \item{"F"}{F-field} 
#' \item{"X"}{Chi-square field"} 
#' \item{"Z"}{Gaussian field}
#' }
#' @param cplot create plot of clusters
#' @param statdir directory where output is saved (by default "./")
#' @return Outputs a statistical value to be used for threshold a statistical field image
#' \describe{
#' \item{"Pset"}{"Probability of obtaining this set of clusters"}
#' \item{"clusters"}{"The number of clusters in the set"}
#' \item{"cPfwe"}{"FWE corrected p value of cluster"}
#' \item{"cPfdr"}{"FDR corrected p value of cluster"}
#' \item{"cke"}{"Number of voxels in the cluster"}
#' \item{"cPuncorr"}{"Uncorrected p value of cluster"}
#' \item{"pPfwe"}{"FWE correct p value of cluster peak maxima"}
#' \item{"pPfdr"}{"FDR correct p value of cluster peak maxima"}
#' \item{"pT"}{"t-score of cluster peak maxima"}
#' \item{"pZ"}{"Z-score of cluster peak maxima"}
#' \item{"pPuncorr"}{"Uncorrected p value of cluster peak maxima"}
#' \item{"xc"}{"x-coordinate of cluster centroid"}
#' \item{"yc"}{"y-coordinate of cluster centroid"}
#' \item{"zc"}{"z-coordinate of cluster centroid"}
#' }
#' 
#' Height threshold: the p-level for uncorrected alpha
#' Extent threshold: The minimum number of adjacent suprathreshold voxels
#' Degrees of freedom
#' FWHM
#' Resels (resolution elements)
#' Search volume
#' 
#' 
#' @description
#' \code{rft.results} 
#'
#'	
#' @References
#' Chumbley J., (2010) Topological FDR for neuroimaging
#' Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power
#' Worsley K.J., (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' 
#' @Author Zachary P. Christensen
#' @kewords rft.pcluster, ants.ec
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
#'  vox<-varmat[,i]
#'  regfit<-lm(vox~var1)
#'  resmat[,i]<-residuals(regfit)
#'  regsum<-summary(regfit)
#'	 regtstat[,i]<-regsum$coefficients[3,3]
#'	 }
#' rdf<-regfit$df.residual
#' timg<-makeImage(mask,regtstat)
#' fwhm<-estPooled.smooth(res,rdf,mask)
#' thresh<-rft.thresh(D,img,pval,ka,fwhm,mask,rdf,"T")
#' results(D,thresh,ka,fwhm,timg,mask,rdf,"T")
#'	
#' @export rft.results
rft.results<-function(StatImg, k, u, resels, df, fieldType, mask, cplot=FALSE, cimg=FALSE, statdir="./"){
  
  D <-StatImg@dimensions
  clustimg <-labelClusters(StatImg, minClusterSize=ka, minThresh=thresh, maxThresh=Inf)
  antsImageWrite(clustimg, file=paste(statdir,"cluster_image.nii.gz",sep=""))
  
  labs <-unique(clustimg[clustimg > 0])
  clusters <-list()
  xcol <-cbind()
  ycol <-cbind()
  zcol <-cbind()
  pPcol <-cbind()
  statcol <-cbind()
  pZcol <-cbind()
  pPfdrcol <-cbind()
  cVolcol <-cbind()
  cPfdrcol <-cbind()
  clust.names <-cbind()
  Ecu <-rft.pval(1, 1, u, resels, df, fieldType)$Ec
  
  for (i in 1:length(labs)){
    labimg <-antsImageClone(StatImg)
    labimg[clustimg !=labs[i]] <-0
    clusters <-lappend(clusters, labimg)
    
    clust.names <-cbind(clust.names, paste("Cluster",i,sep=""))
    cmask <-getMask(clusters[[i]])
    loc <-getCentroids(cmask)
    xcol <-cbind(xcol, loc[1])
    ycol <-cbind(ycol, loc[2])
    zcol <-cbind(zcol, loc[3])
    
    # cluster-level
    k <-rft.resels(cmask,fwhm)[4]
    cVolcol <-cbind(cVolcol,k)
    cpval <-rft.pval(1, k, u, resels, df, fieldType)
    cPfwecol <-cbind(cPfdrcol,cpval$Pcor)
    cPcol <-cbind(cPcol, cpval$Pu)
    
    # voxel-level
    maxpeak <-max(clusters[[i]])
    statcol <-cbind(statcol, maxpeak)
    ppval <-rft.pval(1, 1, maxpeak, resels, df, fieldType)
    punc <-ppval$Ec/Ecu
    pPcol <-cbind(pPcol, punc)
    pZcol <-cbind(pZcol, 1-qnorm(punc))
  }
  # FDR correction of uncorrected values
  pPfdrcol <-p.adjust(pPcol,"BH")
  cPfdrcol <-p.adjust(pPfdrcol,"BH")
  
  # create table
  nclusts <-length(clusters)
  clustable <-matrix(nrow=nclusts, ncol=14)
  colnames(clustable) <-c("Pset","clusters","cPfwe", "cPfdr", "cVol", "cPu", "pPfwe", "pPfdr", "pT", "pZ", "pPu", "xc", "yc", "zc")
  
  clustable[1,1] <-rft.pval(length(labs), min(cVolcol), thresh, resels, df, fieldType)$Pcor
  clustable[1,2] <-length(labs)
  clustable[,3] <-clust.names
  clustable[,4] <-cPfdrcol
  clustable[,5] <-cVolcol
  clustable[,6] <-cPcol
  clustable[,7] <-pPfwecol
  clustable[,8] <-pPfdrcol
  clustable[,9] <-statcol
  clustable[,10] <-pZcol
  clustable[,11] <-pPcol
  clustable[,12] <-xcol
  clustable[,13] <-ycol
  clustable[,14] <-zcol
  write.csv(clustable, file=paste(statdir,"ClusterStats.csv",sep=""))
  if (cplot=="TRUE"){
    
  }
  if (cimg=="TRUE"){
    brain <-renderSurfaceFunction(surfimg=list(mask), alphasurf=.1, funcimg=clustimg, smoothsval=1.5, smoothfval=0)
    id <-par3d("userMatrix")
    sag <-rotate3d(id,-pi/2,1,0,0)
    axi <-rotate3d(id,pi/2,0,0,1)
    cor <-rotate3d(id,-pi/2,0,0,1)
    par3d(userMatrix=id)
    dd <-make3ViewPNG(sag, cor, axi, paste(statdir,"Clusters.png",sep=''))
    par3d(userMatrix=id)
  }
  write.csv(clustable, file=paste(statdir, "RFT_Results.csv", sep=""))
  z <-list(ClusterImg=clustimg, ClusterStats=clustable)
  z
}
