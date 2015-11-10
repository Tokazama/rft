#' @name rft.pcluster
#' @title Calculates probability of 
#' A quick heuristic for thresholding a Statistical field using RFT
#'
#'
#' @param cMask-antsImage mask of the cluster
#' @param bMask-antsImage mask of the brain
#' @param fwhm-full width at half maxima
#' @param stat-statistical value that was used to threshold the image
#' @param df-degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:
#'	"T"- T-field
#'	"F"- F-field
#'	"X"- Chi squar field
#'	"Z"- Gaussian field
#' @return Outputs a statistical value to be used for threshold a SPM
#' @reference Friston K.J., (1994) Assessing the Significance of Focal Activations Using Their Spatial Extent.
#' @reference Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power.
#'
#' @examples
#'
#'  var1<-vardata[,10]
#'  subs<-nrow(varmat)
#'  voxels<-ncol(varmat)
#'  regpval<-matrix(nrow=1,ncol=voxels)
#'  regtstat<-matrix(nrow=1,ncol=voxels)
#'  resmat<-matrix(OL,nrow=subs,ncol=voxels)
#'  for (i in 1:voxels){
#'	  vox<-varmat[,i]
#'	  regfit<-lm(vox~var1)
#'	  ##Extract statistical values
#'	  resmat[,i]<-residuals(regfit)
#'	  regsum<-summary(regfit)
#'	  regtstat[,i]<-regsum$coefficients[3,3]
#'	  }
#'  fwhm<-estPooled.smooth(res,rdf,mask)
#'	negclust<-image2ClusterImages(timg,150,-Inf,-thresh)
#'	negtable<-matrix(ncol=5)
#'	colnames(postable)<-c("Voxels", "Cluster-Probability", "Peak-Height", "Voxel-Probability", "Coordinates")
#'	for (i in 1:length(posclust)){
#'	  cat("Determing negative cluster level statistics:",i,sep=" ")
#'	  cmask<-getMask(negclust[[i]])
#'	  cvoxs<-sum(as.array(cmask))
#'	  pclust<-rft.pcluster(negclust[[i]],mask,fwhm,thresh,df,fieldType)
#'	  peak<-max(posclust[[i]])
#'	  loc<-labelImageCentroids(posclust[[i]])[2]
#'	  resel <-ants.resel(mask,fwhm)
#'	  ec<-ants.ec(stat,fieldType,df)
#'	  pval<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
#'	  negtable[i,]<-c(cvox, pclust, peak, pval, loc$vertices[1],loc$vertices[2],loc$vertices[3])
#'	  clustername<-paste("N-Cluster:",i,sep="")
#'	  rownames(postable[i,])<-c(clustername)
#'	  image<-paste(fileDir,"Nlcuster",i,".nii.gz",sep="")
#'	  antsImageWrite(negclust[[i]],file=image)
#'	  }
#'
#' @export rft.pcluster
rft.pcluster<-function(cMask,bMask,fwhm,stat,df,fieldType){
	voxels<-sum(as.array(mask))
	D<-length(dim(mask))
	if(class(cMask)=="numeric"){
		ka<-cMask
	}else{
		ka<-sum(as.array(mask))
		}
	fwhm<-mean(fwhm)
	
	if (fieldType=="T"){
		EN<-voxels*(1-pt(stat,df))
	}else if(fieldType=="F"){
		EN<-voxels*(1-pf(stat, df[1],df[2]))
	}else if(fieldType=="X"){
		EN<-voxels*(1-pchisq(stat, df[1],df[2]))
	}else if(fieldType=="G"){
		EN<-voxels*(1-qnorm(stat))
	}else{
		stop("A correct fieldtype is required")
		}
	#Based upon the fwhm values or the covariance matrix 
	#of a fields partial derivatives(friston 1994)
	W<-fwhm/((4*log(2))^(1/2))
	
	#The expected number of maxima/clusters (friston 1996; equation 2)
	Em<-bvox*(2*pi)^(-(D+1)/2)*(W^(-D))*(stat^(D-1))*(exp((-stat^2)/2))
	
	#Used to derive the number of expected clusters in a 
	#voxel (friston 1996;equation 3)
	rfB<-(gamma((D/2)+1)*Em/EN)^(2/D)
	
	#Equation 4 (friston 1996) The probability of getting a
	#cluster of the size 'ka' or more above a threshold 'stat'
	#for a single cluster
	Pnk<-exp(-rfB*(ka^(2/D)))
	alpha<-1-exp(-Em*Pnk)
	return(alpha)
	}
