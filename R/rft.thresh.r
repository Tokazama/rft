#' @name rft.thresh
#' @title Calculates 
#' A quick heuristic for thresholding a Statistical field using RFT
#'
#'
#' @param pval-probability of false positive
#' @param ka-minimum desired cluster size
#' @param fwhm-full width at half maxima
#' @param mask-antsImage mask
#' @param df-degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:
#'	"T"- T-field
#'	"F"- F-field
#'	"X"- Chi squar field
#'	"Z"- Gaussian field
#' @return Outputs a statistical value to be used for threshold a SPM
#' @reference Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power
#' @
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
#' @export rft.thresh
rft.thresh<-function(img,pval,ka,fwhm,mask,df,fieldType){
	voxels <-sum(as.array(mask))
	bMask <-mask
	cMask <- ka
	D<-length(dim(mask))
	fwhm<-mean(fwhm)
	alpha<-pval-1
	stat<-10
	
	cat("Determing threshold value based on pval, ka, brain volume")
	while(alpha < pval){
		stat <-stat-.01
		alpha <-rft.pcluster(cMask,bMask,fwhm,stat,df,fieldType)
		}

	posclust<-image2ClusterImages(img,150,stat,Inf)
	lclust<-labelClusters(
	postable<-matrix(ncol=5)
	colnames(postable)<-c("Voxels", "Cluster-Probability", "xc", "yc", "zc")
	for (i in 1:length(posclust)){
		cat("Determing positive cluster level statistics:",i,sep=" ")
		cMask<-getMask(posclust[[i]])
		cvoxs<-sum(as.array(cMask))
		pclust<-rft.pcluster(cMask,mask,fwhm,thresh,df,fieldType)
		peak<-max(posclust[[i]])
		loc<-labelImageCentroids(cMask)[2]
		postable[i,]<-c(cvox, pclust, loc$vertices[1],loc$vertices[2],loc$vertices[3])
		clustername<-paste("P-Cluster:",i,sep="")
		rownames(postable[i,])<-c(clustername)
		image<-paste("Plcuster",i,".nii.gz",sep="")
		}
	
	negclust<-image2ClusterImages(img,150,stat,Inf)
	negtable<-matrix(ncol=5)
	colnames(postable)<-c("Voxels", "Cluster-Probability", "xc", "yc", "zc")
	for (i in 1:length(negclust)){
		cat("Determing positive cluster level statistics:",i,sep=" ")
		cMask<-getMask(negclust[[i]])
		cvoxs<-sum(as.array(cMask))
		pclust<-rft.pcluster(cMask,mask,fwhm,thresh,df,fieldType)
		peak<-max(negclust[[i]])
		loc<-labelImageCentroids(cMask)[2]
		negtable[i,]<-c(cvox, pclust, loc$vertices[1],loc$vertices[2],loc$vertices[3])
		clustername<-paste("P-Cluster:",i,sep="")
		rownames(postable[i,])<-c(clustername)
		image<-paste("Plcuster",i,".nii.gz",sep="")
		}
	cluster.table<-rbind(postable,negtable)
	
	return(list(cluster.table,posclust,negclust,stat))
	}
