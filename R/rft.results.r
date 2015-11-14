#' @name rft.results
#' @title Utilizes RFT to produce cluster and voxel level statistics
#' 
#' 
#' @param img-SPM of class antsImage 
#' @param stat-statistical threshold for SPM
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
#' @description
#'
#'	
#' @reference 
#'	Friston K.J., (1994) Assessing the Significance of Focal Activations Using Their Spatial Extent
#'	Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power
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
#'	
#'
#' @export rft.results
rft.results<-function(img,stat,ka,fwhm,mask,df,fieldType){
	voxels <-sum(as.array(mask))
	D<-mask@dimension
	Mfwhm<-mean(fwhm)
	Clusters<-list()
	posclust <- labelClusters(img, ka, stat, Inf)
	if(max(posclust) < 1){
  		cat("No positive clusters survive threshold
  		")
  	}else{
		labs <- unique(posclust[posclust > 0])
		posclustlist <- list()
		for (i in 1:length(labs)) {
			labimg <- antsImageClone(img)
			labimg[posclust != labs[i]] <- 0
			posclustlist <- lappend(posclustlist, labimg)
			Clusters<-lappend(Clusters, labimg)
			clustername<-paste("PCluster",i,sep="")
			names(Clusters)[length(Clusters)]<-clustername
			}

		postable<-matrix(nrow=length(posclustlist),ncol=7)
		colnames(postable)<-c("Voxels", "Cluster-Probability", "Voxel-Probability","Peak-Height","xc", "yc", "zc")
		cnames<-rep(1,nrow=5)
		for (i in 1:length(posclustlist)){
			cat("Determing positive cluster-level statistics:",i,"
			",sep=" ")
			cMask<-getMask(posclustlist[[i]])
			cvoxs<-sum(as.array(cMask))
			pclust<-rft.pcluster(cMask,mask,Mfwhm,stat,df,fieldType)
			loc<-labelImageCentroids(cMask)[2]
			cat("Determing positive voxel-level statistics:",i,"
			",sep=" ")
			resel<-ants.resel(cMask,fwhm)
			ec<-ants.ec(stat,fieldType,df)
			pvox<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
			postable[i,]<-c(cvoxs,pclust,pvox,max(posclustlist[[i]]),loc$vertices[1],loc$vertices[2],loc$vertices[3])
			clustername<-paste("P-Cluster:",i,sep="")
			cnames[i]<-clustername
			}
		rownames(postable)<-cnames
		Clusters<-lappend(Clusters,postable)
		names(Clusters)[length(Clusters)]<-"PositiveStatistics"
	}
	
	nimg<-img*-1
	negclust <- labelClusters(nimg, ka, stat, Inf)
	if(max(negclust) < 1){
  		cat("No negative clusters survive threshold
  		")
	}else{
		labs <- unique(negclust[negclust > 0])
		negclustlist <- list()
		for (i in 1:length(labs)) {
			labimg <- antsImageClone(nimg)
			labimg[negclust != labs[i]] <- 0
			negclustlist <- lappend(negclustlist, labimg)
			Clusters<-lappend(Clusters,labimg)
			clustername<-paste("NCluster",i,sep="")
			names(Clusters)[length(Clusters)]<-clustername
			}
	
		negtable<-matrix(nrow=length(negclustlist),ncol=7)
		colnames(negtable)<-c("Voxels", "Cluster-Probability", "Voxel-Probability","Peak-Height","xc", "yc", "zc")
		cnames<-rep(1,nrow=5)
		for (i in 1:length(negclustlist)){
			cat("Determing negative cluster-level statistics:",i,"
			",sep=" ")
			cMask<-getMask(negclustlist[[i]])
			cvoxs<-sum(as.array(cMask))
			pclust<-rft.pcluster(cMask,mask,Mfwhm,stat,df,fieldType)
			loc<-labelImageCentroids(cMask)[2]
			cat("Determing negative voxel-level statistics:",i,"
			",sep=" ")
			resel<-ants.resel(cMask,fwhm)
			ec<-ants.ec(stat,fieldType,df)
			pvox<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
			negtable[i,]<-c(cvoxs,pclust,pvox,max(negclustlist[[i]]),loc$vertices[1],loc$vertices[2],loc$vertices[3])
			clustername<-paste("N-Cluster:",i,sep="")
			cnames[i]<-clustername
			}
		rownames(negtable)<-cnames
		Clusters<-lappend(Clusters,negtable)
		names(Clusters)[length(Clusters)]<-"NegativeStatistics"
	}
	
	return(Clusters)
	}
