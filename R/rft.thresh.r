#' @name rft.thresh
#' @title Thesholds SPM producing clusters and statistics
#' 
#'
#' @param img-SPM of class antsImage 
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
#' @description
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
#' @export rft.thresh
rft.thresh<-function(img,pval,ka,cthresh,fwhm,mask,df,fieldType){
	voxels <-sum(as.array(mask))
	bMask <-mask
	cMask <- ka
	D<-mask@dimension
	Mfwhm<-mean(fwhm)
	alpha<-pval-1
	stat<-10
	Clusters<-list()
	cat("Determing threshold value based on pval, ka, brain volume
	")
	while(alpha < pval){
		stat <-stat-.01
		alpha <-rft.pcluster(cMask,bMask,Mfwhm,stat,df,fieldType)
		}
	Clusters<-lappend(Clusters,stat)
	names(Clusters)[1]<-"Threshold Value"
	
	posclust <- labelClusters(img, cthresh, stat, Inf)
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
			}

		Clusters<-lappend(Clusters,posclustlist)
		names(Clusters)[length(Clusters)]<-"PositiveClusters"
		postable<-matrix(nrow=length(posclustlist),ncol=7)
		colnames(postable)<-c("Voxels", "Cluster-Probability", "Voxel-Probability","Peak-Height","xc", "yc", "zc")
		cnames<-rep(1,nrow=5)
		for (i in 1:length(posclustlist)){
			cat("Determing positive cluster level statistics:",i,"
			",sep=" ")
			cMask<-getMask(posclustlist[[i]])
			cvoxs<-sum(as.array(cMask))
			pclust<-rft.pcluster(cMask,mask,Mfwhm,stat,df,fieldType)
			loc<-labelImageCentroids(cMask)[2]
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
	
	negclust <- labelClusters(img, cthresh, -Inf, -stat)
	if(max(negclust) < 1){
  		cat("No negative clusters survive threshold
  		")
	}else{
		labs <- unique(negclust[negclust > 0])
		negclustlist <- list()
		for (i in 1:length(labs)) {
			labimg <- antsImageClone(img)
			labimg[negclust != labs[i]] <- 0
			negclustlist <- lappend(negclustlist, labimg)
			}
	
		Clusters<-lappend(Clusters,negclustlist)
		names(Clusters)[length(Clusters)]<-"NegativeClusters"
		negtable<-matrix(nrow=length(negclustlist),ncol=7)
		colnames(postable)<-c("Voxels", "Cluster-Probability", "Voxel-Probability","Peak-Height","xc", "yc", "zc")
		for (i in 1:length(negclustlist)){
			cat("Determing positive cluster level statistics:",i,"
			",sep=" ")
			cMask<-getMask(negclustlist[[i]])
			cvoxs<-sum(as.array(cMask))
			pclust<-rft.pcluster(cMask,mask,Mfwhm,stat,df,fieldType)
			loc<-labelImageCentroids(cMask)[2]
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
