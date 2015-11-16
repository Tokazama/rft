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
rft.results<-function(StatImg,stat,ka,fwhm,df,fieldType){
	mask<-getMask(StatImg)
	voxels <-sum(as.array(mask))
	D<-mask@dimension
	Mfwhm<-mean(fwhm)
	absimg<-as.antsImage(abs(as.array(StatImg)))
	clusters<-image2ClusterImages(absimg,ka,stat,Inf)
	nclusts<-length(clusters)
	clustable<-matrix(nrow=nclusts,ncol=7)
	colnames(clustable)<-c("Voxels", "Cluster-Probability", "Voxel-Probability","Peak-Height","xc", "yc", "zc")
	clist<-list()
	cnames<-rep(1,nrow=nclusts)
	for (i in 1:nclusts){
		cat("Determing cluster-level statistics:",i,"
		",sep=" ")
		cName<-paste("Cluster",i,sep="")
		cMask<-getMask(clusters[[i]])
		cvoxs<-sum(as.array(cMask))
		loc<-labelImageCentroids(cMask)[2]
		clust<-maskImage(StatImg,cMask,level=1)
		clist<-lappend(clist,clust)
		names(clusters)[[i]]<-cName
		names(clist)[length(clist)]<-cName
		pclust<-rft.pcluster(cMask,mask,Mfwhm,stat,df,fieldType)
		cat("Determing voxel-level statistics:",i,"
		",sep=" ")
		resel<-ants.resel(cMask,fwhm)
		ec<-ants.ec(stat,fieldType,df)
		pvox<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
		clustable[i,]<-c(cvoxs,pclust,pvox,max(clust),loc$vertices[1],loc$vertices[2],loc$vertices[3])
		cnames[i]<-cName
		}
	rownames(clustable)<-cnames
	clist<-lappend(clist,clustable)
	names(clist)[length(clist)]<-"stats"
	return(clist)
	}
