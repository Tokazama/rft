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
rft.results<-function(D,thresh,ka,fwhm,StatImg,mask,df,fieldType){
	cat("Calculating image resels
	",sep="")
	resel<-ants.resel(mask,fwhm)
	voxels <-sum(as.array(mask))
	Mfwhm<-mean(fwhm)
	absimg<-as.antsImage(abs(as.array(StatImg)))
	clist<-list()
	if(thresh > max(absimg)){
		cat("No clusters survive statistical threshold
		")
	}else{
		clusters<-image2ClusterImages(absimg,0,thresh,Inf)
		for (i in 1:length(clusters)){
			cMask<-getMask(clusters[[i]])
			if(sum(as.array(cMask)) > ka){
				clist<-lappend(clist,cMask)
				}
			}
		if(length(clist) < 1){
			cat("No clusters survive ka threshold
			")
		}else{
			nclusts<-length(clusters)
			clustable<-matrix(nrow=nclusts,ncol=7)
			colnames(clustable)<-c("Voxels", "Cluster-Probability", "Voxel-Probability","Peak-Height","xc", "yc", "zc")
			cnames<-rep(1,nrow=nclusts)
			for (i in 1:nclusts){
				cat("Determing cluster-level statistics:",i,"
				",sep=" ")
				cName<-paste("Cluster",i,sep="")
				cMask<-clist[[i]]
				cvoxs<-sum(as.array(cMask))
				loc<-labelImageCentroids(cMask)[2]
				clust<-maskImage(StatImg,cMask,level=1)
				clist<-lappend(clist,clust)
				names(clusters)[[i]]<-cName
				names(clist)[length(clist)]<-cName
				pclust<-rft.pcluster(cMask,mask,Mfwhm,thresh,df,fieldType)
				cat("Determing voxel-level statistics:",i,"
				",sep=" ")
				maxpeak<-max(clust)
				minpeak<-abs(min(clust))
				if (maxpeak > minpeak){
					peak<-max(clust)
					ec<-ants.ec(peak,fieldType,df)
				}else{
					peak<-min(clust)
					ec<-ants.ec(minpeak,fieldType,df)
					}
				pvox<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
				clustable[nrow(clustable)+1,]<-c(cvoxs,pclust,pvox,peak,loc$vertices[1],loc$vertices[2],loc$vertices[3])
				cnames[i]<-cName
				}
			rownames(clustable)<-cnames
			clist<-lappend(clist,clustable)
			names(clist)[length(clist)]<-"stats"
			}
		}
	clist<-lappend(clist,resel)
	names(clist)[length(clist)]<-"resels"
	clist<-lappend(clist,fwhm)
	names(clist)[length(clist)]<-"fwhm"
	return(clist)
	}
