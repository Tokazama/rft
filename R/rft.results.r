#' Utilizes RFT to produce cluster and voxel level statistics
#' 
#' @param D image dimensions
#' @param thresh statistical threshold for SPM
#' @param ka minimum desired cluster size
#' @param fwhm full width at half maxima
#' @param StatImg SPM of class antsImage 
#' @param mask antsImage mask
#' @param df degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:\item{"T"}{T-field} \item{"F"}{F-field} \item{"X"}{Chi-square field"} 
#' \item{"Z"}{Gaussian field}
#' @param resel-resel values for the mask
#' @return Outputs a statistical value to be used for threshold a SPM
#' @description
#' \code{rft.results} 
#'
#'	
#' @reference 
#' Friston K.J., (1994) Assessing the Significance of Focal Activations Using Their Spatial Extent
#' Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power
#' Worsley K.J., (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' Worlsey K.J., (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' 
#' @author Zachary P. Christensen
#' @kewords rft.pcluster, ants.ec
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
rft.results<-function(D, thresh, ka, fwhm, StatImg, mask, df, fieldType, resel, alpha=.05, report=FALSE, imgdir="./"){
	voxels <-sum(as.array(mask))
	Mfwhm <-mean(fwhm)
	negimg <-as.antsImage(as.array(StatImg)*-1)
	clist <-list()
	if (thresh > max(negimg) && thresh > max(StatImg)){
		stop("No voxels survive threshold")
		}
	if (thresh < max(negimg)){
		clusters <-image2ClusterImages(negimg,0,thresh,Inf)
		for (i in 1:length(clusters)){
			cMask <-getMask(clusters[[i]])
			if(sum(as.array(cMask)) > ka){
				clist <-lappend(clist,cMask)
				}
			}
		}
	if (thresh < max(StatImg)){
		clusters <-image2ClusterImages(StatImg,0,thresh,Inf)
		for (i in 1:length(clusters)){
			cMask <-getMask(clusters[[i]])
			if(sum(as.array(cMask)) > ka){
				clist <-lappend(clist,cMask)
				}
			}
		}
		if(length(clist) < 1){
			cat("No clusters survive threshold
			")
		}else{
			nclusts <-length(clusters)
			clustable<-matrix(nrow=nclusts, ncol=7)
			colnames(clustable)<-c("Voxels", "Cluster-Probability", "Voxel-Probability","Peak-Height","xc", "yc", "zc")
			cnames <-c()
			for (i in 1:nclusts){
				cat("Determing cluster-level statistics:",i,"
				",sep=" ")
				cName <-paste("Cluster",i,sep="")
				cMask <-clist[[i]]
				cvoxs <-sum(as.array(cMask))
				loc <-getCentroids(cMask)
				clust <-maskImage(StatImg,cMask,level=1)
				clist <-lappend(clist,clust)
				cnames <-cbind(cnames, cName)
				names(clist)[length(clist)]<-cName
				pclust <-rft.pcluster(D,cMask,mask,Mfwhm,thresh,df,fieldType)
				cat("Determing voxel-level statistics:",i,"
				",sep=" ")
				maxpeak <-max(clust)
				minpeak <-abs(min(clust))
				if (maxpeak > minpeak){
					peak <-max(clust)
					ec <-rft.ec(peak,fieldType,df)
				}else{
					peak <-min(clust)
					ec <-rft.ec(minpeak,fieldType,df)
					}
				pvox <-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
				}
			names(clist) <-cnames
			rownames(clustable) <-cnames
			report <-list()
			if (report=="TRUE"){
				for (sig in 1:nrow(clustable)){
					if (alpha > clustable[sig,2] | alpha > clustable[sig,2]){
						clustimg <-clist[[sig]]
						xloc <-clustable[sig,5]
						yloc <-clustable[sig,6]
						zloc <-clustable[sig,7]
						pngdir <-paste(imgdir, cnames[sig], ".png", sep="")
						png(pngdir)
						new.plot()
						layout(matrix(c(1,1,1), 1, 3,byrow=TRUE))
						title(main=cnames[sig])
						plot(template, clustimg, slices=c(xloc), axis=1, alpha=1 ,bg="white")
						plot(template, clustimg, slices=c(yloc), axis=2, alpha=1 ,bg="white")
						plot(template, clustimg, slices=c(zloc), axis=3, alpha=1 ,bg="white")
						dev.off()
						}
					}
				}
			}
	z <-list("Clusters", "Stats")
	z$Clusters <-clist
	z$Stats <-clustable
	z
	}
