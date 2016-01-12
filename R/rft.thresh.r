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
#'	\item{"voxel"}{using the mask and pval calculates the minimum statistical threshold}
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
rft.thresh <-function(D, StatImg, pval, ka, fwhm, mask, resel, df, fieldType, threshType){
	voxels <-sum(as.array(mask))
	bMask <-mask
	cMask <- ka
	Mfwhm <-mean(fwhm)
	alpha <-pval-1
	stat <-10
	while(alpha < pval){
		stat <-stat-.01
		if (threshType=="cluster"){
			alpha <-rft.pcluster(D, cMask, bMask, Mfwhm, stat, df, fieldType, fdr)
		}else if(threshType=="voxel"){
			ec <-rft.ec(stat, fieldType, df)
			alpha <-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
		}else{
			cat("Must specify appropriate threshType
			")
		}
	}
	x <-StatImg[StatImg > stat]
	if (fdr=="TRUE"){
		if (fieldType=="Z"){
			p <- sort(1 - pnorm(x))
		}else if(fieldType=="T"){
			p <- sort(1 - pt(x, df = df[1]))
		}else if(fieldType=="F"){
			p <- sort(1 - pf(x, df1 = df[1], df2 = df[2]))
		}else if(fieldType=="X"){
			p <-sort(1-pchisq(stat, df[1],df[2]))
		}
		V <- length(p)
		cV <- switch(cV.type, 1, log(V) + 0.5772)
		i <- 1
		while (p[i] <= (i * q)/(V * cV)) i <- i + 1
		i <- max(i - 1, 1)
		if (fieldType=="Z"){
			thresh <- qnorm(1 - p[i])
		}else if(fieldType=="T"){
			thresh <- qt(1 - p[i], df = df[1])
		}else if(fieldType=="F"){
			thresh <- qf(1 - p[i], df1 = df[1], df2 = df[2])
		}else if(fieldType=="X"){
			thresh <-qchisq(1-p[i], df[1],df[2]))
		}
		StatImg[StatImg < thresh] <-0
		}
	z <-list(ThresholdImg=StatImg, RFT.Thresh=stat, FDR.Thresh=thresh)
	return(z)
	}
