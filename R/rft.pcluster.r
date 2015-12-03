#' Cluster level statistics
#'
#' This function calculates the probability of obtaining a cluster the size of 
#' code/{cMask} within a search region of code/{bMask} given the statistical threshold
#' of code/{stat} used to extract said cluster. The degrees of freedom and statistical 
#' field type used to obtain the original statistical map are also required. 
#' 
#' @param cMask - antsImage mask of the cluster
#' @param bMask - antsImage mask of the brain
#' @param fwhm - full width at half maxima
#' @param stat - statistical value that was used to threshold the image
#' @param df - degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:
#'	"T" - T-field
#'	"F" - F-field
#'	"X" - Chi-square field
#'	"Z" - Gaussian field
#' @return Outputs a 
#' @reference Friston K.J., (1994) Assessing the Significance of Focal Activations Using Their Spatial Extent.
#' @reference Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power.
#'
#' @examples
#' 
#' ## estimatation of a single images smoothness
#' outimg1 <-makeImage(c(10,10,10), rnorm(1000))
#' maskimg <-getMask(outimg1)
#' fwhm <-est.Smooth(outimg1,maskimg)
#' ## create arbitrary threshold and degrees of freedom from hypothetical analysis
#' thresh <-mean(outimg1)
#' df <-4
#' clustimg <-thresholdImage(outimg1, thresh, Inf, inval=1, outval=0)
#' pval <-rft.pcluster(3, clustmask, clustimg, thresh, df, fwhm[[1]], fieldType="T")
#'
#' @export rft.pcluster
rft.pcluster<-function(D,cMask,bMask,fwhm,stat,df,fieldType){
	bvox<-sum(as.array(bMask))
	if(class(cMask)=="numeric"){
		ka<-cMask
	}else{
		ka<-sum(as.array(cMask))
		}
	fwhm<-mean(fwhm)
	
	if (fieldType=="T"){
		EN<-bvox*(1-pt(stat,df))
	}else if(fieldType=="F"){
		EN<-bvox*(1-pf(stat, df[1],df[2]))
	}else if(fieldType=="X"){
		EN<-bvox*(1-pchisq(stat, df[1],df[2]))
	}else if(fieldType=="Z"){
		EN<-bvox*(1-qnorm(stat))
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
