#' Cluster level statistics
#'
#' This function calculates the probability of obtaining a cluster the size of 
#' code/{cMask} within a search region of code/{bMask} given the statistical threshold
#' of code/{stat} used to extract said cluster. The degrees of freedom and statistical 
#' field type used to obtain the original statistical map are also required. 
#' 
#' @param u Threshold
#' @param k Cluster size
#' @param c Number of clusters
#' @param mask 
#' @param df degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:
#' \describe{
#' \item{"T"}{T-field} 
#' \item{"F"}{F-field} 
#' \item{"X"}{Chi-square field"} 
#' \item{"Z"}{Gaussian field}
#' }
#' 
#' @return The probability of obtaining the specified cluster
#' @reference 
#' Friston K.J., (1994) Assessing the Significance of Focal Activations Using Their Spatial Extent.
#' Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power.
#' @Author Zachary P. Christensen
#' @note: function currently in beta phase. Waiting for acceptance of peer-reviewed paper
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
rft.pcluster <-function(D, u, k, c, bMask, fwhm, stat, df, fieldType){
	bvox<-sum(as.array(bMask))
	fwhm<-mean(fwhm)
	eps <-.Machine$double.eps
	
	if (fieldType=="T"){
		EN <-bvox*(1-pt(stat,df))
	}else if(fieldType=="F"){
		EN <-bvox*(1-pf(stat, df[1],df[2]))
	}else if(fieldType=="X"){
		EN <-bvox*(1-pchisq(stat, df[1],df[2]))
	}else if(fieldType=="Z"){
		EN <-bvox*(1-qnorm(stat))
	}else{
		stop("A correct fieldtype is required")
		}

	#Based upon the fwhm values or the covariance matrix 
	#of a fields partial derivatives(friston 1994)
	W <-fwhm/((4*log(2))^(1/2))
	
	#The expected number of maxima/clusters (friston 1996; equation 2)
	Em <-bvox*(2*pi)^(-(D+1)/2)*(W^(-D))*(stat^(D-1))*(exp((-stat^2)/2))
	
	#Used to derive the number of expected clusters in a 
	#voxel (friston 1996;equation 3)
	rfB <-(gamma((D/2)+1)*Em/EN)^(2/D)
	
	#Equation 4 (friston 1996) The probability of getting a
	#cluster of the size 'ka' or more above a threshold 'stat'
	#for a single cluster
	Pnk <-exp(-rfB*(ka^(2/D)))
	
	gamma(Em*Pnk)
	# "at high thresholds Ec can be approximated accurately by the average euler characteristic"
	# Applications of Random Fields in Human Brain Mapping-Cao and Worsley 2001

	punc <-1-exp(-Em*Pnk)
	ec <-rft.ec(stat,fieldType,df)
	pstat <-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
	pcor <-1-ppois(c-1,(pstat+eps)*punc)
	z <-list(Pcor=pcor,
		Punc=punc,
		Ec=,
		Ek=
	}
	
	P <-toeplitz( as.numeric(t(EC) * G) ); P[lower.tri(P)] <- 0.0
	P <-mx.exp(P, n)
	P <-P[1,]
