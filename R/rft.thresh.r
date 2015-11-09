#' A quick heuristic for thresholding a Statistical field using RFT
#'
#'
#' @param pval-
#' @param fwhm-
#' @param mask-
#' @param df-
#' @param fieldtype-
#' @return Outputs the estimated fwhm and covariance matrix that was used to estimate it
#' @examples
#'	
#'
#' @export est.smooth
rft.thresh<-function(pval,ka,fwhm,mask,df,fieldtype){
	voxels <-sum(as.array(mask))
	bMask <-mask
	cMask <- ka
	D<-length(dim(mask))
	fwhm<-mean(fwhm)
	alpha<-pval-1
	stat<-10
	
	while(alpha < pval){
		stat <-stat-.01
		alpha <-rft.cluster(cMask,bMask,fwhm,stat,df,fieldtype="T",D)
		}
	return(stat)
	}
