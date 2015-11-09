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
rft.thresh<-function(pval,fwhm,mask,df,fieldtype,ka){
	voxels<-sum(as.array(mask))
	D<-length(dim(mask))
	
	if(ka<2){
		Resel<-voxels/(fwhm^3)
		alpha<-pval -1
		stat<-10
		while(alpha < pval){
			stat<-stat-.01
			ec<-ants.ec(stat,fieldtype,df)
			alpha<-Resel*ec[4]
			}
		}
	}
