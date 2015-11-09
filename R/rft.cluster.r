#'
#'
#' calculates the cluster level p-value using RFT
#'
#'
#'
#' 
#' 
#' 
#'
rft.cluster<-function(cMask,bMask,fwhm,stat,df,fieldtype="T",D){
	bvox<-sum(as.array(bmask))
	cvox<-sum(as.array(cmask))
	#Expected number of voxels above a threshold (friston 1994; equation 3)
	if (fieldtype=="T"){
		EN<-bvox*(1-pt(stat,df))
	}else if(fieldtype=="F"){
		EN<-bvox*(1-pf(stat, df[1],df[2]))
	}else if(fieldtype=="X"){
		EN<-bvox*(1-pchisq(stat, df[1],df[2]))
	}else if(fieldtype=="G"){
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
	
	Pnk<-exp(-rfB*cvox^(2/D))
	alpha<-1-exp(-Em*Pnk)
		
	return(alpha)
	}
