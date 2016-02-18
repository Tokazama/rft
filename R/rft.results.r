#' Produces RFT based cluster and voxel level statistics
#' 
#' @param StatImg statistical field image of class antsImage 
#' @param resels resel values for the mask
#' @param fwhm full width at half maxima
#' @param df degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:
#' \describe{
#' \item{"T"}{T-field}
#' \item{"F"}{F-field}
#' \item{"X"}{Chi-square field"}
#' \item{"Z"}{Gaussian field}
#' }
#' @param RPVImg resels per voxel image
#' @param k minimum desired cluster size
#' @param thresh a numeric value to threshold the statistical field or a list containing specification for calculating an appropriate threshold (see examples for details).The following methods may be used:
#' \describe{
#'	\item{"crft"}{computes a threshold per expected cluster level probability }
#'	\item{"prft"}{uses the mask and pval calculates the minimum statistical threshold}
#'	\item{"cfdr"}{uses an uncorrected threshold at the alpha level and then computes and FDR threshold based on cluster maxima}
#'	\item{"pfdr"}{computes the fdr threshold for the entire field of voxels}
#' }
#' @param n number of images in conjunction
#' @param tex output statistical table in latex format (uses package \code{xtable})
#' @param statdir directory where output is saved (by default "./")
#'
#' @return Outputs a statistical value to be used for threshold a statistical field image
#' \describe{
#' \item{"cP-FWE"}{"cluster-level FWE corrected p-value"}
#' \item{"cPfdr"}{"cluster-level FDR corrected p-value"}
#' \item{"cP"}{"cluster-level uncorrected p-value"}
#' \item{"Voxels"}{"number of voxels in this cluster"}
#' \item{"pP-FWE"}{"peak-level FWE corrected p-value"}
#' \item{"pP-FDR"}{"peak-level FDR corrected p-value"}
#' \item{"pP"}{"peak-level uncorrected p-value"}
#' \item{"MaxStat"}{"maximum satistical value in cluster"}
#' \item{"Z"}{"Z statistic that's comparable to MaxStat"}
#' \item{"xc"}{"x-coordinate of cluster centroid"}
#' \item{"yc"}{"y-coordinate of cluster centroid"}
#' \item{"zc"}{"z-coordinate of cluster centroid"}
#' }
#' 
#' @description 
#' 
#' Statistical fields (in the form of and antsImage) are thresholded using either a specified numeric value or a list (i.e. list(method, pval, pp)). 
#' When a list is used a threshold is estimated according to the specified method and p-value (\code{pval}). For example, \code{list(method="prft",pval=.05)}
#' would estimate a threshold where all remaining peaks in the statistical field are at a p-value of .05 or less. If a false discovery rate method is used 
#' (FDR; pfdr or cfdr) then a an additional parameter is specified \code{pp}. This is an initial threshold that allows peak or cluster FDR values to be 
#' estimated. This approach to thresholding allows users to threshold images without using an arbitrary value. 
#' 
#' \code{rft.pval} is used to compute all random field theory (RFT) based statistics and \code{p.adjust} is used to compute all FDR based statistics.
#' 
#' 
#' @References
#' Chumbley J., (2010) Topological FDR for neuroimaging
#' Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power
#' Worsley K.J., (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' 
#' @Author Zachary P. Christensen
#' @kewords rft.pval
#' @note: function currently in beta phase. Waiting for acceptance of peer-reviewed paper
#' @examples
#' 
#' # specifying thresholds
#' tlist <-list("prft", pval=.05) # threshold to create peak values with p-value of .05
#' tlist <-list("crft", pval=.05) # threshold to create clusters with p-value of .05
#' tlist <-list("pfdr", pval=.05, pp=.001) # initial threshold at p-value of .001 followed by peak FDR threshold at p-value of .05
#' tlist <-list("cfdr", pval=.05, pp=.001) # initial threshold at p-value of .001 followed by cluster FDR threshold at p-value of .05
#'
#'
#' @export rft.results
rft.results <-function(StatImg, resels, fwhm, df, fieldType, RPVImg,  k=0, thresh=list(method="prft",pval=.05,pp=.01),n=1, tex=F, statdir=F){
  if (missing(StatImg)){
    stop("Must specify StatImg")
  }
  if (missing(resels)){
    stop("Must specify resels")
  }
  if (missing(fwhm)){
    stop("Must specify fwhm")
  }
  if (missing(df)){
    stop("Must specify degrees of freedom (df)")
  }
  if (missing(fieldType)){
    stop("Must specify fieldType")
  }
  D <-StatImg@dimension
  vox2res <-1/prod(fwhm)
  k <-k*vox2res
  nvox <-length(as.vector(StatImg[StatImg !=0]))
  if (class(thresh)=="list"){
  pval <-thresh[[2]]
  # find minimum threshold value to acheive desired using RFT
    if (thresh[[1]]=="crft" | thresh[[1]]=="prft"){
      u <-max(StatImg)
      if (threshType=="crft"){
        alpha <-rft.pval(D, 1, k, u, n, resels, df, fieldType)$Pcor
      }else if(threshType=="prft"){
        alpha <-rft.pval(D, 1, 0, u, n, resels, df, fieldType)$Pcor
      }
      if (alpha > pval){
        stop("No voxels survive threshold given the parameters")
      }
      while(alpha < pval){
        u <-u-.01
        if (threshType=="crft"){
          alpha <-rft.pval(D, 1, k, u, n, resels, df, fieldType)$Pcor
        }else if(threshType=="prft"){
          alpha <-rft.pval(D, 1, 0, u, n, resels, df, fieldType)$Pcor
        }
      }
    }else if (thresh[[1]]=="cfdr" | thresh[[1]]=="pfdr"){
      pp <-thresh[[3]]
    # FDR based thresholds  
      if (threshType=="cfdr"){
        if (fieldType=="Z"){
          stat <-qnorm(1-pp)
        }else if(fieldType=="T"){
          stat <-qt(1-pp, df[2])
        }else if(fieldType=="F"){
          stat <-qf(1-pp, df[1], df[2])
        }else if(fieldType=="X"){
          stat <-qchisq(1-pp, df[1],df[2])
        }
        statimg <-image2ClusterImages(StatImg, minClusterSize=1,minThresh=stat,maxThresh=Inf)
        cmax <-c()
        for (i in 1:length(clist)){
          cmax <-cbind(cmax,max(clist[[i]]))
        }
      }else if (threshType=="pfdr"){
        cmax <-as.array(StatImg)
      }
      if (fieldType=="Z"){
        p <-sort(1-pnorm(cmax),decreasing=TRUE)
      }else if(fieldType=="T"){
        p <-sort(1-pt(cmax, df = df[1]),decreasing=TRUE)
      }else if(fieldType=="F"){
        p <-sort(1-pf(cmax, df1 = df[1], df2 = df[2]),decreasing=TRUE)
      }else if(fieldType=="X"){
        p <-sort(1-pchisq(cmax, df[1],df[2]),decreasing=TRUE)
      }
      pfdr <-sort((alpha*(1:length(p))/length(p)),decreasing=TRUE)
      i <- 1
      while (p[i] >= pfdr[i]){
        i <-i+1
      }
      if (fieldType=="Z"){
        u <-qnorm(1-p[i])
      }else if(fieldType=="T"){
        u <-qt(1-p[i], df = df[1])
      }else if(fieldType=="F"){
        u <-qf(1-p[i], df1 = df[1], df2 = df[2])
      }else if(fieldType=="X"){
        u <-qchisq(1-p[i], df[1],df[2])
      }
    }
  }else if(class(thresh)=="numeric"){
    u <-thresh
  }else{
    stop("thresh must be a numeric value or a list of specifications for calculating a threshold")
  }
  
  mask <-getMask(StatImg)
  # extract clusters at threshold
  clust <-labelClusters(StatImg,k,u,Inf)
  if (!(missing(statdir))){
    antsImageWrite(clust,file=paste(statdir))
  }
  labs <-unique(clust[clust > 0])
  nclus <-length(labs) # number of clusters
  stats <-matrix(nrow=nclus, ncol=12)
  # stats[1,2] <-nclus
  # cluster coordinates
  stats[1:nclus,10:12] <-getCentroids(clust)[1:nclus,1:3]
  # # footnote
  # Pz <-rft.pval(D,1,0,u,n,c(1,1,1,1),df,fieldType)
  # Pu <-rft.pval(D,1,0,u,n,resels,df,fieldType)
  # psum <-rft.pval(D,1,k,u,n,resels,df,fieldType)
  # summary <-list(u, # Height threshold
  #                Pz,
  #                Pu,
  #                k/vox2res, # Extent threshold
  #                Pn,
  #                P,
  #                Ek/vox2res, # Expected voxels per cluster
  #                Ec*Pn, # Expected number of clusters
  #                df, # degrees of freedom
  #                nvox, # Volume
  #                resels[4],
  #                fwhm #FWHM
  #                )
  # set-level stat
  Pset <-rft.pval(D, nclus, k, u, n, resels, df, fieldType)
  stats <-as.table(stats)
  Eu <-c()
  for (i in 1:nclus){
    cname <-paste("Cluster",i,sep="")
    cat("Calculating statistics for",cname," \n",sep=" ")
    rownames(stats)[i] <-cname
    cmask <-antsImageClone(mask)
    cmask[clust !=labs[i]] <-0
    stats[i,4] <-sum(as.array(cmask)) # number of voxels in cluster
    
    if (missing(RPVImg)){
      # follows isotropic image assumptions
      K <-stats[i,4]*vox2res
    }else{
      # extract resels per voxel in cluster (for non-isotropic image)
      rkc <-RPVImg[cmask==1]
      lkc <-sum(rkc)/stats[i,4]
      iv <-matrix(rft.resels(cmask,c(1,1,1)),nrow=1)
      iv <-iv %*% matrix(c(1/2,2/3,2/3,1),ncol=1)
      K <-iv*lkc
    }
    # max cluster value
    stats[i,8] <-max(StatImg[cmask==1])
    
    # cluster/peak-level stats
    ppeak <-rft.pval(D, 1, 0, stats[i,8], n, resels, df, fieldType) 
    stats[i,5] <-ppeak$Pcor # fwe p-value (peak)
    stats[i,7] <-ppeak$Punc # uncorrected p-value (peak)
    pclust <-rft.pval(D, 1, K, u, n, resels, df, fieldType)
    stats[i,1] <-pclust$Pcor # fwe p-value (cluster)
    stats[i,3] <-pclust$Punc
    Eu <-c(Eu,rft.pval(D, 1, 0, stats[i,8], n, resels, df, fieldType)$Ec)
    
    # comparable Z stat
    Z <-1-qnorm(stats[i,7])
    Pz <-rft.pval(D, 1, 0, Z, n, c(1,1,1,1), df, fieldType)
    Pu <-rft.pval(D, 1, 0, Z, n, resels, df, fieldType)
    
    stats[i,9] <-Z
  }
  # FDR correction (see Chumbley J., (2010) Topological FDR for neuroimaging)
  EuEz <-Eu/rft.pval(D, 1, 0, u, n, resels, df, fieldType)$Ec
  stats[,2] <-p.adjust(stats[,3],"BH") # cluster FDR
  stats[,6] <-p.adjust(EuEz,"BH") # peak FDR
  colnames(stats) <-c("cP-FWE", "cP-FDR", "cP", "Voxels",  
                      "pP-FWE", "pP-FDR", "pP", "MaxStat", 
                      "Z", "xc", "yc", "zc")
  
  results <-list(Statistics=round(stats,4),ClusterImg=clust,Pset$Pcor)
  
  if (tex=="T"){
    texTable <-xtable(stats)
    results <-lappend(results,texTable)
  }
  if (statdir=="T"){
    write.csv(stats, file=paste(statdir))
  }
  return(results)
  print(stats)
}
