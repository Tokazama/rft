#' @name ants.resel
#' @title Calculates the resels
#'
#'
#' @param mask-statistical value (typically the maxima of a cluster or SPM)
#' @param fwhm- degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @return resel-a vector of the estimated resels
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
#'	negclust<-image2ClusterImages(timg,150,-Inf,-thresh)
#'	negtable<-matrix(ncol=5)
#'	colnames(postable)<-c("Voxels", "Cluster-Probability", "Peak-Height", "Voxel-Probability", "Coordinates")
#'	for (i in 1:length(posclust)){
#'	  cat("Determing negative cluster level statistics:",i,sep=" ")
#'	  cmask<-getMask(negclust[[i]])
#'	  cvoxs<-sum(as.array(cmask))
#'	  pclust<-rft.pcluster(negclust[[i]],mask,fwhm,thresh,df,fieldType)
#'	  peak<-max(posclust[[i]])
#'	  loc<-labelImageCentroids(posclust[[i]])[2]
#'	  resel <-ants.resel(mask,fwhm)
#'	  ec<-ants.ec(stat,fieldType,df)
#'	  pval<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
#'	  negtable[i,]<-c(cvox, pclust, peak, pval, loc$vertices[1],loc$vertices[2],loc$vertices[3])
#'	  clustername<-paste("N-Cluster:",i,sep="")
#'	  rownames(postable[i,])<-c(clustername)
#'	  image<-paste(fileDir,"Nlcuster",i,".nii.gz",sep="")
#'	  antsImageWrite(negclust[[i]],file=image)
#'	  }
#'
#' @export ants.resel
ants.resel <- function(mask, fwhm) {

  P<-sum(as.array(mask))
  dimx <- dim(mask)[1]
  dimy <- dim(mask)[2]
  dimz <- dim(mask)[3]
  rx <- 1/(fwhm[1])
  ry <- 1/(fwhm[2])
  rz <- 1/(fwhm[3])

  Ex <- 0
  Ey <- 0
  Ez <- 0
  Fxy <- 0
  Fxz <- 0
  Fyz <- 0
  cubes <- 0
  
  voxels<-(dimx*dimy*dimz)
  set<-0
  progress <- txtProgressBar(min = 0, max = voxels, style = 3)
  for (i in 1:(dimx)){
      for (j in 1:(dimy)){
          for (k in 1:(dimz)){
              if(mask[i,j,k]==1){
                  Ex <-ifelse(mask[i+1,j,k]==1,Ex+1,Ex)
                  Ey <-ifelse(mask[i,j+1,k]==1,Ey+1,Ey)
                  Ez <-ifelse(mask[i,j,k+1]==1,Ez+1,Ez)
                  Fxy <-ifelse(mask[i+1,j,k]==1 && mask[i,j+1,k]==1 && mask[i+1,j+1,k]==1,Fxy+1,Fxy)
                  Fxz <-ifelse(mask[i+1,j,k]==1 && mask[,j,k+1]==1 && mask[i+1,j,k+1]==1,Fxz+1,Fxz)
                  Fyz <-ifelse(mask[i,j+1,k]==1 && mask[i,j,k+1]==1 && mask[i,j+1,k+1]==1,Fyz+1,Fyz)
                  cubes <-ifelse(mask[i,j,k]==1 && mask[i+1,j,k]==1 && mask[i,j+1,k]==1 && mask[i+1,j+1,k]==1 && mask[i,j,k+1]==1 && mask[i+1,j,k+1]==1 && mask[i,j+1,k+1]==1 && mask[i+1,j+1,k+1]==1,cubes+1,cubes)
                  set<-set+1
                  setTxtProgressBar(progress, set)
                }
            }
        }
    }
  r1 <- (P-(Ex+Ey+Ez)+(Fyz+Fxz+Fxy)-cubes)
  r2 <- (((Ex-Fxy-Fxz+cubes)*rx)+((Ey-Fxy-Fyz+cubes)*ry)+((Ez-Fxz-Fyz+cubes)*rz))
  r3 <- (((Fxy-cubes)*rx*ry)+((Fxz-cubes)*rx*rz)+((Fyz-cubes)*ry*rz))
  r4 <- (cubes*rx*ry*rz)
  resel<-c(r1,r2,r3,r4)
  return(resel)
}
