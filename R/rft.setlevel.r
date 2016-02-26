#' @param x  statistical field in matix form (1 x nvox)
#'
#' @references
#' Barnes et al., (2013) Set-level threshold-free tests on the intrinsic volumes of SPMs
#'
#' @export rft.setlevel

rft.setlevel <-function(x, r, sr, df, mask,fwhm,resels){
  D <-mask@dimension
  if (class(x)=="numeric"){
    x <-sparseMatrix(x,nrow=1)
  }
  vox <-ncol(x)
  testEC <-sparseMatrix(nrow=4,ncol=vox)
  res2lkc <-matrix(c(1, (4*log(2))^(1/2), 4*log(2), (4*log(2))^(3/2)), ncol=1)
  
  for (i in x[x !=0]){
    testEC[,i] <-rft.euler(x[,i],df, fieldType)*res2lkc
  }
  
  # might be a quicker way to make testEC as follows
  testUniques <-unique(x)
  for (i in testUniques){
    testEC[,x[x==testUniques[i]]] <-rft.euler(testUniques[i], df, fieldType)
  }
  
  # threshold matrix
  h <-seq(floor(min(x)),ceiling(cmax(x)),by=.2)
  nh <-length(h)
  trialEC <-sparseMatrix(nrow=4,ncol=nh)
  for (i in 1:nh){
    trialEC[,i] <-rft.euler(h[i], df, fieldType="Z")/res2lkc
  }
  
  # get resel estimates based on smoothness
  lkcres <-resels/t(res2lkc)
  
  # no idea what this is yet
  ecR0 <-sparseMatrix(nrow=n+1,ncol=1)
  alleuler2_spm <-sparseMatrix(nrow=n+1,ncol=nh)
  
  n <-nrow(r)
  
  # combine beta matrix and residuals for one loop
  imat <-rbind(x,r)
  for (thr in 1:nh){
    thresh <-h[thr]
    
    # work around negative values (if thresh==0 then everything stays the same)
    threshmat <-imat
    if (thresh > 0){
      threshmat[threshmat < thresh] <-0
    }else if(thresh < 0){
      threshmat[threshmat > thresh] <-0
    }
    
    for (i in 1:n+1){
      img <-makeImage(mask,threshmat[i,])
      resels <-rft.resels(img,fwhm)
    }
    
    
    
  }
  
  # LKC based on average EC through basic regression
   for (i in 1:n+2){
     if (i==n+2){
       Q <-colMeans(imat[2:n+1,]) # average of residuals
     }else{
       Q <-sparseMatrix(imat[i,],nrow=1)
     }
     
     if (i==1){
       Qdash <-Q-LKC0 %*% testEC
       LKC_est <-ginv(testEC[,2:4]) %*% Qdash
     }else{
       Qdash <-Q-LKC0 %*% trialEC
       LKC_est <-ginv(trialEC[,2:4]) %*% Qdash
     }
     
     
     
   }
  
  
  # multivariate test
  ## create design matrix
  gX <-sparseMatrix(nrow=n+1,ncol=2)
  gX[1,1] <-1
  gX[2:n+1,2] <-1
  
  ## beta parameters 
  ### 1st row = test/statistical field LKC estimates
  ### 2nd row = mean of residual field LKC estimates
  
  
  }

