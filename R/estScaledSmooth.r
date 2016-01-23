estScaledSmooth<-function(res,df,mask){
	D <-mask@dimensions
	n <-nrow(res)
	scale <-(1/n)*(df[1]/df[2])
	sr <-res/sqrt(res*res)
	cat("Estimating fwhm/smoothing")
	progress <- txtProgressBar(min = 0, max = subs, style = 3)
	for (i in 1:subs){
	for (i in 1:n){
		img <-makeImage(mask,sr[i,])
		imgar[i,] <-as.array(img)
	}

	close(progress)
	
	dimx <-dim(mask)[1]
	dimy <-dim(mask)[2]
	dimz <-dim(mask)[3]
	d1 <-array(0, dim=c(n,dim(mask)+2))
	d2 <-array(0, dim=c(n,dim(mask)+2))
	m1 <-array(0, dim=c(n,dim(mask)+2))
	m2 <-array(0, dim=c(n,dim(mask)+2))
	maskar <-as.array(mask)
	voxels <-sum(maskar)
	imgar <-as.array(img)
	fwhm <-matrix(0L,nrow=1,ncol=3)
	lambda <-matrix(0L,nrow=3,ncol=3)

	#calculate partial derivatives x
	d1[,2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-imgar
	d2[,1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-imgar
	m1[,2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-maskar
	m2[,1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-maskar
	m3 <-(m1 + m2) == 2
	dx <-((d1 - d2)[m3 == 1])
	
	#calculate partial derivatives y
	d2 <-array(0, dim = dim(img) + 2)
	m2 <-array(0, dim = dim(img) + 2)
	d2[,2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-imgar
	m2[,2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-maskar
	m3 <-(m1 + m2) == 2
	dy <-((d1 - d2)[m3 == 1])

	#calculate partial derivatives z
	d2 <-array(0, dim = dim(img) + 2)
	m2 <-array(0, dim = dim(img) + 2)
	d2[,2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-imgar
	m2[,2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-maskar
	m3 <-(m1 + m2) == 2
	dz <-((d1 - d2)[m3 == 1])
	
	# 
	dx2 <-dx*dx
	dy2 <-dy*dy
	dz2 <-dz*dz
	dxy <-dx*dy
	dxz <-dx*dz
	dyz <-dy*dz
	
	# this should sum all of them into one image create variance covariance
	Vxx <-rowSums(dx2[1:n,])*scale
	Vyy <-rowSums(dy2[1:n,])*scale
	Vzz <-rowSums(dz2[1:n,])*scale
	Vxy <-rowSums(dxy[1:n,])*scale
	Vxz <-rowSums(dxz[1:n,])*scale
	Vyz <-rowSums(dyz[1:n,])*scale
	
	lambda <-cbind(as.matrix(Vxx,ncol=1),as.matrix(Vyy,ncol=1),as.matrix(Vzz,ncol=1))
	fwhm_img <-Vxx*Vyy*Vzz+Vxy*Vyz*Vxz*2-Vxx*Vyz*Vyz-Vxy*Vxy*Vzz-Vxz*Vyy*Vxz
	
	lambda <-sqrt(lambda/(4*log(2)))
	fwhm_img <-sqrt(fwhm_img/(4*log(2))^D)
	
	lambda <-colSums(lambda)/voxels
	fwhm_img <-sum(fwhm_img)/voxels
	
	RESEL <-fwhm_img^(1/D)*(lambda/prod(lambda)^(1/D))
	fwhm <-1/RESEL
	fwhm
	}
