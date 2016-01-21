rft.aov <-function(x, y, conmat, mask, findSmooth="TRUE", findResels="TRUE", statdir="./"){
	z <-.lm.fit(x,y,tol=tol)
	p <-z$rank
	df <-c(n-p, p-1)
	r <-z$residuals
	fv <-y-r
	b <-z$coefficients
	rss <-sum(r^2)
	mrss <-rss/df[1]
	p1 <-1L:p

	mss <-ss/df[2]

	fstat <-mss/mrss
	 
	 
	if (findSmooth=="TRUE"){
		cat("Estimating FWHM. \n")
		fwhm <-estScaledSmooth(res,mask,df)
		}
	if (findResels=="TRUE"){
		cat("Calculating resels. \n")
		resels <-rft.resels(mask, fwhm)
		}
	ans <-list(design.matrix=x,
		F.imgs=F.imgs,
		conmat=conmat,
		coefficients=coef,
		residuals=r,
		df=df,
		fwhm=fwhm,
		resels=resels)
	ans
}
