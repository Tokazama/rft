rftModelFormula <-function(formula, conmat, data){
	mf <-model.frame(mf, data=data)
	y <-model.response(mf, "numeric")
	x <-model.matrix(mt-1, data=mf)
	x <-cbind(1,x)
	colnames(x)[1] <-"Intercept"
	conmat <-cbind(0,conmat)
	colnames(conmat) <-"Intercept"
	if (missing(conmat)){
		
	}else{
		if (colnames(conmat)=="NULL"){
			colnames(conmat) <-colnames(x)
			for (r in 1:nrow(conmat)){
				if (sum(conmat[r,]=0)){
					for (c in 2:ncol(conmat)){
						posname <-paste("", sep="")
						negnames <-paste("", sep="")
						if (conmat[r,c] > 0){
							posname <-paste(posname, colnames(conmat)[c], sep="")
						}else if(conmat[r,c] < 0){
							negname <-paste(negname, colnames(conmat)[c], sep="")
						}
					}
				rownames(conmat)[r] <-paste(posname, ">", negname, sep="")
				}
			}
		}
	}
	if (intercept=="FALSE"){
		x <-x[,2:ncol(x)]
		conmat <-conmat[,2:ncol(conmat)]
		}
	z <-list(y=y, x=x, conmat=conmat)
	}
