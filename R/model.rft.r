
model.rft <-function(variables,conditions,modelType){
	ncon <-ncol(conditions)
	nvar <-ncol(variables)
	conlev <-1
	if (modelType=="regression"){
		if (class(variables) !="matrix" | class(variables) !="data.frame"){
			variables <-as.matrix(variables)
			}
		formula <-~ variables - 1
		dm <-model.matrix(formula)
		dm <-cbind(dm,1)
	}else if(modelType=="anova" | modelType=="ancova"){
		if (class(conditions) !="matrix" | class(conditions) !="data.frame"){
			conditions <-as.matrix(variables)
			}
		if (ncol(conditions) > 1){
			formula <-~
			for (n in 1:ncon){
				formula <-formula*conditions[,n]
				}
			dm <-model.matrix(formula-1)
		}else{
			dm <-model.matrix(~conditions-1)
		}
	}else{
		stop("Must specify appropriate modelType")
	}
	if (modelType=="ancova"){
		dm <-cbind(variables,dm)
		}
	dm <-cbind(dm,1)
	colnames(dm)[ncol(dm)] <-"Intercept"
	design.matrix <-list()
	design.matrix <-lappend(design.matrix,dm)
	design.matrix <-lappend(design.matrix,DF)
	design.matrix <-lappend(design.matrix,colnames(DM))
	names(design.matrix) <-c("DesignMatrix", "DegreesOfFreedom", "D.M.Names")
	return(design.matrix)
	}
