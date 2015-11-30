
model.rft <-function(variables,conditions,controls,modelType){
	if (modelType=="anova"){
		conditions=variables
		}
	if (modelType=="regression"){
		if (class(variables) !="matrix" | class(variables) !="data.frame"){
			variables <-as.matrix(variables)
			}
		nvar <-ncol(variables)
		formula <-~ variables - 1
		dm <-model.matrix(formula)
		dm <-cbind(dm,1)
		DF <-nrow(dm)-ncol(dm)
	#Creates design matrix based on conditions
	}else if(modelType=="anova" | modelType=="ancova"){
		if (missing(conditions)){
			stop("Both anova and ancova modelTypes must have conditions specified")
			}
		if (class(conditions) !="matrix" | class(conditions) !="data.frame"){
			conditions <-as.matrix(conditions)
			}
		ncon <-ncol(conditions)
		if (ncon > 1){
			formula <-~
			for (n in 1:ncon){
				formula <-formula:conditions[,n]
				}
			dm <-model.matrix(formula-1)
		}else{
			dm <-model.matrix(~conditions-1)
		}
		dm <-cbind(dm,1)
		nsub <-nrow(dim)
		dmcol <-ncol(dm)
		DF <-nsub-dmcol
	}else{
		stop("Must specify appropriate modelType")
	}
	#This creates the interaction between the each specified variable and condition
	if (modelType=="ancova"){
		if (missing(variables)){
			stop("Ancova modelType must have variables specificied")
			}
		nvar <-ncol(variables)
		tmpmat <-matrix(1L,nrow=nsub,ncol=1)
		for (i in 1:nvar){
			var <-variables[,i]
			varmat <-matrix(rep(var,dmcol),ncol=dmcol)
			tmpmat <-cbind((dm *varmat),tmpmat)
			}
		}
	if (missing(controls)){
		colnames(dm)[ncol(dm)] <-"Intercept"
	}else {
		if (class(controls) !="matrix" | class(controls) !="data.frame"){
			controls <-as.matrix(controls)
			}
		dm <-cbind(controls,dm)
		colnames(dm)[ncol(dm)] <-"Intercept"
	}
	design.matrix <-list()
	design.matrix <-lappend(design.matrix,dm)
	design.matrix <-lappend(design.matrix,DF)
	design.matrix <-lappend(design.matrix,colnames(dm))
	names(design.matrix) <-c("DesignMatrix", "DegreesOfFreedom", "D.M.Names")
	return(design.matrix)
	}
