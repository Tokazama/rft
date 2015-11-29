
rft.matrix <-function(variables,conditions,modelType){
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
		conlist <-list()
		for (n in 1:ncon){
			formula <-~ conditions[,i] - 1
			dm <-model.matrix(formula)
			conlist <-lappend(conlist, dm)
			}
		if (ncon > 1){
			for (n in 2:ncon){
			
			name1 <-names(
	if (modelType=="ancova"|modelType=="anova"){
		congroups <-list()
		nsub <-nrow(conditions)
		#This pads the subjects as zeros and ones
		for (i in 1:ncon){
			cfactor <-as.factor(conditions[,i])
			lnames <-levels(cfactor)
			nlev <-length(lnames)
			condmat <-matrix(nrow=nsub,ncol=nlev)
			ngroup <-0
			for (lev in 1:nlev){
				tmpgroup <-as.numeric(cfactor)
				tmpgroup[tmpgroup !=lev] <-0
				tmpgroup[tmpgroup==lev] <-1
				condmat[,lev] <-tmpgroup
				}
			ngroup <-ngroup + nlev
			conlev <-conlev * nlev
			colnames(condmat) <-lnames
			congroups[i] <-lappend(congroups,condmat)
			}
		DM <-matrix(1L,nrow=nsub)
		colnames(DM) <-"Intercept"
		#This combines the conditions into a design matrix
		congroups1 <-congroups[[1]]
		if (ncon > 1){
			for (i in 2:ncon){
				condmat <-congroups[[i]]
				for (n in 1:ncol(condmat)){
					nlev <-ncol(congroups1)
					tmpmat <-matrix(rep(condmat[,i],nlev),ncol=nlev) * congroups1
					nam1 <-colnames(condmat)[n]
					nam2 <-colnames(congroups1)
					namelist <-list()
					for (nam in 1:length(nam2)){
						tmpname <-print(nam2,":",nam1[nam],sep="")
						namelist <-lappend(namelist,tmpname)
						}
					names(tmpmat) <-namelist
					congroups1 <-tmpmat
					}
				}
			}
		}
		if (modelType=="ancova"){
			tmpmat <-matrix(rep(variables,ncol(congroups1)),ncol=ncol(congroups1))
			dm <-tmpmat * congroups1
			DM <-cbind(dm,DM)
			DF <-(nsub-ngroup)
		}else if(modelType=="anova"){
			DM <-cbind(congroups1,DM)
			DF <-(nsub-ngroup)
		}else{
			DM <-cbind(variables,DM)
			}
		
		if (modelType=="regression"){
			DM <-cbind(variables,1)
			nvar <-ncol(DM)
			DF <-(nsub-nvar)
		}
	fullList <-list()
	fullList <-lappend(fullList,DM)
	fullList <-lappend(fullList,DF)
	fullList <-lappend(fullList,colnames(DM))
	return(fullList)
	}
