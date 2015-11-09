
rft.results<-function(timg,fwhm,ka,pval){
	cat("Determing clustering threshold using RFT, selected pval, and ka")
	thresh<-rft.thresh(pval,ka,fwhm,mask,df,fieldtype)

	posclust<-image2ClusterImages(timg,150,thresh,Inf)
	postable<-matrix(ncol=5)
	colnames(postable)<-c("Voxels", "Cluster-Probability", "Peak-Height", "Voxel-Probability", "Coordinates")
	for (i in 1:length(posclust)){
		cat("Determing positive cluster level statistics:",i,sep=" ")
		cmask<-getMask(posclust[[i]])
		cvoxs<-sum(as.array(cmask))
		pclust<-rft.pcluster(posclust[[i]],mask,fwhm,thresh,df,fieldtype)
		peak<-max(posclust[[i]])
		loc<-labelImageCentroids(posclust[[i]])[2]
		resel <-ants.resel(mask,fwhm)
		ec<-ants.ec(stat,fieldtype,df)
		pval<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
		postable[i,]<-c(cvox, pclust, peak, pval, loc$vertices[1],loc$vertices[2],loc$vertices[3])
		clustername<-paste("P-Cluster:",i,sep="")
		rownames(postable[i,])<-c(clustername)
		}

	negclust<-image2ClusterImages(timg,150,-Inf,-thresh)
	negtable<-matrix(ncol=5)
	colnames(postable)<-c("Voxels", "Cluster-Probability", "Peak-Height", "Voxel-Probability", "Coordinates")
	for (i in 1:length(posclust)){
		cat("Determing negative cluster level statistics:",i,sep=" ")
		cmask<-getMask(negclust[[i]])
		cvoxs<-sum(as.array(cmask))
		pclust<-rft.pcluster(negclust[[i]],mask,fwhm,thresh,df,fieldtype)
		peak<-max(posclust[[i]])
		loc<-labelImageCentroids(posclust[[i]])[2]
		resel <-ants.resel(mask,fwhm)
		ec<-ants.ec(stat,fieldtype,df)
		pval<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
		negtable[i,]<-c(cvox, pclust, peak, pval, loc$vertices[1],loc$vertices[2],loc$vertices[3])
		clustername<-paste("N-Cluster:",i,sep="")
		rownames(postable[i,])<-c(clustername)
		}
	cluster.table<-rbind(postable,negtable)
	
	return(cluster.table)
	}
