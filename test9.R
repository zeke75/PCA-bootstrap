r=500
pvalue=array(,dim=c(2,r))

for (k in 1:r){
	
	# generate simulated data
	n=100
	B=500
	d=500
	mu=rep(0,d)
	sigma=diag(d)
	for (i in 1:d){
		for (j in 1:d) {
			if (i>=20 && j>=20 && (i!=j)){sigma[i,j]=.995}
		}
	}
	xdata=mvrnorm(n, mu, sigma)
	
	# simulate target distribution 
	x=mvrnorm(5000,rep(0,d),sigma)
	xx=apply(x,1,max)
	
	# PCA
	s = 1/n*t(xdata)%*%xdata
	eigenvectors = eigen(s)$vectors
	eigenvalues = eigen(s)$values
	
	scores = xdata%*%eigenvectors      # new data after pca
	total.var = sum(eigenvalues)       # total variance
	prop.var = eigenvalues/total.var   # proportion of variance
	cum.var = cumsum(prop.var)         # culmulative proportion
	
	
	# bootstrapping the firt p principle components
	pj = solve(eigenvectors)
	center = apply(xdata,2,mean)
	bdata = array(,dim=c(n,d,B))
	a = array(,dim=c(n,d,B))
	bs = rep(0,B)
	
		
	for (method in 1:2){
		if (method==1){p=sum(cum.var < .995)+1} else(p=d)
		for (i in 1:B) {
			a[,1:p,i] = scores[sample(1:n,replace=T),1:p]
			a[,,i] = cbind(a[,1:p,i],scores[,-(1:p)])
			
			# mapping back
			bdata[,,i]=a[,,i]%*%(pj)
			bs[i]=max(apply(bdata[,,i],2,mean)-center)*sqrt(n)
		}
		
		# Two-sample Kolmogorov-Smirnov test
		# bs & xx have same distribution if p-value>0.01
		ks=ks.test(bs,xx)
		pvalue[method,k]=ks$p.value
	}
}


print(sum(pvalue[1,]>.01))
print(sum(pvalue[2,]>.01))

print(sum(pvalue[1,]>pvalue[2,]))


