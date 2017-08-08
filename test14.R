#|| ||_1

r=4
pvalue=array(,dim=c(2,r))
D=array(,dim=c(2,r))

par(mfrow=c(2,2))

n=100
B=1000
d=500
	
mu=rep(0,d)
sigma=diag(d)
for (i in 1:d){
	for (j in 1:d) {
			if (i>=d/10 && j>=d/10 && (i!=j)){sigma[i,j]=.995}
	}
}

x=mvrnorm(100000,rep(0,d),sigma)
norm_vec <- function(x) sum(abs(x))
xx=apply(x,1,norm_vec)
#xx=apply(abs(x),1,max)

for (k in 1:r){
	
	plot(ecdf(xx))
	
	# generate simulated data
	xdata=mvrnorm(n, mu, sigma)
	
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
		if (method==1){p=sum(cum.var < .995)+1; print(p)} else(p=d)
		for (i in 1:B) {
			a[,1:p,i] = scores[sample(1:n,replace=T),1:p]
			a[,,i] = cbind(a[,1:p,i],scores[,-(1:p)])
			
			# mapping back
			bdata[,,i]=a[,,i]%*%(pj)
			bs[i]=sum(abs(apply(bdata[,,i],2,mean)-center))*sqrt(n)
		}
		
		
		
		# Two-sample Kolmogorov-Smirnov test
		# bs & xx have same distribution if p-value>0.01
		ks=ks.test(bs,xx)
		
		lines(ecdf(bs),col=45+method)
		
		pvalue[method,k]=ks$p.value
	    D[method,k]=ks$statistic
	}
}


print(sum(pvalue[1,]>.01))
print(sum(pvalue[2,]>.01))
print(sum(pvalue[1,]>pvalue[2,]))

print(pvalue)
print(D)

print(norm_vec(center))
print(max(abs(var(xdata)-sigma)))