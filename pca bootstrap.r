## PCA bootstrap

bpca = function(data,B,per){
library(bootstrap)
library(MASS)

#e.g. xdata=mvrnorm(n, c(1,2), matrix(c(1,.99,.99,1),2,2))
n = nrow(data)
d = ncol(data)

pcdata = princomp(data)$scores     # new data after pca
v = (princomp(data)$sdev)^2        # proportion of variance
prop = cumsum(v/sum(v))            # culmulative proportion
p = sum(prop < per)+1               # first p components capture 95% variance

pj=princomp(data)$loadings         # projection matrix
center=apply(data,2,mean)

a=array(,dim=c(n,d,B))

ndata=array(,dim=c(n,d,B))

for (i in 1:B) {
	# bootstrapping the firt p principle components
	a[,1:p,i] = pcdata[sample(1:n,replace=T),1:p]
	a[,,i] = cbind(a[,1:p,i],pcdata[,-(1:p)])
    
    # mapping back
    ndata[,,i]=a[,,i]%*%solve(pj)+center
}

return(ndata)

#theta = function(x, pcdata){mean(pcdata[x,1:p])}
#result = bootstrap(1:n,B,theta,pcdata)
#pc = result$thetastar
#thetahat = cbind(t(pc),rep(0,B))

}
