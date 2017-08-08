## PCA bootstrap

bPCA = function(data,B,percentage){
library(bootstrap)
library(MASS)

#e.g. xdata=mvrnorm(n, c(1,2), matrix(c(1,.99,.99,1),2,2))
n = nrow(data)
d = ncol(data)

s = 1/n*t(data)%*%data
eigenvectors = eigen(s)$vectors
eigenvalues = eigen(s)$values

scores = data%*%eigenvectors       # new data after pca
total.var = sum(eigenvalues)       # total variance
prop.var = eigenvalues/total.var   # proportion of variance
cum.var = cumsum(prop.var)         # culmulative proportion
p=1
#p = sum(cum.var < percentage)+1    # first p components capture 95% variance

#pcdata = princomp(data)$scores     # new data after pca
#v = (princomp(data)$sdev)^2        # proportion of variance
#prop = cumsum(v/sum(v))            # culmulative proportion
#p = sum(prop < per)+1              # first p components capture 95 variance

#pj=princomp(data)$loadings         # projection matrix

pj = solve(eigenvectors)
center = apply(data,2,mean)
ndata = array(,dim=c(n,d,B))
a = array(,dim=c(n,d,B))


for (i in 1:B) {
	# bootstrapping the firt p principle components
	a[,1:p,i] = scores[sample(1:n,replace=T),1:p]
	a[,,i] = cbind(a[,1:p,i],scores[,-(1:p)])
    
    # mapping back
    ndata[,,i]=a[,,i]%*%(pj)+center
}

princ1 = array(,dim=c(d,B))            

for(i in 1:B) {
	scov=1/n*t(ndata[,,i])%*%ndata[,,i]
    princ1[,i]=eigen(scov)$vectors[,1]
    
    if (sum(princ1[,i])>0) {princ1[,i]=10*princ1[,i]}
    else {princ1[,i]=-10*princ1[,i]}
}

princ= apply(princ1,1,mean)
plot(princ,,"l")

#return(ndata)

#theta = function(x, pcdata){mean(pcdata[x,1:p])}
#result = bootstrap(1:n,B,theta,pcdata)
#pc = result$thetastar
#thetahat = cbind(t(pc),rep(0,B))

}
