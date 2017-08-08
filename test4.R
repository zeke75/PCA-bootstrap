## statistic: sqrt(n)*(||x_bar||^2-||mu||^2), d=2, p=1

n=1000
B=500
method="bootstrap PCA"

par(mfrow=c(1,3))

mu=c(1,2)
sigma=matrix(c(1,.95,.95,1),2,2)
xdata=mvrnorm(n, mu, sigma)
d=ncol(xdata)

s = 1/n*t(xdata)%*%xdata
eigenvectors = eigen(s)$vectors
eigenvalues = eigen(s)$values

scores = xdata%*%eigenvectors      # new data after pca
total.var = sum(eigenvalues)       # total variance
prop.var = eigenvalues/total.var   # proportion of variance
cum.var = cumsum(prop.var)         # culmulative proportion

if (method == "bootstrap PCA"){
	p = sum(cum.var < .95)+1       # first p components capture 95% variance
	} else {p=d};
# bootstrapping the firt p principle components
pj = solve(eigenvectors)
center = apply(xdata,2,mean)
bdata = array(,dim=c(n,d,B))
a = array(,dim=c(n,d,B))
bs = rep(0,B)
w=rep(0,B)

for (i in 1:B) {
	a[,1:p,i] = scores[sample(1:n,replace=T),1:p]
	a[,,i] = cbind(a[,1:p,i],scores[,-(1:p)])
    # mapping back
    bdata[,,i]=a[,,i]%*%(pj)
    bs[i]=(sum((apply(bdata[,,i],2,mean))^2)-sum(center^2))*sqrt(n)
}

hist(bs)
hist(rnorm(500,0,sqrt(35.2)))
var(bs)