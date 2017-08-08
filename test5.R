## statistic: sqrt(n)*(||x_bar||^2-||mu||^2), d large

n=100
B=500
d=100

par(mfrow=c(3,1))

mu=rep(1,d)
sigma=0.1*diag(d)+matrix(rep(0.9,d^2),d,d)
xdata=mvrnorm(n, mu, sigma)

s = 1/n*t(xdata)%*%xdata
eigenvectors = eigen(s)$vectors
eigenvalues = eigen(s)$values

scores = xdata%*%eigenvectors      # new data after pca
total.var = sum(eigenvalues)       # total variance
prop.var = eigenvalues/total.var   # proportion of variance
cum.var = cumsum(prop.var)         # culmulative proportion

method="bootstrap PCA"

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
hist(rnorm(500,0,sqrt(36000)))
var(bs)