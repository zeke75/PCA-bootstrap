# statistic: sqrt(n)*(x_bar-mu)

test3 = function(B,xdata,method){

#n=100
#B=50
#method="standard PCA"

#par(mfrow=c(1,2))
#xdata=mvrnorm(n, c(0,0), matrix(c(1,.99,.99,1),2,2))
n=nrow(xdata)
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
	} else {p=d}
	
# bootstrapping the firt p principle components
pj = solve(eigenvectors)
center = apply(xdata,2,mean)
bdata = array(,dim=c(n,d,B))
a = array(,dim=c(n,d,B))
bs = array(,dim=c(d,B))
w=rep(0,B)

for (i in 1:B) {
	a[,1:p,i] = scores[sample(1:n,replace=T),1:p]
	a[,,i] = cbind(a[,1:p,i],scores[,-(1:p)])
    
    # mapping back
    bdata[,,i]=a[,,i]%*%(pj)
    bs[,i]=(apply(bdata[,,i],2,mean)-center)*sqrt(n)
}

for (i in 1:B){
	w[i]=mecdf(t(bs),c(bs[1,i],bs[2,i]))
}
#scatterplot3d(bs[1,],bs[2,],w)
##
#w2=rep(0,B)
#for(i in 1:B){
#	w2[i]=pmvnorm(lower=-Inf, upper=c(bs[1,i],bs[2,i]), mean=c(0,0),
#matrix(c(1,0.99,0.99,1),2,2))
#}
#scatterplot3d(bs[1,],bs[2,],w2)
##
plot(t(bs))
print(cor(t(bs)))
}


par(mfrow=c(1,2))
mu=c(0,0)
sigma=matrix(c(1,.97,.97,1),2,2)
xdata=mvrnorm(100, mu, sigma)
test3(200,xdata,"bootstrap PCA")
test3(200,xdata,"standard bootstrap")



