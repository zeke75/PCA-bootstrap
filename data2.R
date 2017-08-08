#data2

n=1000
B=500
d=100

mu=rep(0,d)
sigma=diag(d)
for (i in 1:d){
	for (j in 1:d) {
		if (i>=20 && j>=20 && (i!=j)){sigma[i,j]=1}
	}
}
xdata=mvrnorm(n, mu, sigma)

s = 1/n*t(xdata)%*%xdata
eigenvectors = eigen(s)$vectors
eigenvalues = eigen(s)$values

scores = xdata%*%eigenvectors      # new data after pca
total.var = sum(eigenvalues)       # total variance
prop.var = eigenvalues/total.var   # proportion of variance
cum.var = cumsum(prop.var)         # culmulative proportion