#data1

n=1000
B=500
d=100

mu=rep(0,d)
sigma=diag(d)
sigma[1,1]=d^2
xdata=mvrnorm(n, mu, sigma)


s = 1/n*t(xdata)%*%xdata
eigenvectors = eigen(s)$vectors
eigenvalues = eigen(s)$values

scores = xdata%*%eigenvectors      # new data after pca
total.var = sum(eigenvalues)       # total variance
prop.var = eigenvalues/total.var   # proportion of variance
cum.var = cumsum(prop.var)         # culmulative proportion

