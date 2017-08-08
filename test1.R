# Bivariate data, standard Bootstrap

library(bootstrap)
library(scatterplot3d)

n=100
B=50
par(mfrow=c(1,2))

xdata=mvrnorm(n, c(0,0), matrix(c(1,.99,.99,1),2,2))
theta=function(x,xdata){(apply(xdata[x,],2,mean)-apply(xdata,2,mean))*sqrt(n)}

result=bootstrap(1:n,B,theta,xdata)
bs=result$thetastar

w=rep(0,B)
for(i in 1:B){
	w[i]=mecdf(t(bs),c(bs[1,i],bs[2,i]))
	}
scatterplot3d(bs[1,],bs[2,],w)

w2=rep(0,B)
for(i in 1:B){
	w2[i]=pmvnorm(lower=-Inf, upper=c(bs[1,i],bs[2,i]), mean=c(0,0), matrix(c(1,0.99,0.99,1),2,2))
}
scatterplot3d(bs[1,],bs[2,],w2)

error=max(abs(w-w2))

