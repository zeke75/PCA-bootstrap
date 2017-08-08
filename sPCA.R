library(MASS)

n=1024
p=2048
par(mfrow=c(2,1))

q = seq(0, 1, length = p)
rho=0.7*dbeta(q,1500,3000)+0.5*dbeta(q,1200,900)+0.5*dbeta(q,600,160)
rho=10*rho/sqrt(sum(rho^2))
plot(rho,,"l")

xdata=array(,dim=c(n,p))

for (i in 1:n){
	xdata[i,]=rnorm(1,0,1)*rho+rnorm(p,0,1)
}

s=1/n*t(xdata)%*%xdata
princ1=10*eigen(s)$vectors[,1]
if (max(princ1)<.5) princ1=-princ1
plot(princ1,,"l")