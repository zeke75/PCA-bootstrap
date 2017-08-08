P=matrix(c(.95,.1,.05,.9),2,2)
H=matrix(c(1/6,1/10,1/6,1/10,1/6,1/10,1/6,1/10,1/6,1/10,1/6,1/2),2,6)
p=c(1,0)
o=c(3,1,5,1,1,6,2,4,6,4,4,6,6,4,4,2,4,5,3,1,1,3,2,1,6,3,1,1,6,4,1,5,2,1,3,3,6,2,5,1,4,4,5,4,3,6,3,1,6,5,6,6,2,6,5,6,6,6,6,6,6,5,1,1,6,6,4,5,3,1,3,2,6,5,1,2,4,5,6,3,6,6,6,4,6,3,1,6,3,6,6,6,3,1,6,2,3,2,6,4)
n=100

delta=array(,dim=c(n,2))
psi=array(,dim=c(n,2))
jstar=array(,dim=n)
guess=array(,dim=n)

for (j in 1:2){
	delta[1,j]=p[j]*H[j,o[1]]
}

for (k in 1:(n-1)){
	for (j in 1:2){
		delta[k+1,j]= max(delta[k,1]*P[1,j],delta[k,2]*P[2,j])*H[j,o[k+1]]
		psi[k+1,j]=which.max(c(delta[k,1]*P[1,j],delta[k,2]*P[2,j]))	
	}
}

jstar[n]=which.max(c(delta[n,1],delta[n,2]))

for (k in (n-1):1){
	jstar[k]=psi[k+1,jstar[k+1]]
}

for (k in n:1){
	if (jstar[k]==1){guess[k]="F"}
	else guess[k]="L"
}

print(data.frame(o,delta[,1],delta[,2],guess))
