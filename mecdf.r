mecdf<-function(data,z) {

n=nrow(data)
nb.temp=0

for(j in 1:n) {
	if (prod(data[j,]<=z)==1){
		nb.temp=nb.temp+1
		}
}
	nb=nb.temp/n
return(nb)
}
