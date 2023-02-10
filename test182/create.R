require('pomp.orig')
n<-20
all<-list(len=n)
loglik<-as.null()
j<-1
for(i in 1:100)
{
if(file.exists(sprintf('mle_ncov%3.3d.rda',i))){
load(sprintf('mle_ncov%3.3d.rda',i))
if(i<21) {
	loglik<-c(loglik,mle[[1]]@loglik)
	all[[i]]<-mle[[1]]
	}
if(i>n & mle[[1]]@loglik>loglik[ifelse(i%%n==0,n,i%%n)])
	all[[ifelse(i%%n==0,n,i%%n)]]<-mle[[1]]
}
}
save(all,file='allregion.rda')
