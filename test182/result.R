i<-1
file<-sprintf("best_ncov%3.3d.csv",i)
while(file.exists(file)==F & i<401)
file<-sprintf("best_ncov%3.3d.csv",i<-i+1)
if(i<401)
{
require('pomp.orig')
source('ncovmodel.R')
y<-read.csv(file)

loglik.max<-y$loglik
j<-1
for(i in c(1:1100))
{
file<-sprintf("best_ncov%3.3d.csv",i)
if(file.exists(file))
{
y1<-read.csv(file)
if(ncol(y1)==ncol(y))
y<-rbind(y,y1)
if(loglik.max<y1$loglik)
{
loglik.max<-y1$loglik
j<-i
}
}
}
colnames(y)[1]<-'town'
#print(y$loglik)
z<-read.csv('best_ncov0.csv',row.names=1)
tw<-row.names(z)
for(i in 1:length(tw))
{
w<-y[y[,1]==tw[i],]
if(dim(w)[1]>0)
z[i,]<-w[which.max(w$loglik),-1]
}
write.csv(z,'current_ncov.csv')
print(z[sort(z$loglik,dec=T,index=T)$ix,])
x<-read.csv('best_ncov0.csv',row.names=1)
y<-read.csv('current_ncov.csv',row.names=1)
print(matrix(round(cbind(y$loglik-x$loglik),4),ncol=2))
print(matrix(round(cbind(ifelse(y$loglik-x$loglik>0,y$loglik-x$loglik,0)),4),ncol=2))
y[x$loglik>y$loglik,]<-x[x$loglik>y$loglik,]
#print(sort(round(cbind(sort(y$loglik-x$loglik)),4)))
write.csv(y,'best_ncov.csv')
}
print(summary((y$BE.0*y$bpop)[1:20]))
