dir <- getwd()
dir <- substr(dir,nchar(dir)-6,nchar(dir))

file <- paste('~/fig3','_',dir,'.pdf',sep="")
  pdf(file,width=10,height=5)
###
par(mar=c(3,4,2,4),mfrow=c(1,2),las=1,xaxs='i',yaxs='i')
  require('pomp.orig')
  require('pracma')
  d1 <- as.numeric(as.Date(c('2020-1-31','2021-12-5'))-as.Date("2019-12-31"))/365.25+2020


  load('mle_ncov001.rda')
  nd<-ncol(mle[[1]]@data)
  x <- read.csv('best_ncov.csv',row=1)
mm<-read.csv('model_parameters.csv',r=1)
ad<-sum(mm$type=='par')-length(grep('log.beta',rownames(mm)))
  y <- x
  np<-x$nm+ad+4
print(nd)
  x$bic<- -2*x$loglik + np*log(nd)
  x$aicc<- -2*x$loglik+2*np*nd/(nd-np-1)
  y$bic<-x$bic
  y$aicc<-x$aicc
 j<-c(which.max(x$loglik[1:10]),10+which.max(x$loglik[11:20]))
 print(j)
jj<-j
jj<-c(1,12)
p<-0
Ro2<-as.null()
for(j in c(jj)){
p<-p+1
if(p==4)plot(c(0,1),c(0,1),type='n',ann=F,axe=F)
	print(j)
country<-'MANAUS'

	for(i in 1){
	y[i,] <- x[j[i],]
  load('allregion.rda')
        all[[i]] <- all[[j[i]]]
	}
  dyn.load('ncovmodel.so')
  source('ncovmodel.R')
  mle <- all[[1]]
  m<-mle@pred.mean
  params <- as.numeric(y[1,])
  names(params) <- colnames(y)
  mle@coef <- par.trans(params[match(names(coef(mle)),names(params))])
  n <- 1000
  s <- simulate(mle,nsim=n,seed=100)
  print("succ")

  times <- mle@data[1,]
  death <- mle@data[3,]

  simul2 <- apply(sapply(s,function(x)x@data[3,]),1,median)
  BC <- apply(sapply(s,function(x)x@data[2,]),1,median)
  SH <- apply(sapply(s,function(x)x@data[4,]),1,median)
  s3 <- apply(sapply(s,function(x)x@data[3,]),1,function(x)quantile(x,0.025))
  s4 <- apply(sapply(s,function(x)x@data[3,]),1,function(x)quantile(x,0.975))

    matplot(times,sqrt(cbind(death,simul2,s3,s4)/y$bpop[1]*1e6),col='grey',lty=1,type='l',xlim=d1,ylim=c(0,35),
  ann=F,axe=F,frame=T,xlab='',ylab='')
  mtext(side=3,line=0,cex=0.8,bquote(psi==.(round(y$psi[1],2))))
  polygon(c(times,times[length(times):1]),sqrt(c(s3,s4[length(s3):1])/y$bpop[1]*1e6),col='lightgrey',border=NA)
  points(times,sqrt(death/y$bpop[1]*1e6),type='p',pch=21,bg='white',lty=1,col='red',lwd=1.5,cex=0.75)
  lines(times,sqrt(simul2/y$bpop[1]*1e6),col='black',lwd=2)
 month <- as.numeric(as.Date(c(paste(2020,seq(4,12,by=3),1,sep="-"),paste(2021,seq(1,12,by=3),1,sep="-")))-as.Date('2019-12-31'))/365.25+2020
  axis(1,at=month,lab=month.abb[c(seq(4,12,by=3),seq(1,12,by=3))],cex.axis=0.85,tck=-0.01,padj=-1.2)
  axis(1,at=month[seq(1,7,by=4)],lab=2020:2021,padj=0.2)
  axis(2,hadj=0.6,tck=-0.02)
  mtext(side=2,line=1.8,text=expression(sqrt(weekly~deaths~per~mil)),cex=0.75,las=0,col='black')

  mtext(side=3,line=0.5,adj=-0.1,text=letters[ifelse(p>3,1,0)+p],font=2,cex=1)

        col<-c('darkgreen','brown','red','black','blue')
  if(p==1)
        legend(times[38],28,bty='n',lty=c(1,1,NA,1,1),lwd=1.5,pch=c(2,NA,21,NA,3),pt.bg=c(NA,NA,'white',NA,NA),pt.cex=c(1,1,0.75,1,1),
	       col=col,text.col=col, cex=0.8,legend=c('cum infect','immunized','reported deaths','simulation median','transmission rate'),seg.len=3)
  a <- as.numeric(y[1,match('log.beta',names(y))+c(1:y$nm[1])-1])
    cut <- 2021-(3.36-y$eta[1])/33.6
  time2 <- seq(min(times),cut,by=1/366)
  my <-cubicspline(seq(0,1,len=y$nm[1]/2),a[1:(y$nm[1]/2)],(time2-min(time2))/(max(time2)-min(time2)),endp2nd=TRUE)
  time2 <- seq(cut,max(times),by=1/366)
  my <-c(my,cubicspline(seq(0,1,len=y$nm[1]/2),a[(y$nm[1]/2+1):y$nm[1]],(time2-min(time2))/(max(time2)-min(time2)),endp2nd=TRUE))
  time2<- c(seq(min(times),cut,by=1/366),seq(cut,max(times),by=1/366))

  par(new=T)
  R0<-exp(my)*(1/y$gammah[1])
  plot(time2,R0, xlim=d1,type='l',lwd=0.5,cex=0.8,ylim=c(0,8),col='blue',ann=F,axe=F,frame=F,xlab='',ylab='')
  if(p==1)
  write.csv(cbind(time2,exp(my)),file='beta.csv',quote=F)
  points(time2[seq(1,length(R0),by=15)],R0[seq(1,length(R0),by=15)],pch=3,lwd=0.5,col='blue')
  abline(h=1,lty=2,col='red')
  axis(4,tck=-0.02,hadj=0.2,at=0:5,col.lab='blue',col.ticks='blue')
  mtext(side=4,line=2.5,adj=0.3,text=expression(italic(R)[0](t)),cex=0.75,las=0,col='blue')
  mtext(side=4,line=2,adj=0.9,text='immunized',cex=0.75,las=0,col='brown')
  mtext(side=4,line=2.7,adj=0.9,text='cum infect',cex=0.75,las=0,col='darkgreen')
abline(h=5,col='lightgrey')
axis(4,tck=-0.02,hadj=0.2,at=seq(5,8,len=6),lab=c('',seq(0,1,len=6)[-1]),col.axis='brown') 
points(times,5+3*(0.99-SH/y$bpop[1]),type='l',lwd=3,col='brown')
points(times,5+3*(BC/(y$bpop[1])),type='l',lwd=1,col='darkgreen')
points(times[seq(1,89,by=5)],5+3*(BC[seq(1,89,by=5)]/(y$bpop[1])),pch=2,col='darkgreen')
abline(h=0.7*3+5,col='red',lty=2)
abline(h=0.6*3+5,col='red',lty=2)
abline(h=0.4*3+5,col='red',lty=2)
abline(h=0.2*3+5,col='red',lty=2)
abline(v=(as.Date("2020-10-1")-as.Date("2019-12-31"))/366+2020,lty=2)
abline(v=(as.Date("2020-7-1")-as.Date("2019-12-31"))/366+2020,lty=2)

print(y[1,])
if(p %in% c(1,2)) mtext(side=3,line=-1,adj=0.05,ifelse(p<2,'Manaus','Iquitos'))
if(FALSE){
for(p in c(3:4)){
if(p==3) par(mar=c(3,4,2,4),new=T,fig=c(0.75,1,0.5,1))
if(p==4) par(mar=c(4,4,1,4),new=T,fig=c(0.75,1,0,0.5))
x<-read.csv('best_ncov.csv',r=1)[c(1:10)+10*(p-3),]
x<-x[sort(x$psi,dec=T,index=T)$ix,]
#print(x)
plot(x$psi,x$loglik,type='b',log='',ylab='log likelihood',xlab=expression(psi),ylim=c(ifelse(p==3,-8,-7),1)+round(max(x$loglik)),lwd=2)
abline(h=max(x$loglik)-0.5*qchisq(0.95,1),lty=2)
mtext(side=1,line=3,adj=0.7,'')
par(new=T)
plot(x$psi,x$theta*x$theta*100,type='b',col='red',log='',pch=22,ann=F,axe=F,frame=T,ylim=c(0.25, 1.5))
points(x$psi,x$alpha*x$theta*x$theta*100,type='b',col='blue',pch=23)
axis(4)
mtext(side=4,line=3,las=0,'IFR (%)')
  mtext(side=3,line=0.5,adj=-0.15,text=letters[ifelse(p==3,4,8)],font=2,cex=1)
if(p==3)
legend("bottomleft",bty='n',pch=21:23,col=c('black','red','blue'),legend=c('log likelihood','IFR pre','IFR post'),lty=1,lwd=c(2,1,1))

mtext(side=3,line=-1,adj=0.05,ifelse(p==3,'Manaus','Iquitos'))
}
}
}
dev.off()

