require('pomp.orig')

dir <- getwd()
dir <- substr(dir,nchar(dir)-6,nchar(dir))

  d1 <- as.numeric(as.Date(c('2020-1-31','2021-12-5'))-as.Date("2019-12-31"))/365.25+2020

file <- paste('~/fig7','_',dir,'.pdf',sep="")
pdf(file)
at<-as.numeric(as.Date(paste(rep(2020:2021,each=12),1:12,1,sep="-"))-as.Date('2019-12-31'))/366+2020
at1<-as.numeric(as.Date(paste(2020:2021,c(3,1),1,sep="-"))-as.Date('2019-12-31'))/366+2020
n<-1000
y<-read.csv('best_ncov.csv',r=1)
load('allregion.rda')
  source('ncovmodel.R')
  dyn.load('ncovmodel.so')
par(mfrow=c(2,1),mar=c(3,4,1,1),las=1)
p<-0
for(j in c(1:5,11:15))
{

  if(j %in% c(1,11))
  {
 s<-as.null()
p<-p+1
  }
  mle <- all[[j]]
  m<-mle@pred.mean
  params <- as.numeric(y[j,])
  names(params) <- colnames(y)
  mle@coef <- par.trans(params[match(names(coef(mle)),names(params))])
  n <- 200
  s <- append(s,simulate(mle,nsim=n,seed=100))

    if(j %in% c(5,15))
    {
  times <- mle@data[1,]
  print(length(s))
  simul2 <- apply(sapply(s,function(x)x@data[2,]),1,median)

    s3 <- apply(sapply(s,function(x)x@data[2,]),1,function(x)quantile(x,0.025))
    s4 <- apply(sapply(s,function(x)x@data[2,]),1,function(x)quantile(x,0.975))

  matplot(times,(cbind(simul2,s3,s4)/y$bpop[j]),col='grey',lty=1,type='l',xlim=d1,ylim=c(0,1.5),
  ann=F,axe=F,frame=T,xlab='',ylab='')
  polygon(c(times,times[length(times):1]),(c(s3,s4[length(s3):1])/y$bpop[j]),col='lightgrey',border=NA)
  lines(times,(simul2/y$bpop[j]),col='black',lwd=2)
 month <- as.numeric(as.Date(c(paste(2020,seq(4,12,by=3),1,sep="-"),paste(2021,seq(1,12,by=3),1,sep="-")))-as.Date('2019-12-31'))/365.25+2020
  axis(1,at=month,lab=month.abb[c(seq(4,12,by=3),seq(1,12,by=3))],cex.axis=0.85,tck=-0.01,padj=-1.2)
  axis(1,at=month[seq(1,7,by=4)],lab=2020:2021,padj=0.2)
  axis(2,hadj=0.6,tck=-0.02)
  mtext(side=2,line=1.8,text='cum infect per person',cex=0.75,las=0,col='black')
#  abline(h=c(0.4,0.6),lty=2,col='lightgrey')

  mtext(side=3,line=0,adj=-0.1,text=letters[p],font=2,cex=1)
  mtext(side=3,line=-1,adj=0.05,ifelse(j<10,'Manaus','Iquitos'))
  grid(lty=2,col='lightgrey')
    }
}
dev.off()
