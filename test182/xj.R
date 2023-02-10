Npar<-10000
require(pomp.orig)

Args <- commandArgs(trailingOnly=TRUE)
print(Args)
ci<-eval(parse(text = Args))
nreps <- 1
#imagefile <- sprintf("serljob%3.3d.rda", ci) 
id <- 'fit1f1_ncov'
mpfile <- 'mp_ncov.csv'
bestfile <- sprintf("best_ncov%3.3d.csv", ci)  
mlefile <- sprintf("mle_ncov%3.3d.rda", ci) 
filtnames <- c('BS','BE','BI','BT','BR','BR1','BR2','BT1','BT2','BD','CC','BC')
#pfilt<-TRUE
pfilt<-FALSE

x<-read.csv('sari_city2.csv',as.is=T)[c(67:689),]
{
x$time<-as.numeric(as.Date(x[,28])-as.Date('2019-12-31'))/366+2020
x<-x[,c(29,10)]
}

if(ifelse(ci %% 20 ==0, 20, ci%%20)>10)
{
x<-read.csv('peru_excess.csv',r=1)[1093+(67:689),]
x$time<-as.numeric(as.Date(x[,1])-as.Date('2019-12-31'))/366+2020
x<-x[,c(3,2)]
}
city<-'MANAUS'
v<-read.csv('vacc20.csv',r=1)
v<-v[,which(names(v)==city)]
d<-as.numeric(as.Date('2021-1-17')-as.Date('2020-2-27'))
v<-c(rep(0,d),v)
x<-cbind(x[seq(7,nrow(x),by=7),1],apply(matrix(x[,2],ncol=7,byrow=T),1,function(x)sum(x,na.rm=T)))
pop<-2219580
#pop<-2260000
if(ifelse(ci %% 20 ==0, 20, ci%%20)>10)
pop<-426000
last<-function(x)x[length(x)]
v<-v/100
vac<-diff(c(v,last(v)))/(1-v)
import<-NA

x<-as.data.frame(x)
names(x)<-c('time','deaths')
x$confirm<-round(x[,2]*100)
x$dead<-round(x[,2])
x$susc<-round(0.99*pop)


cases<-x[,c('time','confirm','dead','susc')]
row.names(cases)<-1:nrow(cases)
print(head(cases))

parameters <- read.csv(file=mpfile,row.names=1,as.is=TRUE)
districts <- names(parameters)[-c(1,2)]
districts <- districts[ifelse(ci %% 20 ==0, 20, ci%%20)] 

print(districts)
#print(cbind(rownames(parameters),parameters[,districts]))
ifelse(pfilt==TRUE, test.sigma <- parameters$sd/1e20, test.sigma <- parameters$sd)
names(test.sigma) <- rownames(parameters)

source('ncovmodel.R')

parameters[grep('pop',rownames(parameters)),2+ifelse(ci %% 20 ==0, 20, ci%%20)]<- pop
parameters[match('eta',rownames(parameters)),2+ifelse(ci %% 20 ==0, 20, ci%%20)]<- 1e-12
if(ifelse(ci %% 20 ==0, 20, ci%%20)>10)
parameters[match('eta',rownames(parameters)),2+ifelse(ci %% 20 ==0, 20, ci%%20)]<- 3.36

parameters[grep('susc.ini',rownames(parameters)),2+ifelse(ci %% 20 ==0, 20, ci%%20)]<- 0.99
parameters[grep('log.beta',rownames(parameters)),'type']<-'fixed'
nm<-14
parameters[grep('nm',rownames(parameters)),districts]<-nm
parameters[grep('kappah',rownames(parameters)),districts]<-365/12
parameters[grep('log.beta',rownames(parameters))[1:nm],'type']<-'par'
parameters[grep('log.beta',rownames(parameters))[-c(1:nm)],-c(1:2)]<- 3

print(parameters[,c('type',districts)])

ndists <- length(districts)
njobs <- nreps*ndists
#seeds <- lavaseed(n=njobs)
seeds<-ceiling(10^runif(n=njobs,min=7,max=9))
tic<-Sys.time()
joblist <- vector(mode='list',length=njobs)
j <- 0
for (d in districts) {
  mle <- district.pomp(cases,vac,import,rownames(parameters),parameters$type,parameters['nm',d])
  test.params <- parameters[[d]]
  names(test.params) <- rownames(parameters)
  theta.x <- as.matrix(particles(mle,Np=nreps,center=par.trans(test.params),sd=test.sigma)[rownames(parameters),])
  for (r in 1:nreps) {
    j <- j+1
    x <- theta.x[,r]
    names(x) <- rownames(theta.x)
	print(names(x))
    joblist[[j]] <- list(
                         mle=mle,
                         district=d,
                         seed=seeds[j],
                         theta=x,
                         sigma=test.sigma,
                         ivpnames= rownames(parameters)[parameters$type=='ivp'],
                         estnames=rownames(parameters)[parameters$type=='par'],
                         filtnames=filtnames
                         )
  }
}
dyn.load('ncovmodel.so')
done <- 0
for(i in 1:nreps){
joblist[[i]] <-
with(joblist[[i]],
                    {
                      save.seed <- .Random.seed
                      set.seed(seed)
                      x <- mif(
                               mle,
                               Nmif=1,
                               ivps=ivpnames,
                               pars=estnames,
                               stvs=filtnames,
                               alg.pars=list(Np=Npar,CC=4,T0=nrow(cases),cooling.factor=0.95),
                               start=theta,
                               rw.sd=sigma,
                               max.fail=1000,
                               weighted=T,
                               warn=F
                               )
                      result <- list(
                                     mle=x,
                                     district=district,
                                     seed=seed,
                                     rngstate=.Random.seed
                                     )
                      .Random.seed <<- save.seed
                      result
                    }
                    )
}

joblist <- joblist[!sapply(joblist,is.null)]
#save.image(file=imagefile)
done <- done+1

while (done < ifelse(pfilt==TRUE,1,10)) {
for(i in 1:nreps){
joblist[[i]] <-
with(joblist[[i]],
                    {
                        save.seed <- .Random.seed
                        .Random.seed <<- rngstate
                        if (done < 2) {
                          x <- continue(mle,Nmif=11,weighted=T,max.fail=1000)
                        } else {
                          x <- continue(mle,Nmif=11,weighted=T,max.fail=1000)
                        }
                        result <- list(
                                       mle=x,
                                       district=district,
                                       seed=seed,
                                       rngstate=.Random.seed
                                       )
                        .Random.seed <<- save.seed
                        result
                      }
                      )
}
  joblist <- joblist[!sapply(joblist,is.null)]
 # save.image(file=imagefile)
  done <- done+1
}

for(i in 1:nreps){
joblist[[i]] <-
with(joblist[[i]],
                    {
                      save.seed <- .Random.seed
                      .Random.seed <<- rngstate
                      ff <- vector(mode='list')
                      for (n in 1:10)
                        ff[[n]] <- pfilter(mle,max.fail=1000)
                      result <- list(
                                     mle=mle,
                                     district=district,
                                     seed=seed,
                                     rngstate=rngstate,
                                     nfail=sapply(ff,function(x)x$nfail),
                                     loglik=sapply(ff,function(x)x$loglik)
                                     )
                      .Random.seed <<- save.seed
                      result
                    },
                    joblist=joblist
                    )
}
joblist <- joblist[!sapply(joblist,is.null)]
#save.image(file=imagefile)


x <- do.call(
             rbind,
             lapply(
                    joblist,
                    function (x) {
                      cbind(
                            data.frame(
                                       district=x$district,
                                       loglik=mean(x$loglik),
                                       loglik.sd=sd(x$loglik),
                                       nfail.max=max(x$nfail),
                                       nfail.min=min(x$nfail)
                                       ),
                            as.data.frame(t(par.untrans(coef(x$mle))))
                            )
                    }
                    )
             )

best.ind <- tapply(1:nrow(x),x$district,function(k)k[which.max(x$loglik[k])])
best <- x[best.ind,]
rownames(best) <- as.character(best$district)
best$district <- NULL
best <- best[order(rownames(best)),]
write.csv(best,file=bestfile)

mle <- lapply(
              joblist[best.ind],
              function(x)x$mle
              )
names(mle) <- sapply(joblist[best.ind],function(x)x$district)
mle <- mle[order(names(mle))]
save(list='mle',file=mlefile)
toc <- Sys.time()
print(toc-tic)


dyn.unload('ncovmodel.so')

quit()

