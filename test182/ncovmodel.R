logit <- function(p){log(p/(1-p))}
expit <- function(r){1/(1+exp(-r))}
expp<-function(gamma,n1,dt){365*dt*n1/(1.0-exp(-n1*gamma*dt))}
logg<-function(gamma,n1,dt){-log(1.0-(365*dt*n1)/gamma)/(n1*dt)}


par.trans <- function (params) {
  x <- params
  r <- length(dim(x))
  if (r > 1) {
    x <- apply(x,2:r,par.trans)
  } else {
    x['sigmah'] <- log(params['sigmah'])
    x['theta'] <- logit(params['theta'])
    x['p0'] <- logit(params['p0'])
    x['phi'] <- log(params['phi'])
    x['gammah'] <- log(params['gammah'])
    x['psi'] <- logit(params['psi'])
    x['kappah'] <- log(params['kappah'])
    x['eta'] <- log(params['eta'])
    x['alpha'] <- log(params['alpha'])
    x['tau'] <- log(params['tau'])
  }
  x
}

par.untrans <- function (params) {
  x <- params
  r <- length(dim(x))
  if (r > 1) {
    x <- apply(x,2:r,par.untrans)
  } else {
    x['sigmah'] <- exp(params['sigmah'])
    x['theta'] <- expit(params['theta'])
    x['p0'] <- expit(params['p0'])
    x['phi'] <- exp(params['phi'])
    x['gammah'] <- exp(params['gammah'])
    x['psi'] <- expit(params['psi'])
    x['kappah'] <- exp(params['kappah'])
    x['eta'] <- exp(params['eta'])
    x['alpha'] <- exp(params['alpha'])
    x['tau'] <- exp(params['tau'])
  }
  x
}

district.pomp <- function (ncovcases, vac, import,  paramnames, types, nb) {

  pomp(
       nstep=7,
       data=rbind(
         time=ncovcases$time,
         cases.b=ncovcases[,'confirm'],
         death.b=ncovcases[,'dead'],
         susc.b=ncovcases[,'susc']
         ),
       t0=as.double(2*ncovcases$time[1]-ncovcases$time[2]),
       rprocess = function (X, t1, t2, nstep, tbp, bp, tbasis, basis, ...) {
         nvar <- nrow(X)
         np <- ncol(X)
         stateindex <- match(c('BS','BE','BI','BT','BR','BR1','BR2','BT1','BT2','BD','CC','BC'),rownames(X))-1
         parindex <- match(c('sigmah','theta','p0','phi','gammah','psi','kappah','eta','log.beta','bpop','alpha'),rownames(X))-1
         result <- .C("basic_sirs_pois",
                      X = as.double(X),
                      t1 = as.double(t1),
                      t2 = as.double(t2),
                      nstep = as.integer(nstep),
		      t_start=as.double(2*ncovcases$time[1]-ncovcases$time[2]),
		      t_end=as.double(ncovcases$time[length(ncovcases$time)]),
                      nvar = as.integer(nvar),
                      np = as.integer(np),
                      stateindex = as.integer(stateindex),
                      parindex = as.integer(parindex),
		      seas_dim = as.integer(nb),
		      vac = as.double(vac),
		      impt = as.double(import),
                      NAOK = TRUE,
                      DUP = FALSE,
                      PACKAGE = "ncovmodel"
                      )$X
         array(result,dim=dim(X),dimnames=dimnames(X))
       },
       dmeasure = function (X, Y, ...) {
         n <- dim(X)
         np <- n[2]
         ncovindex <- match(c('BC','BD','BS','tau'),rownames(X))-1
         .C("negbin_dmeasure",
            n = as.integer(n),
            index = as.integer(ncovindex),
            X = as.double(X),
            y1 = as.double(Y[2]),
            y2 = as.double(Y[3]),
            y3 = as.double(Y[4]),
            f = double(np),
            DUP = FALSE,
            NAOK = TRUE,
            PACKAGE = "ncovmodel"
            )$f
       },
       rmeasure = function (X, time, ...) {
         n <- dim(X)
         nv <- n[1]
         np <- n[2]
         ncovindex <- match(c('BC','BD','BS','tau'),rownames(X))-1
         matrix(
                .C("negbin_rmeasure",
                   n = as.integer(n),
                   index = as.integer(ncovindex),
                   X = as.double(X),
                   cases = double(3*np),
                   DUP = FALSE,
                   NAOK = TRUE,
                   PACKAGE = "ncovmodel"
                   )$cases,
                3,
                np
                )
       },
       particles = function (Np, center, sd, parnames, fixnames, ivpnames, statenames, bp, ...) {
         pop <- c(center['bpop'])
	nm<-center['nm']
         infe.ini <- center['infe.ini'] 
         susc.ini <- center['susc.ini'] 
         X <- matrix(
                     data=0,
                     nrow=1+length(statenames)+length(parnames)+length(fixnames)+length(ivpnames),
                     ncol=Np,
                     dimnames=list(
                       c(statenames,parnames,fixnames,ivpnames,c('BC')),
#since only one state BC is added, nrow only add 1, other if two states added, nrow add 2
                       NULL
                       )
                     )

	center['theta']<-ifelse(center['theta']>logit(0.006^0.5),logit(0.006^0.5),center['theta'])
	center['theta']<-ifelse(center['theta']<logit(0.002^0.5),logit(0.002^0.5),center['theta'])
	center['alpha']<-ifelse(center['alpha']<log(0.5),log(0.5),center['alpha'])
	center['alpha']<-ifelse(center['alpha']>log(1.8),log(1.8),center['alpha'])
	center[grep('log.beta',names(center))]<-ifelse(center[grep('log.beta',names(center))]>log(587),log(587),
						center[grep('log.beta',names(center))])

#for(i in 1:6)
#{
#        names2<-paste('log.beta',i,sep="")
#        names1<-paste('log.beta','',sep="")
#                if(center[match(names1,names(center))]< center[match(names2,names(center))])
#        center[match(c(names2),names(center))]<- center[match(names1,names(center))]
#}

for(i in 1:(nm/2-1))
{
        names2<-paste('log.beta',center['nm']/2+i,sep="")
        names1<-paste('log.beta',i,sep="")
                if(center[match(names2,names(center))]< center[match(names1,names(center))])
        center[match(c(names1,names2),names(center))]<- 0.5*center[match(names1,names(center))]+
                                                                0.5*center[match(names2,names(center))]
}
	names2<-paste('log.beta',center['nm']/2,sep="")
	names1<-paste('log.beta',center['nm']/2-1,sep="")
	        if(center[match(names2,names(center))]< center[match(names1,names(center))])
        center[match(c(names1,names2),names(center))]<- 0.5*center[match(names1,names(center))]+
                                                                0.5*center[match(names2,names(center))]
#	center[grep('log.beta',names(center))]<-ifelse(center[grep('log.beta',names(center))]<log(60.833),log(60.833),
#						center[grep('log.beta',names(center))])
#	center[grep('log.beta',names(center))[nm]]<-center[grep('log.beta',names(center))[nm-1]]

         X[parnames,] <- rnorm(
                               n=Np*length(parnames),
                               mean=center[parnames],
                               sd=sd[parnames]
                               )
         X[fixnames,] <- center[fixnames]
         
         X[ivpnames,] <- apply(
                               matrix(
                                      center[ivpnames]*rlnorm(
                                                              n=Np*length(ivpnames),
                                                              meanlog=-0.5*sd[ivpnames]^2,
                                                              sdlog=sd[ivpnames]
                                                              ),
                                      nrow=length(ivpnames),
                                      ncol=Np
                                      ),
                               2,
                               function(x)
                                 {
				x[1]<-susc.ini
                                x[2]<-ifelse(x[2]>1/infe.ini,1/infe.ini,x[2])
                                x[3]<-ifelse(x[3]>1/infe.ini,1/infe.ini,x[3])
			        x[c(2:3)]<-mean(x[c(2:3)])
				x[4]<-x[3]*0.1
				x[9]<-x[3]*0.1
				x[10]<-x[3]*0.1
				x[5]<-0.0
				x[11]<-0
				x[6:8]<-(1-sum(x[c(1:4,9,10)]))/3
                                x
				}
                                )
         X[statenames[c(1:11)],] <- round(pop*X[ivpnames[c(1:11)],]) # must round to nearest integer!
         X
       },
       parnames=paramnames[types=='par'],
       ivpnames=paramnames[types=='ivp'],
       fixnames=paramnames[types=='fixed'],
       statenames=sub('.0','',paramnames[types=='ivp'])
       )
}
