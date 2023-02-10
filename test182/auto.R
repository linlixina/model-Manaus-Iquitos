## sharcwrap.R
system('rm -f ncovmodel.o ncovmodel.so')
system('R CMD SHLIB ncovmodel.c')

source('transfer_initial.R')

n<-0
m<-20
j<-c(1:20)
Njob<-as.null()
for(i in 0:n)Njob<-c(Njob,j+m*i)
print(Njob)
for (ijob in Njob) {
  system.command <- sprintf("R CMD BATCH --vanilla '--args %d' xj.R serljob%3.3d.Rout", ijob, ijob);
  print(system.command);
  system(system.command,wait=FALSE);
}
