system('rm best_ncov0.csv -f')
system('cp best_ncov.csv best_ncov0.csv')
best<-read.csv('best_ncov0.csv',row.names=1)
townname<-rownames(best)
para<-read.csv("model_parameters.csv",row.names=1)
para$type<-as.factor(para$type)
parameters<-as.data.frame(matrix(nrow=length(para[,'sd']),ncol=nrow(best)+2))
rownames(parameters)=row.names(para)
colnames(parameters)=c('sd','type',townname)
for(town in townname)
{
parameters[,town]<-t(best[town,row.names(para)])
}
parameters[,'sd']<-para[,'sd']
parameters[,'type']<-(levels(para[,'type'])[para[,'type']])
#print(parameters)
write.csv(parameters,'mp_ncov.csv')
