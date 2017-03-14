library(lme4)
library(xtable)
library(MASS)
###########

# Fits the 0 to 0 comparison without face as a variable.

resFace0=read.table('distances_0_0.csv',header=T,sep=',')

##Create Indicator for which Fiducial
Cube=rep(c(1:5),each=32)

###############
factors=c('Solenoidal.Fraction','Virial.Parameter','k', 'Mach.Number','Plasma.Beta')
des0=resFace0[,factors]

colnames(des0)=c('sf','vp','k','m','pb')
res0=resFace0[,!(names(resFace0)%in%c('Ind',factors,'Fiducial','Designs', 'X'))]
stats=colnames(res0)

data=cbind(Cube, des0, res0[,stats])

apply(data,2,mean)

##############
numStat=length(stats)
indStat=1:numStat

TVals=matrix(0,nrow=31,ncol=numStat)

for(i in 1:numStat)
{
    print(stats[i])
    mod=paste(stats[i],'~sf*vp*k*m*pb+(1|Cube)',sep='')
    ReMod=lmer(mod,REML=F,data=data)

    tvals = matrix(summary(ReMod)$coeff[-1,3,drop=F])

    if (length(tvals) != length(TVals[,i]))
    {
        print("Fit to statistic failed! Output is all zeros!")
        TVals[,i] = rep(0, length(TVals[,i]))
    } else {
        TVals[,i]=matrix(summary(ReMod)$coeff[-1,3,drop=F])
    }
}

TVals=data.frame(TVals)
colnames(TVals)=stats
rownames(TVals)=rownames(summary(ReMod)$coeff[-1,3,drop=F])

write.table(TVals,'ResultsFactorial.csv',sep=',')
write.table(data, 'DataforFits.csv', sep=',')
