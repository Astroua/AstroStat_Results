
# Read in cmd line args
# Should contain 1) # of iterations
args = commandArgs(TRUE)

startTime = Sys.time()

FidDes00 = read.csv('distances_0_0.csv', header = T)
FidFid00 = read.csv('fiducials_0_0.csv', header = T)

# Remove unneeded columns
FidDes00$X = FidDes00$Fiducial = FidDes00$Designs = FidDes00$Plasma.Beta = FidDes00$k = FidDes00$Mach.Number = FidDes00$Solenoidal.Fraction = FidDes00$Virial.Parameter = NULL

# Ignore any results that include VCS_Break
FidDes00$VCS_Break = FidFid00$VCS_Break = FidDes00$Tsallis = FidFid00$Tsallis = NULL

FidFid00$X = FidFid00$Fiducial.1 = FidFid00$Fiducial.2 = NULL
FidFid22$X = FidFid22$Fiducial.1 = FidFid22$Fiducial.2 = NULL

# Ensure we use a common set of statistics across all
stats = intersect(colnames(FidDes00), colnames(FidFid00))

y = rbind(FidFid00[,stats], FidDes00[,stats])
x = c(rep(0, length(FidFid00[,1])),
      rep(1, length(FidDes00[,1])))

nperm = as.numeric(args[1])
nstats = length(stats)
pVals = rep(0, nstats)

for(i in 1:nperm)
{
    ys = y[sample(nrow(y)),]
    for(j in 1:nstats)
    {
        permdata = cbind(ys[j], x)
        data = cbind(y[j], x)
        names(data) = names(permdata) = c("y", "x")
        if(summary(lm(y ~ x, data = permdata))$coefficients[2] > summary(lm(y ~ x, data = data))$coefficients[2])
        {
            pVals[j] = pVals[j] + 1
        }
    }
    if (i %% nperm/10 == 0)
    {
        print(paste(i, '/', nperm, sep=""))
    }
}
names(pVals) = stats
pVals = pVals/nperm
print(pVals)
print(Sys.time() - startTime)
write.table(pVals, file=paste("pValues", nperm, "face_00.csv"), sep=",", col.names = F, row.names=T)
