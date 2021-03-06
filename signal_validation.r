# Read in cmd line args
# Should contain 1) # of iterations
args = commandArgs(TRUE)

startTime = Sys.time()

FidDes00 = read.csv('distances_0_0.csv', header = T)
FidDes22 = read.csv('distances_2_2.csv', header = T)

FidDes00$X = FidDes00$Fiducial = FidDes00$Designs = FidDes00$Plasma.Beta = FidDes00$k = FidDes00$Mach.Number = FidDes00$Solenoidal.Fraction = FidDes00$Virial.Parameter = NULL
FidDes22$X = FidDes22$Fiducial = FidDes22$Designs = FidDes22$Plasma.Beta = FidDes22$k = FidDes22$Mach.Number = FidDes22$Solenoidal.Fraction = FidDes22$Virial.Parameter = NULL

# Ignore any results that include VCS_Break or Tsallis
FidDes00$VCS_Break = FidDes22$VCS_Break = NULL
FidDes00$Tsallis = FidDes22$Tsallis = NULL

# Find the common stats between both files
stats = intersect(colnames(FidDes00), colnames(FidDes22))

y_all = rbind(FidDes00[,stats],FidDes22[,stats])

fct = c(rep(c(1:32), 5), rep(33:64, 5))
fct = as.factor(fct)

nstats = length(stats)

r2 = rep(0,nstats)
for(k in 1:nstats)
{
    y=y_all[,k]
    r2[k]=summary(lm(y~fct))$r.sq
}

sr2=sort.list(r2,dec=T)
sts=names(y_all)

y_all = y_all[,sr2]

nperm = as.numeric(args[1])
s = nstats
pfv = matrix(0,nperm,s)
p2 = rep(0,nstats)

for(k in 1:s)
{
    y=y_all[,k]
    r2[k]=summary(lm(y~fct))$r.sq
    for(j in 1:nperm)
    {
        sp=sample(1:length(y),length(y),replace=T)
        pfv[j,k] = summary(lm(y[sp]~fct[sp]))$r.sq
    }
    p2[k] = sum(pfv[,k]>=.9)
}

p2 = p2 / nperm

out_matrix = matrix(0, nstats, 3)
out_matrix[, 1] = sts[sr2]
out_matrix[, 2] = r2
out_matrix[, 3] = p2

out_table = as.data.frame(out_matrix)
colnames(out_table) = c("Statistic", "R2", "p2")

write.table(out_table, file=paste('metric_validation_', nperm, "_faces_00_22.csv"), sep=",")
