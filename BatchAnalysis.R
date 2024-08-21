
#source('AnalyseExperiment.R')
source('AnalyseExperimentShort.R')

knn = 10

folder<-'D:/Projects/SMLM/Anoushka/Pair Correlation Data/CA1SR'
dataname <-'Day1CA1SRPCF'

pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
png(paste('D:/Projects/SMLM/Anoushka/Plots','/',dataname, '.png', sep=""))
plot(pb$x, pb$y,pch = ".")
dev.off()


AnalyseExperimentShort(folder, dataname, resfolder, c(30000,70000, 0, 80000), knn)
##############################################################################

folder<-'D:/Projects/SMLM/Anoushka/Pair Correlation Data/CA1SLM'
resfolder<-paste(folder,rfolder, sep="")
dir.create(resfolder)
dataname <-'Day9GLM632CA1SLMPCF'
print(dataname)
pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
png(paste('D:/Projects/SMLM/Anoushka/Plots','/',dataname, '.png', sep=""))
plot(pb$x, pb$y,pch = ".")
dev.off()

AnalyseExperimentShort(folder, dataname, resfolder, c(20000,100000, 0, 60000), knn)

##############################################################################
dataname <-'Day13GLM633CA1SLMPCF'
print(dataname)
pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
folder<-'D:/Projects/SMLM/Anoushka/Pair Correlation Data/CA1SR'
dataname <-'Day1CA1SRPCF'
rfolder<-'/Resk'
resfolder<-paste(folder,rfolder, sep="")
dir.create(resfolder)

pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
png(paste('D:/Projects/SMLM/Anoushka/Plots','/',dataname, '.png', sep=""))
plot(pb$x, pb$y,pch = ".")
dev.off()

AnalyseExperimentShort(folder, dataname, resfolder, c(0,80000, 0, 80000), knn)
##############################################################################



folder<-'D:/Projects/SMLM/Anoushka/Pair Correlation Data/CA3SO'
resfolder<-paste(folder,rfolder, sep="")
dir.create(resfolder)
dataname <-'Day3GLM629CA3SOPCF'
print(dataname)
pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
folder<-'D:/Projects/SMLM/Anoushka/Pair Correlation Data/CA1SR'
dataname <-'Day1CA1SRPCF'
rfolder<-'/Resk'
resfolder<-paste(folder,rfolder, sep="")
dir.create(resfolder)

pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
png(paste('D:/Projects/SMLM/Anoushka/Plots','/',dataname, '.png', sep=""))
plot(pb$x, pb$y,pch = ".")
dev.off()

AnalyseExperimentShort(folder, dataname, resfolder, c(20000,80000, 0, 100000), knn)
##############################################################################

dataname <-'Day7GLM632CA3SOPCF'

pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
plot(pb$x, pb$y,pch = ".")
print(dataname)
AnalyseExperimentShort(folder, dataname, resfolder, c(10000,80000, 0, 80000), knn)
##############################################################################

folder<-'D:/Projects/SMLM/Anoushka/Pair Correlation Data/DGMO'
resfolder<-paste(folder,rfolder, sep="")
dir.create(resfolder)
dataname <-'Day4GLM629DGMOPCF'
print(dataname)
pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
png(paste('D:/Projects/SMLM/Anoushka/Plots','/',dataname, '.png', sep=""))
plot(pb$x, pb$y,pch = ".")
dev.off()

AnalyseExperimentShort(folder, dataname, resfolder, c(40000,80000, 0, 80000), knn)
##############################################################################



folder<-'D:/Projects/SMLM/Anoushka/Pair Correlation Data/DGPO'
resfolder<-paste(folder,rfolder, sep="")
dir.create(resfolder)
dataname <-'Day3GLM629DGPOPCF'
print(dataname)
pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
plot(pb$x, pb$y,pch = ".")

AnalyseExperimentShort(folder, dataname, resfolder, c(35000,85000, 0, 100000), knn)
##############################################################################

dataname <-'Day12GLM629DGPOPCF'
print(dataname)
pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
plot(pb$x, pb$y,pch = ".")

AnalyseExperimentShort(folder, dataname, resfolder, c(20000,80000, 0, 80000), knn)
##############################################################################

folder<-'D:/Projects/SMLM/Anoushka/Pair Correlation Data/CA3SR'
resfolder<-paste(folder,rfolder, sep="")
dir.create(resfolder)
dataname <-'DAY12GLM633CA3SRPCF'
print(dataname)
pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
plot(pb$x, pb$y,pch = ".")

AnalyseExperimentShort(folder, dataname, resfolder, c(10000,80000, 0, 80000), knn)
##############################################################################

folder<-'D:/Projects/SMLM/Anoushka/Pair Correlation Data/CA1SO'
resfolder<-paste(folder,rfolder, sep="")
dir.create(resfolder)
dataname <-'Day2GLM633CA1SOPCF'
print(dataname)
pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
png(paste('D:/Projects/SMLM/Anoushka/Plots','/',dataname, '.png', sep=""))
plot(pb$x, pb$y,pch = ".")
dev.off()

AnalyseExperimentShort(folder, dataname, resfolder, c(40000,80000, 0, 80000), knn)
##############################################################################


dataname <-'Day7GLM633CA1SOPCF'
print(dataname)
pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
png(paste('D:/Projects/SMLM/Anoushka/Plots','/',dataname, '.png', sep=""))
plot(pb$x, pb$y,pch = ".")
dev.off()

AnalyseExperimentShort(folder, dataname, resfolder, c(0,100000, 0, 80000), knn) # many problems
##############################################################################
