
#source('AnalyseExperiment.R')
source('AnalyseExperimentShort.R')

knn = 10

folder<-'D:/CA1SR'
resfolder<-paste(folder,'/Res/', sep="")
dir.create(resfolder)
dataname <-'Day1CA1SRPCF'

pb<-read.csv(paste(folder,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))
png(paste(resfolder,'/',dataname, '.png', sep=""))
plot(pb$x, pb$y,pch = ".")
dev.off()


AnalyseExperimentShort(folder, dataname, resfolder, c(30000,70000, 0, 80000), knn)
