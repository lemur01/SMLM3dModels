AnalyseExperimentShort<-function(dr, dataname, resdr, cropLim, knn){
# performs a simpler analysis than AnalyseExperiment
  
library(scatterplot3d)
library(spatstat)
library(spptest)
library(ggplot2)  
  

source('Thomas3d.R')

nn = knn
nrsim<-49
delta<-10000 # size of analyzed images patches


pb<-read.csv(paste(dr,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))

# filter z

print(dr)
print(dataname)
print(resdr)
pb<-pb[!is.na(pb$x),]
pb<-pb[!is.na(pb$y),]
pb<-pb[!is.na(pb$z),]

c<-list()
cc<-list()

# Analyze a part of the image
pb<-pb[pb$y>cropLim[3],]
pb<-pb[pb$y<cropLim[4],]
pb<-pb[pb$x<cropLim[2],]
pb<-pb[pb$x>cropLim[1],]

pb$x<-pb$x-min(pb$x)
pb$y<-pb$y-min(pb$y)

fracx<-max(pb$x)/delta
fracy<-max(pb$y)/delta

lambda<-matrix(0, fracx, fracy)

# Analyse a patch of the image - for computational efficiency and homegeneity analysis
# if stopped on error, modify here and continue the run
for (k in 0:(fracx))
{
for (m in 0:(fracy))
{
  print("k")
	print(k)
	print("m")
	print(m)

	id <- (pb$x > k*delta)&(pb$x <= (k+1)*delta)&(pb$y > m*delta)&(pb$y <= (m+1)*delta)
	w<-box3(c( k*delta,(k+1)*delta), c( m*delta,(m+1)*delta),c(min(pb$z[id]),max(pb$z[id])))
	P<-pp3(pb$x[id], pb$y[id],pb$z[id],win = w)  # owin =w

	PClean<- split(nnclean(P,nn))
	
	png(paste(resdr,'/Plot2d_', k, '_', m, '_',dataname, '.png', sep=""))
	plot(PClean$noise$data$x,PClean$noise$data$y)
  	points(PClean$feature$data$x,PClean$feature$data$y,col = 19)
  	dev.off()
	
	Pmono<-unmark(PClean$noise)
	png(paste(resdr,'/Clean', k, '_', m, '_',dataname, '.png', sep=""))
	plot(PClean, pch = ".", main = dataname)
	dev.off()
	
	df0 <- data.frame(Pmono)
	write.table(df0,paste(resdr,'/Monomers_',k, '_', m, '_', dataname, '.csv', sep=""), sep = ",")
	
	monopc <- pcf3est(Pmono, rmax =1000, nrval = 10001)
	df <- data.frame(monopc)
	write.table(df,paste(resdr,'/pcmono_',k, '_', m, '_', dataname, '.csv', sep=""))
	
	Pnanocl<-unmark(PClean$feature)
 	clpc<- pcf3est(Pnanocl, rmax =1000, nrval = 10001)
	df1 <- data.frame(clpc)
 	write.table(df1,paste(resdr,'/pccl_',k, '_', m, '_', dataname, '.csv', sep=""))

  	df2 <- data.frame(Pnanocl)
  	write.table(df2,paste(resdr,'/Nanocl_',k, '_', m, '_', dataname, '.csv', sep=""))
  }
}


}
