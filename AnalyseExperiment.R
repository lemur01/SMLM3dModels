AnalyseExperiment<-function(dr, dataname, resdr, cropLim, nn){
# crop data as given in cropLim
# separate monomers from clusters based on k=th neares neighbour (given by nn) 
# Model fit: for monomers: CSR and Thomas, for clusters: Thomas and double (nested) Thomas based on envelope tests and pair correlation functions  
  
library(scatterplot3d)
library(spatstat)
library(GET) # library(spptest)

source('Thomas3d.R')

# number of simuations for envelope tests
nrsim<-49
# size of the analyzed patch - for computational efficiency and homogeneity testing
delta<-1000
totR<-NULL

# read 3d data
pb<-read.csv(paste(dr,'/',dataname,'.csv',sep=""),  col.names = c("x", "y", "z"))

# filter z

print(dr)
print(dataname)
print(resdr)


c<-list()
cc<-list()

# Analyze a part of the image - of size delta x delt x range z (without top and bottom part of size 100)

#range x 20000-100000,  y 0-80000
pb<-pb[pb$y>cropLim[3],]
pb<-pb[pb$y<cropLim[4],]
pb<-pb[pb$x<cropLim[2],]
pb<-pb[pb$x>cropLim[1],]


pb<-pb[pb$z>min(pb$z)+500,]
pb<-pb[pb$z<max(pb$z)-500,]


pb$x<-pb$x-min(pb$x)
pb$y<-pb$y-min(pb$y)

fracx<-max(pb$x)/delta
fracy<-max(pb$y)/delta

for (k in 0:(fracx-1)) 
{
for (m in 0:(fracy-1)) 
{
	print(paste(c("k=", k), collapse = " "))
  	print(paste(c("m=", m), collapse = " "))
  

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
	
	# monomers
	monopc <- pcf3est(Pmono)

  	emonoCSR<-envelope(Pmono, fun=K3est, nsim=nrsim,savefuns = TRUE)
  	resmonoCSR <- rank_envelope(emonoCSR)
  	png(paste(resdr,'/monoCSR',k, '_', m, '_', dataname, '.png', sep=""))
  	plot(resmonoCSR)
  	dev.off()
  
  	s<-100
  	prm<-NonIter1GFit(monopc$r, monopc$trans)
  	cmono<-mincontrast(monopc, gThomas3d, c(prm$kappa,prm$s*prm$s),ctrl=list(q = 1/4, p = 2, rmin=NULL, rmax=1000))
 
 	#cmono<-mincontrast(monopc, gThomas3d, c(1/((4*pi*s*s)^(3/2)*(monopc$trans[2]-1)),s*s),ctrl=list(q = 1/4, p = 2, rmin=NULL, rmax=1000))
 
 	npmono<-dim(Pmono$data)
 	mumono <-npmono[1]/volume(Pmono$domain)/cmono$par[1]
 	rho0<-cmono$par[1]
 	if (rho0*volume(w)>1000){
   		png(paste(resdr,'/Fitcmono_',k, '_', m, '_', dataname, '.png', sep=""))
   		plot(cmono)
   		dev.off()
   		presmonoT3 <- -1
   		pintresmonoT3<-c(-1, -1)
 	}
 	else{
  		T3d <- lapply(seq(1,nrsim), simThomas3d, cmono$par[1], mumono,sqrt(abs(cmono$par[2])), c( k*delta,(k+1)*delta, m*delta,(m+1)*delta,min(pb$z[id]),max(pb$z[id])) ) 
  		monoT3d<-envelope(Pmono,fun=K3est, nsim = nrsim,simulate = T3d, savefuns = TRUE)
  		resmonoT3 <- rank_envelope(monoT3d)
  		presmonoT3 <- resmonoT3$p
  		pintresmonoT3<-resmonoT3$p_interval
  
  		df1 <- data.frame(monoT3d)
  		write.table(df1,paste(resdr,'/emonoT3d_',k, '_', m, '_', dataname, '.csv', sep=""))
   
  		png(paste(resdr,'/monoT3_', k, '_', m, '_', dataname, '.png', sep=""))
  		plot(resmonoT3)
  		dev.off()
  		png(paste(resdr,'/Fitcmono_',k, '_', m, '_', dataname, '.png', sep=""))
  		plot(cmono)
  		dev.off()
 	}
 
 
	#nanoclusters
	Pnanocl<-unmark(PClean$feature)
	clpc<- pcf3est(Pnanocl)
	
  	print("nanoclusters")	
	s1<-300
	s2<-100
	prm<-NonIter1GFit(clpc$r, clpc$trans)
	
	# Thomas process model (single cluster)
	c1Thomas<-mincontrast(clpc, gThomas3d, c(prm$kappa,prm$s*prm$s),ctrl=list(q = 1/4, p = 2, rmin=NULL, rmax=1000))
	npcl<-dim(Pnanocl$data)
	mu <-npcl[1]/volume(Pnanocl$domain)/c1Thomas$par[1]
	rho1<-c1Thomas$par[1]
	sigma <-sqrt(c1Thomas$par[2])
	
	if (rho1*volume(w)>100000){
	  png(paste(resdr,'/FitclT_',k, '_', m, '_', dataname, '.png', sep=""))
	  plot(c1Thomas)
	  dev.off()
	  presclT3 <- -1
	  pintresclT3<-c(-1, -1)
	}
	else{
	  T1c_3d <- lapply(seq(1,nrsim), simThomas3d, c1Thomas$par[1], mu, sqrt(abs(c1Thomas$par[2])), c( k*delta,(k+1)*delta, m*delta,(m+1)*delta,min(pb$z[id]),max(pb$z[id])) ) 

	  clT3d<-envelope(Pnanocl,fun=K3est, nsim = nrsim, simulate = T1c_3d, savefuns = TRUE)
	  resclT3 <- rank_envelope(clT3d)
	  presclT3 <- resclT3$p
	  pintresclT3<-resclT3$p_interval
	  
	  png(paste(resdr,'/FitclT_',k, '_', m, '_', dataname, '.png', sep=""))
	  plot(c1Thomas)
	  dev.off()
	  png(paste(resdr,'/clT_',k, '_', m, '_', dataname, '.png', sep=""))
	  plot(resclT3)
	  dev.off()
  	}
  
	# double clusterThomas process model
	c2Thomas<-mincontrast(clpc, gDoubleThomas3d,c(sqrt(prm$kappa),prm$s*prm$s,sqrt(prm$kappa/10),prm$s*prm$s*4),  ctrl=list(q = 1/4, p = 2, rmin=NULL, rmax=1000))
	
	# A^2 = 1/(rho2*(4*pi*(s1+s2))^(d/2))
	# B^2 = 1/(rho2*mu2*(4*pi*s2)^(d/2))
	A<-c2Thomas$par[1]^2
	B<-c2Thomas$par[3]^2
	sigma2 <-min(sqrt(c2Thomas$par[2]),sqrt(c2Thomas$par[4]))
	sigma1 <- sqrt(max(c2Thomas$par[2],c2Thomas$par[4])-sigma2^2)
	if (c2Thomas$par[2]>c2Thomas$par[4]){
	  rho2 <- (1/A)*1/(4*pi*(sigma1^2+sigma2^2))^(3/2)
	  mu2 <-(1/B)*1/(rho2*(4*pi*sigma2^2)^(3/2))
	  
	}
	else{
	  rho2 <- (1/B)*1/(4*pi*(sigma1^2+sigma2^2))^(3/2)
	  mu2 <-(1/A)*1/(rho2*(4*pi*sigma2^2)^(3/2))
	  
	}  
	nu2 <-npcl[1]/volume(Pnanocl$domain)/(mu2*rho2)
	
	#T2c_3d <- lapply(seq(1,nrsim), simDoubleThomas3d, rho2, mu2, nu2, sigma1, sigma2, c( k*delta,(k+1)*delta, m*delta,(m+1)*delta,min(pb$z[id]),max(pb$z[id])) ) 
	#cl2T3d<-envelope(Pnanocl,fun=K3est, nsim = nrsim,simulate = T2c_3d, savefuns = TRUE)
	#rescl2T3 <- rank_envelope(cl2T3d)
	png(paste(resdr,'/FitclDT_',k, '_', m, '_', dataname, '.png', sep=""))

	plot(c2Thomas)
	dev.off()
	
	#res<-cbind(k, m, resmonoCSR$p, resmonoCSR$p_interval,mumono,sqrt(abs(cmono$par[2])), resmonoT3$p, resmonoT3$p_interval,  rho1, mu, sigma, resclT3$p, resclT3$p_interval, rho2, mu, nu, sigma, sigma2, rescl2T3$p, rescl2T3$p_interval)
	res<-cbind(k, m, resmonoCSR$p, resmonoCSR$p_interval[1], resmonoCSR$p_interval[2], cmono$par[1], mumono, sqrt(abs(cmono$par[2])), presmonoT3, pintresmonoT3[1], pintresmonoT3[2],  rho1, mu, sigma, presclT3 , pintresclT3[1], pintresclT3[2] , rho2, mu2, nu2, sigma1, sigma2, 0, 0)
	
	df0 <- data.frame(Pmono)
	write.table(df0,paste(resdr,'/Monomers_',k, '_', m, '_', dataname, '.csv', sep=""), sep = ",")
	
	df0 <- data.frame(Pnanocl)
	write.table(df0,paste(resdr,'/Nanocl_',k, '_', m, '_', dataname, '.csv', sep=""))
	
	
	df1 <- data.frame(monopc)
	write.table(df1,paste(resdr,'/pcmono_',k, '_', m, '_', dataname, '.csv', sep=""))
	
	df2 <- data.frame(clpc)
	write.table(df2,paste(resdr,'/pccl_',k, '_', m, '_', dataname, '.csv', sep=""))
	
	
  cat(res, file= paste(resdr,'/Params_', delta,'_',dataname,'.csv', sep=""), append=TRUE, sep = ",")
  cat("", file= paste(resdr,'/Params_', delta,'_',dataname,'.csv', sep=""), append=TRUE, sep = "\n")

	
  }
}


}
