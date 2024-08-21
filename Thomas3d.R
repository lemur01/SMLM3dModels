simThomas3d<-function(id, k, mu, s1, rg)
{
  # simulate 3d Thomas process
  
  gen0<-rpoispp3(k, domain =  box3(c(rg[[1]],rg[[2]]),c(rg[[3]],rg[[4]]),c(rg[[5]],rg[[6]])), nsim = 1)
  
  while(dim(gen0$data)[1]==0){
    gen0<-rpoispp3(k, domain =  box3(c(rg[[1]],rg[[2]]),c(rg[[3]],rg[[4]]),c(rg[[5]],rg[[6]])), nsim = 1)
  }
  nr1<-rpois(length(gen0$data$x),mu)
  x1 <- as.numeric(rep.int(gen0$data$x, nr1 ))+ rnorm( sum(nr1), mean = 0, s1)
  y1 <- as.numeric(rep.int(gen0$data$y, nr1 ))+ rnorm( sum(nr1), mean = 0, s1)
  z1 <- as.numeric(rep.int(gen0$data$z, nr1 ))+ rnorm( sum(nr1), mean = 0, s1)
  
  P<-pp3(x1,y1,z1,win=  box3(c(rg[[1]],rg[[2]]),c(rg[[3]],rg[[4]]),c(rg[[5]],rg[[6]])))
  return (P)
}


gThomas3d<-function(par,r)
{
  # pair correlation function of a Thomas process   
  # gz(z) = 1+(1/k)*1/(4*pi*(s^2))^-d/2*exp(-z*z/(4*(s^2)))
  
  rho1<-par[1] # = (1/k)*1/(4*pi*(s2^2+s1^2))^(-d/2)
  s<-par[2]
  
  d = 3.0
  # The full form is:
  # return(1+(rho2)*exp(-r*r/(4*s2*s2))+(rho1)*exp(-r*r/(4*(s1*s1+s2*s2))))
  
  # simplified form for fitting
  return(1+1/(rho1*(4*pi*s)^(d/2))*exp(-r*r/(4*s)))
  
 
}


gThomas3dMod<-function(par,r)
{ 
  # gz(z) = 1+a*a*((1/k)*1/(4*pi*(s^2))^-d/2*exp(-z*z/(4*(s^2)))-1)
  
  rho1<-par[1] # = (1/k)*1/(4*pi*(s2^2+s1^2))^(-d/2)
  s<-par[2]
  a<-par[3]
  
  d = 3.0
  #return(1+(rho2)*exp(-r*r/(4*s2*s2))+(rho1)*exp(-r*r/(4*(s1*s1+s2*s2))))
  #return(1+a*a*(1/(rho1*(4*pi*s)^(d/2))*exp(-r*r/(4*s))))
  return(1+A*A/(4*pi*s)^(d/2)*exp(-r*r/(4*s)))
  
  #	k<-par[1] # = (1/k)*1/(4*pi*(s2^2+s1^2))^-d/2
  #	s1<-par[2]
  #	mu<-par[3] # =(1/mu*k)*1/(4*pi*s2^2)^-d/2
  #	s2<-par[4]
  
  #	return(1+1/k*1/((4*pi*(s2^2+s1^2))^(-d/2))*exp(-r*r/(4*(s1^2+s2^2)))+1/(mu*k)*1/((4*pi*s2^2)^(-d/2))*exp(-r*r/(4*s2^2)))
}

simDoubleThomas3d<-function(id, k, mu, nu, s1, s2,rg)
{
  # simulate double Thomas process in 3D
  w <-box3(c(rg[[1]],rg[[2]]),c(rg[[3]],rg[[4]]),c(rg[[5]],rg[[6]]))
  gen0<-rpoispp3(k, domain =  w, nsim = 1)
  
  while(dim(gen0$data)[1]==0){
    gen0<-rpoispp3(k, domain =  w, nsim = 1)
  }
 
  nr1<-rpois(length(gen0$data$x),mu)
  #print(nr1)
  x1 <- as.numeric(rep.int(gen0$data$x, nr1 ))+ rnorm( sum(nr1), mean = 0, s1)
  y1 <- as.numeric(rep.int(gen0$data$y, nr1 ))+ rnorm( sum(nr1), mean = 0, s1)
  z1 <- as.numeric(rep.int(gen0$data$z, nr1 ))+ rnorm( sum(nr1), mean = 0, s1)
  
  nr2<-rpois(length(x1),nu)
  print(nr2)
  x2 <- as.numeric(rep.int(x1, nr2 ))+ rnorm( sum(nr2), mean = 0, s2)
  y2 <- as.numeric(rep.int(y1, nr2 ))+ rnorm( sum(nr2), mean = 0, s2)
  z2 <- as.numeric(rep.int(z1, nr2))+ rnorm( sum(nr2), mean = 0, s2)
  
 
  Q<-pp3(x2,y2,z2,win=  box3(c(rg[[1]],rg[[2]]),c(rg[[3]],rg[[4]]),c(rg[[5]],rg[[6]])))
 
  return (Q)
}


gDoubleThomas3d<-function(par,r)
{ # pair correlation function of a Thomas process 
  # gz(z) = 1+1/(mu*k)*1/(4*pi*s2^2)^-d/2*exp(-z*z/(4*s2*s2))+(1/k)*1/(4*pi*(s2^2+s1^2))^-d/2*exp(-z*z/(4*(s2^2+s1^2)))
  

  A<-par[1] # = (1/k)*1/(4*pi*(s2^2+s1^2))^(-d/2)
  s1<-par[2]
  #rho2<-par[3] # =(1/(mu*k))*1/(4*pi*s2^2)^(-d/2)
  B<-par[3] # =(1/(mu*k))*1/(4*pi*s2^2)^(-d/2)
  s2<-par[4]
  
  
  d = 3.0

  return(1+A^2*exp(-r*r/(4*s1))+B^2*exp(-r*r/(4*s2)))
  
  #	k<-par[1] # = (1/k)*1/(4*pi*(s2^2+s1^2))^-d/2
  #	s1<-par[2]
  #	mu<-par[3] # =(1/mu*k)*1/(4*pi*s2^2)^-d/2
  #	s2<-par[4]
  
  #	return(1+1/k*1/((4*pi*(s2^2+s1^2))^(-d/2))*exp(-r*r/(4*(s1^2+s2^2)))+1/(mu*k)*1/((4*pi*s2^2)^(-d/2))*exp(-r*r/(4*s2^2)))
}



gDoubleThomas3dMod<-function(par,r)
{ 
  # pair correlation function of a Thomas process - slightly different parametrisation
  # gz(z) = 1+1/(mu*k)*1/(4*pi*s2^2)^-d/2*exp(-z*z/(4*s2*s2))+(1/k)*1/(4*pi*(s2^2+s1^2))^-d/2*exp(-z*z/(4*(s2^2+s1^2)))
  

  A<-par[1] # = (1/k)*1/(4*pi*(s2^2+s1^2))^(-d/2)
  s1<-par[2]
  #rho2<-par[3] # =(1/(mu*k))*1/(4*pi*s2^2)^(-d/2)
  B<-par[3] # =(1/(mu*k))*1/(4*pi*s2^2)^(-d/2)
  s2<-par[4]
  alpha<-par[5]
  d = 3.0
  return(1+alpha^2*(A^2*exp(-r*r/(4*s1))+B^2*exp(-r*r/(4*s2))))

}