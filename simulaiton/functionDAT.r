#--------------------------------------Simplex distribution ----------------------------------------------
rIG<-function(N,mu,lambda)
{
  y<-rchisq(N,1)
  #y<-(rnorm(N,0,1))^2
  a<-(mu^2/(2*lambda))*y
  b<-4*mu*lambda*y+mu^2*y^2
  x1<-mu+a-(mu/(2*lambda))*sqrt(b)
  u<-runif(N)
  x<-rep(0,N)
  for(i in 1:N)
  {
    if(u[i]<mu/(mu+x1[i])){x[i]<-x1[i]}
    else{x[i]<-mu^2/x1[i]}
  }
  return(x)
}

rMIG<-function(N,mu,sigma2,prob)
{
  x2<-rIG(N,mu,1/sigma2)
  z<-rbinom(N,1,prob)
  chi<-rchisq(N,1)
  x3<-sigma2*mu^2*chi
  yy<-x2+z*x3
  return(yy)
}

rsimplex<-function(N,mu,sigma)
{
  yy<-rMIG(N,mu/(1-mu),sigma^2*(1-mu)^2,mu)
  xx<-yy/(1+yy)
  return(xx)
}

psimplex<-function(b,mu,sigma)
{
  inter=function(x)
  {
    (2*pi*sigma^2*x^3*(1-x)^3)^(-0.5)*exp(-(x-mu)^2/(2*sigma^2*x*(1-x)*mu^2*(1-mu)^2))
  }
  
  FF=rep(0,length(b))
  for(i in 1:length(b))
  {
    FF[i]=integrate(inter,0,b[i])$value
  }
  return(FF)
}

dsimplex<-function(x,mu,sigma)
{
  (2*pi*sigma^2*x^3*(1-x)^3)^(-0.5)*exp(-(x-mu)^2/(2*sigma^2*x*(1-x)*mu^2*(1-mu)^2))
}


#-----------------------------------------mixed model----------------------------------------------
rinflatedmodel.beta<-function(n,initialPi0,initialPi1,alfa,beta)
{
  pio1<-rbinom(n,1,initialPi0+initialPi1)
  pio12<-pio1
  pio12[which(pio1==1)]<-rbinom(sum(pio1),1,initialPi1/(initialPi0+initialPi1))
  pio12[which(pio1==0)]<-rbeta(sum(pio1==0),alfa,beta)
  #pio12[which(pio1==0)]<-rsimplex(sum(pio1==0),alfa, beta)
  return(pio12)
}

rinflatedmodel.simplex<-function(n,initialPi0,initialPi1,alfa,beta)
{
  pio1<-rbinom(n,1,initialPi0+initialPi1)
  pio12<-pio1
  pio12[which(pio1==1)]<-rbinom(sum(pio1),1,initialPi1/(initialPi0+initialPi1))
  pio12[which(pio1==0)]<-rsimplex(sum(pio1==0),alfa, beta)
  return(pio12)
}

#----------------------------discrete distribution with MM------------------------------------------------
logMLE<-function(n_read,x_read,n,H,pai)
{
  phi<-phi.est(H)$phi
  tau<-phi.est(H)$tau
  M<-phi.est(H)$M
  
  Pai<-matrix(0,n,H)
  ppai<-t(matrix(rep(pai,each=n),H,n,byrow=TRUE))
  Phi<-t(matrix(rep(phi,each=n),H,n,byrow=TRUE))
  X_read<-t(array(rep(x_read,each=H),c(length(phi),n)))
  N_read<-t(array(rep(n_read,each=H),c(length(phi),n)))
  
  Pai<-ppai*(Phi^X_read)*((1-Phi)^(N_read-X_read))
  
  el<-sum(log(choose(n_read,x_read))+log(rowSums(Pai)))
  return(el)
}



discreteMM.est<-function(n_read,x_read,n,H){
  
  phi<-phi.est(H)$phi
  tau<-phi.est(H)$tau
  M<-phi.est(H)$M
  #inital value
  
  pai=pai.ini=rep(1/H,H) 
  
  
  error=3
  
  while( error>10^(-4) )
  {
    Pai<-matrix(0,n,H)
    nu<-matrix(0,n,H)
    
    ppai<-t(matrix(rep(pai,each=n),H,n,byrow=TRUE))
    Tau<-t(matrix(rep(tau,each=n),M,n,byrow=TRUE))
    Phi<-t(matrix(rep(phi,each=n),H,n,byrow=TRUE))
    X_read<-t(array(rep(x_read,each=H),c(length(phi),n)))
    N_read<-t(array(rep(n_read,each=H),c(length(phi),n)))
    
    Pai<-ppai*(Phi^X_read)*((1-Phi)^(N_read-X_read))
    nu<-Pai/(t(matrix(rep(rowSums(Pai),each=H),H,n)))  
    
    pai.est<-colSums(nu)/n
    if(all(is.na(pai.est))) pai.est<-pai.ini
    
    
    error<-sum(abs(pai.est-pai))
    #ell[k+1]=logMLE(n_read,x_read,n,H,pai.est)
    #error=abs(ell[k+1]-ell[k])/(1+abs(ell[k]))
    #k = k+1
    
    pai<-pai.est
  }
  return(list(pai=pai,phi=phi))
}


boot_discreteMM.est<-function(comb)
{
  x_read<- comb[,1]
  n_read<-comb[,2]
  if(length(n_read)==ntotal8c.new) H=H_select(n8c.new_each,"8cell") else H=H_select(n4c.new_each,"4cell")
  if(length(n_read)==ntotal8c.new) phi=phi.est(H_select(n8c.new_each,"8cell"))$phi else phi=phi.est(H_select(n4c.new_each,"4cell"))$phi
  
  #H<-25
  #phi<-phi8c
  pai.boot<-discreteMM.est(n_read,x_read,n=length(n_read),H)$pai
  mu.boot<-sum(pai.boot*phi)
  
  return(mu.boot)
}

sim_boot_discreteMM.est<-function(comb)
{
  x_read<- comb[,1]
  n_read<-comb[,2]
  if(length(n_read)==ntotal8c.new) H=30 else H=25
  if(length(n_read)==ntotal8c.new) phi=phi.est(30)$phi else phi=phi.est(25)$phi
  
  #H<-25
  #phi<-phi8c
  pai.boot<-discreteMM.est(n_read,x_read,n=length(n_read),H)$pai
  mu.boot<-sum(pai.boot*phi)
  
  return(mu.boot)
}
#----------------------------------------------------------------------------------
rdiscrete.appro<-function(n,phi,pai)
{
  g<-integer(n)
  FP<-cumsum(pai)
  u<-runif(n)
  FP<-c(0,FP,1)
  for(j in 1:(length(phi)))
  {
    ind<-u>FP[j] & u<=FP[j+1]
    g[ind]<-phi[j]
  }
  return(g)
}

phi.est<-function(H)
{
  #length<-1/H
  #tau<-seq(0,1,by=length)
  #if(H==30) a=c(-0.0172,1.0164) else a=c(-0.0208,1.0205)
  if(H==30){
    a=c(-0.017,1.0172)
  }else{
    if(H==12){
      a=c(-0.0452,1.0454)
    }else{
      if(H==20){
        a=c(-0.0262,1.0263)
      }else{
        if(H==10){#H=15
          #a=c(-0.0355,1.0357)
          a=c(-0.0553,1.0555)
        }else{
          if(H==7){
            a=c(-0.0831,1.0833)
          }else{#H=25
            a=c(-0.0206,1.0208)
          }
        }
      }
    }
  }
  tau<-seq(a[1],a[2],by=(a[2]-a[1])/H)
  M<-length(tau)
  phi<-(tau[-M]+tau[2:M])/2 
  return(list(phi=phi,tau=tau,M=M))
}

H_select<-function(cell_length,type)
{
  if(type == "8cell"){
    if(cell_length < 25){
      H = 12
    }else{
      if(cell_length >= 25 & cell_length < 40){H = 20}else{H = 30}
    }
  }else{
    if(cell_length < 15){
      H = 7
    }else{
      if(cell_length >= 15 & cell_length < 20){H = 10}else{H = 20}
    }
  } 
  return(H)
}

