
rm(list=ls())
source("D:/single cell/functionDAT.r")
setwd("D:/single cell/scFMA/")
## 有差异实验
treadn<- read.delim("simulation/cpg_chr21_all_n.txt")
dim(treadn) # 645756 76

## value setting

n1<-48
n2<-25
H8c<-30
H4c<-25

#-------------------------------------------generate simululation data---------------------------------------------
res.8c<-matrix(0,1000,(H8c+1))
res.4c<-matrix(0,1000,(H4c+1))
pio=vector(mode="list")

testRegion=vector(mode="list")
for(sim in 1:1000)
{ 
  #generate simulation data
  qujian<-10
  
  eachbinn48c<-treadn[sample(dim(treadn)[1],qujian,replace =FALSE),4:76]
  
  for (j in 1:73) {
    nnumb<-c(treadn[!is.na(treadn[,(j+3)]),(j+3)])
    eachbinn48c[,j]<-sample(nnumb,qujian,replace =TRUE)
  }
  
  
  eachbinx48c<-eachbinn48c
  pio.8c<-NULL
  pio.4c<-NULL
  for (i in 1:qujian) {
    ###  8 cell
    initialPi08c<-0.3
    initialPi18c<-0
    initialalfa8c<-0.6
    initialbeta8c<-4
    binn8c<-eachbinn48c[i,1:n1]
    binn8c1<-binn8c[!is.na(binn8c)]
    pio1.8c<-rbinom(length(binn8c1),1,initialPi08c+initialPi18c)
    pio12.8c<-pio1.8c
    pio12.8c[which(pio1.8c==1)]<-rbinom(sum(pio1.8c),1,initialPi18c/(initialPi08c+initialPi18c))
    pio12.8c[which(pio1.8c==0)]<-rsimplex(sum(pio1.8c==0),initialalfa8c, initialbeta8c)
    
    binx8c<-binn8c1
    binx8c[pio12.8c==0]<-0 
    nn<-binx8c[pio12.8c>0 & pio12.8c<1]
    binx8c[pio12.8c>0 & pio12.8c<1]<-mapply(rbinom, n = rep(1, length(nn)), size =nn, prob =pio12.8c[pio12.8c>0 & pio12.8c<1])#rbinom(length(nn),nn,pio12[pio12.8c>0 & pio12.8c<1])
    binn8c[!is.na(binn8c)]<-binx8c
    eachbinx48c[i,1:n1]<-binn8c
    
    ###  4 cell
    initialPi04c<-0.3
    initialPi14c<-0
    initialalfa4c<-0.6
    initialbeta4c<-4
    binn4c<-eachbinn48c[i,(n1+1):(n1+n2)]
    binn4c1<-binn4c[!is.na(binn4c)]
    pio1.4c<-rbinom(length(binn4c1),1,initialPi04c+initialPi14c)
    pio12.4c<-pio1.4c
    pio12.4c[which(pio1.4c==1)]<-rbinom(sum(pio1.4c),1,initialPi14c/(initialPi04c+initialPi14c))
    pio12.4c[which(pio1.4c==0)]<-rsimplex(sum(pio1.4c==0),initialalfa4c, initialbeta4c)
    
    
    binx4c<-binn4c1
    binx4c[pio12.4c==0]<-0
    nn<-binx4c[pio12.4c>0 & pio12.4c<1]
    binx4c[pio12.4c>0 & pio12.4c<1]<-mapply(rbinom, n = rep(1, length(nn)), size =nn, prob =pio12.4c[pio12.4c>0 & pio12.4c<1])#rbinom(length(nn),nn,pio12.4c[pio12.4c>0 & pio12.4c<1])
    binn4c[!is.na(binn4c)]<-binx4c
    eachbinx48c[i,(n1+1):(n1+n2)]<-binn4c
    
    pio.8c<-c(pio.8c,pio12.8c)
    pio.4c<-c(pio.4c,pio12.4c)
    
    
  }
  
  testRegion[[sim]]=vector("list",2)
  names(testRegion[[sim]])=c("x","n")
  testRegion[[sim]][[1]]=eachbinx48c
  testRegion[[sim]][[2]]=eachbinn48c
  
  eachbinx48c<-as.matrix(eachbinx48c)
  eachbinn48c<-as.matrix(eachbinn48c)
  eachbinn8c<-eachbinn48c[,1:n1]
  eachbinx8c<-eachbinx48c[,1:n1]
  eachbinn4c<-eachbinn48c[,(n1+1):(n1+n2)]
  eachbinx4c<-eachbinx48c[,(n1+1):(n1+n2)]
  
  
  pio[[sim]]=vector("list",2)
  names(pio[[sim]])=c("pio.8c","pio.4c")
  pio[[sim]][[1]]<-pio.8c
  pio[[sim]][[2]]<-pio.4c
  
  #--------------------------parameter estimate------------------------------------------------
  start<-proc.time()[1]
  resulteach8c<-matrix(0,qujian,H8c)
  resulteach4c<-matrix(0,qujian,H4c)
  ###8 cell 
  
  for(i in 1:qujian)
  {
    
    eight_n<-(eachbinn48c[i,1:n1])
    eight_x<-(eachbinx48c[i,1:n1])
    
    
    discreteMM.est8c<-discreteMM.est(eight_n,eight_x,n1,H8c)
    pai8c<-discreteMM.est8c$pai  #function(n_read,x_read,n,H)
    phi8c<-discreteMM.est8c$phi
    
    resulteach8c[i,]<-pai8c
    
    
    ###4 cell 
    
    four_n<-(eachbinn48c[i,(n1+1):(n1+n2)])
    four_x<-(eachbinx48c[i,(n1+1):(n1+n2)])
    
    
    discreteMM.est4c <- discreteMM.est(four_n,four_x,n2,H4c) 
    pai4c<- discreteMM.est4c$pai
    phi4c<- discreteMM.est4c$phi
    
    resulteach4c[i,]<-pai4c 
  }
  
  end<-proc.time()[1]
  time<-end-start
  result8c<- c(time,apply(resulteach8c,2,mean))  #c(time1,pai8c)
  result4c<-c(time,apply(resulteach4c,2,mean))#c(time2,pai4c) 
  
  res.8c[sim,]<-result8c
  res.4c[sim,]<-result4c
  
}

res<-list(res.8c=res.8c,res.4c=res.4c,testRegion=testRegion,pio=pio,phi8c=phi8c,phi4c=phi4c)
save(res, file="simulation/nodifference.simplex_est.Rdata")


