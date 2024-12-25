##数据过滤程序

##输入数据：mat_n,mat_x(分别表示总读数数据矩阵和甲基化读数矩阵)
##输出数据：染色体号_4cVS8c_data.Rdata，染色体号_8cVS8c_data.Rdata （例如chr5_4c8c_data.Rdata）
##输出数据形式均为列表形式，列表中存放了过滤位点后并划分好区域后生成的列表testRegion和将mat_n和mat_x按照位点位置排好序且过滤位点后的矩阵treadn和treadx

rm(list=ls())

#过滤函数
filter<-function(mat_n,mat_x,sample8c,sample4c,per){
  #mat_n:total reads matrix
  #mat_x:methylation reads matrix
  #sample8c:Column of 8 cell samples 8细胞样本所在列
  #sample4c:Column of 4 cell samples 4细胞样本所在列
  #per:过滤位点时所设置的缺失值占比值最小值，占比超过per则过滤掉此位点
  
  ##对原数据mat_n和mat_x按照位点位置先进行排序
  mat_x=mat_x[order(as.numeric(rownames(mat_x))),]
  mat_n=mat_n[order(as.numeric(rownames(mat_n))),]
  
  ##将4cell和8cell样本的x和n数据分别表示出来
  cell4_x=mat_x[,sample4c]
  cell4_n=mat_n[,sample4c]
  cell8_x=mat_x[,sample8c]
  cell8_n=mat_n[,sample8c]
  
  ##生成testRegion
  index=c() ###记录已过滤的位点
  precite=as.numeric(rownames(mat_n[1,])) ###记录某区域的第一个位点的位置，这里初始化precite为输入数据中第一个位点的位置
  j=1 ###j记录
  testRegion=vector(mode="list")
  groupindex=c()
  testRegion1=vector(mode="list")
  m=1
  for (i in 1:dim(mat_n)[1]) {
    if((sum(is.na(cell8_n[i,]))/dim(cell8_n)[2])>per | (sum(is.na(cell4_n[i,]))/dim(cell4_n)[2])>per){
      ###NA值大于样本数50%时记录该位点是第几个
      index=c(index,i)
    }else{
      ###如果位点被保留，则进行区域划分，以一个区域最大长度为300bp且一个区域至少包括3个位点为标准划分区域
      if((as.numeric(rownames(mat_n[i,]))-precite)<=300){
        groupindex=c(groupindex,i)
        testRegion[[j]]=vector("list",2)
        names(testRegion[[j]])=c("x","n")
        testRegion[[j]][[1]]=mat_x[groupindex,]
        testRegion[[j]][[2]]=mat_n[groupindex,]
      }else{
        if(j==length(testRegion)){
          if(dim(testRegion[[j]]$n)[1]>=3){
            testRegion1[[m]]=vector("list",2)
            names(testRegion1[[m]])=c("x","n")
            testRegion1[[m]]=testRegion[[j]]
            m=m+1
          }
        }
        precite=as.numeric(rownames(mat_n[i,])) #由于>300bp，故更新区域
        j=j+1
        groupindex=c(i) #区域中选中的位点重新更新为新位点
      }
    }
    if(i%%1000==0)  cat(i," ")
  }
  
  if(length(index)!=0){
    ###在原mat_x和mat_n中删除过滤掉的位点信息，生成新的mat_n和mat_x
    cell8_n=cell8_n[-index,]
    cell8_x=cell8_x[-index,]
    cell4_n=cell4_n[-index,]
    cell4_x=cell4_x[-index,]
    mat_x_sort=mat_x[-index,]
    mat_n_sort=mat_n[-index,]
  }
  
  return(list(treadn=mat_n_sort,treadx=mat_x_sort,testRegion=testRegion1,index=index))
}

#读取数据
setwd("D:/single cell/scFMA/realdata/")
mat_n=read.delim("mat_n")
mat_x=read.delim("mat_x")

#8cell,4cell样本所在列号
sample8c=c(2:12,16,17,19,21:25,27,29,33,35,36,38:40,43:45,47:49,52:55,57:60,62:64,66,69,70,73)
sample4c=c(1,13:15,18,20,26,28,30:32,34,37,41,42,46,50,51,56,61,65,67,68,71,72)
per=0.5 #设置占比值

chr5_data=filter(mat_n,mat_x,sample8c,sample4c,per)
save("chr5_data", file="D:/single cell/scFMA/realdata/chr5_4cVS8c.Rdata")


