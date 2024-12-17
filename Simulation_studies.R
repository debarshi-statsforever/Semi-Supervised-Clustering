source("allfunctions.R")
library(parallel) 
library(doParallel)
library(foreach)
library(EMCluster)

nrep = 100 ; ncores = 10 



perform_eval <- function(k){
  
  set.seed(seed=481+k) ; nclass = 4
  mydata = datagen(nclass,100000,-0.5) ; imageorg = mydata[,1:2] ; image = scale(imageorg)
  
  
  kmfit = kmeans(image, nclass, iter.max=1000, algorithm="Lloyd")
  lkmfit = lkmeans(image,nclass)
  ldakmfit = ldakmeans(image,nclass)
  
  
  cids_km = kmfit$cluster ; cids_lkm = lkmfit$class  ;  cids_ldakm = ldakmfit$class
  
 
  
  result_km = as.numeric(RRand(mydata[,3],cids_km)) 
  result_lkm = as.numeric(RRand(mydata[,3],cids_lkm))
  result_ldakm = as.numeric(RRand(mydata[,3],cids_ldakm))
  
  result = c(result_km,result_lkm,result_ldakm)
  
  return(result)
  
  
 
  
}



registerDoParallel(ncores)  



output = foreach(i=1:nrep,.errorhandling = "pass") %dopar% { 
  
  perform_eval(i)
  
} 

stopImplicitCluster()

output = do.call(rbind,output)

summary = colMeans(output)

result = matrix(summary,ncol=3)
colnames(result) = c("Kmeans","LKmeans","LDkmeans") 
rownames(result) = c("RandIndex","Adj-RandIndex","E-index")
result
