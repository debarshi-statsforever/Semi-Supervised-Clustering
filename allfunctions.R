lkmeans <- function(data,nclass){
  
  result = list()
  
  n = nrow(data)
  nsub = ceiling(n*0.1)
  
  s = sample(1:n,nsub,replace=F)
  xtrain = data[s,]
  dummy = kmeans(xtrain,nclass,iter.max=1000,algorithm = "Lloyd")
  mu = dummy$centers
  
  #cids <- predict_kmeans(data, mu)
  
  library(clue)
  cids = cl_predict(dummy,data)
  
  
  result$class = as.numeric(cids)
  result$means = mu  
  return(result)
  
  
} 




ldakmeans <- function(image,nclass){
  
  result = list()
  

  
  n = nrow(image)
  nsub = ceiling(n*0.1)
  
  s = sample(1:n,nsub,replace=F)
  xtrain = image[s,]
  dummy = kmeans(xtrain,nclass,iter.max=1000,algorithm = "Lloyd")
  ytrain = as.factor(dummy$cluster)
  train = as.data.frame(cbind(ytrain,xtrain))
  colnames(train) = c("labels","x","y")
  library(MASS)
  fit = lda(labels~.,data=train,prior = rep(1/nclass,nclass))
  xtest = as.data.frame(image)
  colnames(xtest) = c("x","y")
  ypreds = predict(fit,xtest)$class
  mu = fit$means
  
  result$class = as.numeric(ypreds)
  result$centers = mu
  
  return(result)
}


datagen <- function(nclass=4, gsize,rho){
  nclass = nclass ; gsize = gsize  ; n = nclass * gsize
  ids = rep(1:nclass,each = gsize) 
  means1 = c(5,5,-5,-5)  ; means2 = c(5,-5,5,-5) ; rho = rho
  sig = matrix(c(1,rho,rho,1),ncol=2)  ; mult = chol(sig) 
  
  data = list() 
  
  for(i in 1:nclass){ 
    mu1 = means1[i] ; mu2 = means2[i] 
    obs = cbind(rnorm(gsize,mean=mu1),rnorm(gsize,mean=-mu2))
    data[[i]] = obs
  } 
  
  dfraw = do.call(rbind,data) 
  
  for(i in 1:nrow(dfraw)){
    get = dfraw[i,]
    fill = t(mult)%*%get
    dfraw[i,] = fill
  }
  
  
  df = cbind(dfraw,ids)
  colnames(df) = c("x","y","ids")
  
  return(df)
  
}

# start_time <- Sys.time()
# kmfit = kmeans(image, nclass, iter.max=1000, algorithm="Lloyd")
# end_time <- Sys.time()
# 
# end_time - start_time