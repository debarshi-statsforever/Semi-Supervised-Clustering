sourceCpp("calc_parallel.cpp") 


get_quantized_image <- function(image,cids,nclass) {
  
  imgn <- image
  dim(imgn) <- c(length(cids),3)
  sm = matrix(nrow=nclass,ncol=3)
  
  for(i in 1:nclass){
    subset <- which(cids == i)
    sm[i,] <- colMeans(imgn[subset,])
    
  }
  
  
  cids_mat = cids
  dim(cids_mat) <- c(dim(image)[1],dim(image)[2])
  cids_cpp = cids_mat - 1
  q_img <- get_qimg_cpp(image, cids_cpp, sm)
  q_img
  
  
}




lkmeans <- function(image,nclass){
  
  result = list()
  
  data = image
  width = dim(data)[1] ; height = dim(data)[2]
  dim(data) = c(width*height,3)
  n = nrow(data)
  nsub = ceiling(n*0.01)
  s = sample(1:n,nsub)
  
  #Demo codes for a stratified sample with columns as strata
  # rowsamplesize = ceiling(height*0.01)
  # 
  # samples = list()
  # 
  # for(i in 1:width){
  #     samples[[i]] = sample(1:height,rowsamplesize) +(i-1)
  # }
  # 
  # 
  # s = do.call(rbind,samples) ; dim(s) = c(rowsamplesize*width)
  
  
  xtrain = data[s,]
  dummy = kmeans(xtrain,nclass,iter.max=1000,algorithm = "Lloyd")
  mu = dummy$centers
  
  library(clue)
  cids = cl_predict(dummy,data)
  
  
  result$class = as.numeric(cids)
  result$means = mu  
  return(result)
  
  
} 




ldakmeans <- function(image,nclass){
  
  result = list()
  
  
  
  data = image
  width = dim(data)[1] ; height = dim(data)[2]
  dim(data) = c(width*height,3)
  n = nrow(data)
  nsub = ceiling(n*0.01)
  s = sample(1:n,nsub)
  
  
  #Demo codes for a stratified sample with columns as strata
  # rowsamplesize = ceiling(height*0.01)
  # 
  # samples = list()
  # 
  # for(i in 1:width){
  #   samples[[i]] = sample(1:height,rowsamplesize) +(i-1)
  # }
  # 
  # 
  # s = do.call(rbind,samples) ; dim(s) = c(rowsamplesize*width)
 
  xtrain = data[s,]
  dummy = kmeans(xtrain,nclass,iter.max=1000,algorithm = "Lloyd")
  ytrain = as.factor(dummy$cluster)
  train = as.data.frame(cbind(ytrain,xtrain))
  colnames(train) = c("labels","x","y")
  library(MASS)
  fit = lda(labels~.,data=train,prior = rep(1/nclass,nclass))
  xtest = as.data.frame(data)
  colnames(xtest) = c("x","y")
  ypreds = predict(fit,xtest)$class
  mu = fit$means
  
  result$class = as.numeric(ypreds)
  result$centers = mu
  
  return(result)
}

