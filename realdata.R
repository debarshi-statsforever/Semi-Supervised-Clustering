library(tiff)
library(pals) 
library(Rcpp)
library(jpeg)
library(png)
sourceCpp("calc_parallel.cpp") 
source("datafunctions.R")

nclass = 16


load_data <- function() {
  #x0 <- readTIFF(source="durga.tiff")
  x0 <- readJPEG(source="galaxy.jpg")
  shift <- 1
  newden <- 257
  newden <- max(255 + shift + 1, newden) 
  x1 <- x0 + (shift/255)
  x2 <- x1 * (255/newden)
  x2
}



view_image <- function(img) {
  dev.new()
  grid::grid.raster(img)
}

image <- load_data()

image1 <- image ; p = dim(image)[1]*dim(image)[2]
dim(image1) <- c(p,dim(image)[3])
dim(image1)


#Doing kmeans

start_time <- Sys.time()
kmfit = kmeans(image1, nclass, iter.max=1000, algorithm="Lloyd")
cids_km = kmfit$cluster
end_time <- Sys.time()

end_time - start_time

#Doing lkmeans

start_time <- Sys.time()
lkmfit = lkmeans(image,nclass)
cids_lkm = lkmfit$class
end_time <- Sys.time()

end_time - start_time



#Doing ldakmeans


start_time <- Sys.time()
ldakmfit = ldakmeans(image,nclass)
cids_ldakm = ldakmfit$class
end_time <- Sys.time()

end_time - start_time



quant_km = get_quantized_image(image,cids_km,16)
quant_lkm = get_quantized_image(image,cids_lkm,16)
quant_ldakm = get_quantized_image(image,cids_ldakm,16)

view_image(quant_km)
view_image(quant_lkm) 
view_image(quant_ldakm) 



