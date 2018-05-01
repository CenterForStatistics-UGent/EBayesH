### log2 transformation of the gene expressions data
transFunc=function(x){
  # log2 transform automatically
  qx <- as.numeric(stats::quantile(x, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {x[which(x <= 0)] <- NaN
  x <- log2(x) }
  return(x)
}
