### Creates z-values
#'@export
zfunc<-function(x,y){
  (mean(x[y==1])-mean(x[y==0]))/stats::sd(x)
}
