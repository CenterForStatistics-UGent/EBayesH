### Creates u-values
ufunc<-function(x)
{
  (x-mean(x))/stats::sd(x)
}
