logmhat<-function(z,df=3,breaks=120){
  #estimating f
  bre<-seq(min(z),max(z),length=breaks)
  x<-(bre[-1]+bre[-length(bre)])/2

  zz<-pmin(pmax(z,min(bre)),max(bre))
  y<-graphics::hist(zz,breaks=bre,plot=F)$counts

  fit<-stats::glm(y~splines::ns(x,df=df),stats::poisson)$fit
  logm<-log(fit)

  logm.<-sfsmisc::D1ss(x, logm, z)
  x[1]<-bre[1]
  x[length(x)]=bre[length(bre)]
  #browser()
  fit<-stats::approx(x,fit,z)$y
  out<-list(z=z,m=fit/max(fit),logmhat=logm.)
  return(out)
}
