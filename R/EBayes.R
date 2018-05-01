#' This function computes empirical Bayes estimates of effect sizes
#' and their variances across studies.
#'
#' \code{EBayes}  returns empirical Bayes estimates of effect sizes and their variances.
#'
#' @param Z  is a matrix with columns containing the z-values  from the different studies
#' @param Y is a list containing the binary response variable for each study.
#'  The order of the z-values in the columns of  \code{Z} must be maintained.
#'  That is if the z-values of study  \code{k} are in column  \code{k} of \code{Z},
#'  then response values of study  \code{k} should be the kth element of the list  \code{Y}.
#'@param df is the degrees of freedom of the spline fit used to non-parametrically estimate the marginal density. See
#'Efron (2009). We set the default to 7. This will work for most cases.
#'@param breaks are the number of bins used in the spline fit see   Efron (2009).
#'The default is 120 and will also work in most cases if
#'the number of predictors is not below 1000.
#'@examples
#'res <- EBayes(Z,Y)
#'head(res)
#' @export

EBayes<-function(Z,Y,df=7,breaks=120)
{
  c_sVec<-1/sapply(Y,function(x)sum(1/table(x)),simplify = TRUE)
  c_s<-sum(c_sVec)
  summand<-Z*c_sVec

  J<-length(c_sVec)
  Deltahat<-(1/c_s)*rowSums(summand)

  z_Delta_Sq<-apply(Z,2,function(x,y)(x-Deltahat)^2)

  temp<-z_Delta_Sq*c_sVec
  tauSqhat<-(1/(J-1))*rowSums(temp)

  dlogm<-apply(Z,2,function(x)logmhat(x,df=df,breaks = breaks)$logmhat)

  EB_tauSq<-abs((1/((J-1)*tauSqhat))*rowSums((Deltahat-Z)*dlogm))
  EB_Delta<-(Deltahat*EB_tauSq)+(1/c_s)*rowSums(dlogm)

  list(EB_Delta=EB_Delta,EB_tauSq=EB_tauSq)
}