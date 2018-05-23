#' This function computes empirical Bayes estimates of effect sizes
#' and their variances across studies.
#'
#' \code{EBayes}  returns empirical Bayes estimates of effect sizes and their variances.
#'
#' @param Z is a matrix with columns containing the z-values for each predictor  from the different studies.
#' the dimension of Z should be "number of predictors by number of studies"
#' @param Y is a list containing the binary response variables for each study. Each element of the list
#' should be of length equal to the number of obsersation in each study.
#'  The order of the z-values in the columns of  \code{Z} must be maintained.
#'  That is if the z-values of study  \code{k} are in column  \code{k} of \code{Z},
#'  then response values of study  \code{k} should be the kth element of the list  \code{Y}.
#'@param df is the degrees of freedom of the spline fit used to non-parametrically
#'estimate the marginal density (see Efron (2009)).
#'The default value is  df=7. This will work in  most cases.
#'@param breaks are the number of bins used in the spline fit (see   Efron (2009)).
#'The default is 120 and will also work in most cases if
#'the number of predictors is not below 1000.
#'@details
#'Details are found in Mbah et al. (2018).
#'With multiple high-dimensional datasets,
#'this function computes empirical Bayes
#'estimates of each predictor's overall
#'effect size \code{beta}
#'which is the ratio between the predictor's
#'mean and its between studies variance.  It also computes the predictor variances across studies (heterogeneity).
#'The predictor effect sizes can then be used in a linear prediction rule.
#'@return
#'\item{EB_beta}{a vector of empirical Bayes
#'estimates of  predictor effects}
#'\item{EB_tauSq}{a vector of variances across studies for each predictor.}
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
  EB_beta<-(Deltahat*EB_tauSq)+(1/c_s)*rowSums(dlogm)

  list(EB_beta=EB_beta,EB_tauSq=EB_tauSq)
}
