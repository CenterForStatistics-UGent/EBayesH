#' This function computes the misclassification error \code{alpha}.
#'
#' \code{alpha} returns  alpha_pos and alpha_neg and their values after correlation
#' correction
#'
#' @param Delta  Empirical Bayes  effect size estimates from \code{EBayes}
#' @param tauSq  Empirical Bayes variance estimates from \code{EBayes}
#' @param Y  list of response variables as explained in \code{EBayes}
#' @param X  A merge of the orginal predictors values used in computing the z-values.
#' The  dataset with response variable as first entry of the list \code{Y} stays ontop
#' and the second dataset is that with response variable as second entry of the list \code{Y} and so on.
#' @param num_Pred is the number of predictor from which the prediction error is computed.
#' @export

alpha<-function(Delta,tauSq,Y,X,num_Pred)
{
  len<-length(X)
  beta_Delta<-((Delta)^2)*(1/tauSq)

  ord<-rev(order(abs(Delta)))


  beta_Delta<-beta_Delta[ord]
  Delta<-Delta[ord]
  tauSq<-tauSq[ord]

  logOdds<-sapply(Y,function(x)log(table(x)[2]/table(x)[1]))
  all_X<-NULL
  for(j in 1:len)
    all_X<-rbind(all_X,t(X[[j]])[,ord[1:num_Pred]])



  alpha_naive_pos<-matrix(rep(NA,num_Pred*len),ncol = len)
  alpha_cor_pos<-matrix(rep(NA,num_Pred*len),ncol = len)

  alpha_naive_neg<-matrix(rep(NA,num_Pred*len),ncol = len)
  alpha_cor_neg<-matrix(rep(NA,num_Pred*len),ncol = len)




  rhoC<-rep(NA,num_Pred)

  for(amt in c(1:num_Pred))
  {
    rhos<-rep(NA,amt-1)
    for(ff in 1:(amt-1))
    {
      if(amt==1)
      {
        rhos[ff]<-0
      }
      else
      {
        cors<-apply(as.matrix(all_X[,1:amt][,-c(1:ff)]),2,stats::cor,all_X[,ff])
        varTerm<-sqrt(1/tauSq[ff])*sqrt(1/(tauSq[1:amt][-c(1:ff)]))
        betaTerm<-Delta[ff]*Delta[1:amt][-c(1:ff)]
        rhos[ff]<-sum(betaTerm*varTerm*cors)
      }

    }
    rhoC[amt]<-2*sum(rhos)

    den_cor<-sqrt(sum(beta_Delta[1:amt])+2*sum(rhos))
    den_naive<-sqrt(sum(beta_Delta[1:amt]))
    num<-0.5*sum(beta_Delta[1:amt])


    alpha_cor_pos[amt,]<-stats::pnorm((-logOdds-num)/den_cor)
    alpha_cor_neg[amt,]<-stats::pnorm((logOdds-num)/den_cor)

    alpha_naive_pos[amt,]<-stats::pnorm((-logOdds-num)/den_naive)
    alpha_naive_neg[amt,]<-stats::pnorm((logOdds-num)/den_naive)

  }
  list(alpha_cor_pos=alpha_cor_pos,
       alpha_cor_neg=alpha_cor_neg,
       alpha_naive_pos=alpha_naive_pos,
       alpha_naive_neg=alpha_naive_neg,
       rho=rhoC)
}
