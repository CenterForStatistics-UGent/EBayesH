#' This function computes the misclassification error rataes \code{alpha_+} and  \code{alpha_-}.
#'
#' \code{alpha} returns  naive estimates of \code{alpha_+} and \code{apha_-} \code{alpha_naive_pos} and
#' \code{alpha_naive_neg}
#' and their values after correlation correction.
#'
#' @param beta  Empirical Bayes estimates of the  overall effect variance ratio from \code{EBayes}.
#' @param tauSq  Empirical Bayes variance estimates from \code{EBayes}.
#' @param Y  list of response variables as explained in \code{EBayes}.
#' @param X  A merge of the orginal predictors values
#' used in computing the z-values.
#' The  dataset with response variable as first entry
#' of the list \code{Y} stays ontop
#' and the second dataset is that with response variable as
#' second entry of the list \code{Y} and so on.
#' @param num_Pred is the number of predictor from which the
#' prediction error is computed.
#' @examples
#' \dontrun{
#'### Download the datasets from GEO####
#'
#'gseNo.<-"GSE21374"
#'gset <- getGEO(gseNo., GSEMatrix =TRUE, getGPL=FALSE)
#'if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
#'gset <- gset[[idx]]
#'
#'file<-paste(getwd(),"gset",gseNo.,".rda",sep="")
#'save(gset,file=file)
#'
#'gseNo.<-"GSE36059"
#'gset <- getGEO(gseNo., GSEMatrix =TRUE, getGPL=FALSE)
#'if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
#'gset <- gset[[idx]]
#'
#'file<-paste(getwd(),"gset",gseNo.,".rda",sep="")
#'save(gset,file=file)
#'
#'gseNo.<-"GSE48581"
#'gset <- getGEO(gseNo., GSEMatrix =TRUE, getGPL=FALSE)
#'if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
#'gset <- gset[[idx]]
#'
#'file<-paste(getwd(),"gset",gseNo.,".rda",sep="")
#'save(gset,file=file)
#'
#'gseNo.<-"GSE50058"
#'gset <- getGEO(gseNo., GSEMatrix =TRUE, getGPL=FALSE)
#'if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
#'gset <- gset[[idx]]
#'
#'file<-paste(getwd(),"gset",gseNo.,".rda",sep="")
#'save(gset,file=file)
#'
#'
#'######################## READ IN DATA ###################################
#'## Read  in the 4 downloaded datasets
#'
#' #Dataset name
#' gseNo.<-"GSE21374"
#' #location
#' file<-paste(getwd(),"/gset",gseNo.,".rda",sep="")
#' #load data
#' assign(paste("gset",gseNo.,sep=""), get(load(file=file)))
#'
#' gseNo.<-"GSE36059"
#' file<-paste(getwd(),"/gset",gseNo.,".rda",sep="")
#' assign(paste("gset",gseNo.,sep=""), get(load(file=file)))
#'
#' gseNo.<-"GSE48581"
#' file<-paste(getwd(),"/gset",gseNo.,".rda",sep="")
#' assign(paste("gset",gseNo.,sep=""), get(load(file=file)))
#'
#' gseNo.<-"GSE50058"
#' file<-paste(getwd(),"/gset",gseNo.,".rda",sep="")
#' assign(paste("gset",gseNo.,sep=""), get(load(file=file)))
#'
#'#############################################################################
#'
#'
#'
#'
#'## Note that different studies name things differently some
#'## classify patients as reject and non-reject while others
#'## subdivide the reject group into  tcmr, abmr and mixed.
#'## Patients that had nephrectomy are not eligible and  are
#'## always deleted.
#'
#'
#' ### Get pnenotype (kidney rejection status) for GSE21374
#' phenoGSE21374<-phenoData(gsetGSE21374)
#' #patient ids
#' rowNames<-rownames(phenoGSE21374@data)
#' #See the reject non-reject proportions
#' table(phenoGSE21374@data[,"characteristics_ch1.3"])
#' #Code rejection status as 0 and 1
#' response<-ifelse(phenoGSE21374@data[,"characteristics_ch1.3"]=="rejection/non rejection: nonrej",0,1)
#' GSE21374Res<-data.frame(rowNames,response)#Patients and their rejection status
#'
#'
#' #For GSE21374
#' phenoGSE36059<-phenoData(gsetGSE36059)
#' rowNames<-rownames(phenoGSE36059@data)
#' table(phenoGSE36059@data[,"characteristics_ch1"])
#' #These patients are not eligible so we take them out
#' removeGSE36059<-c(1:length(rowNames))[phenoGSE36059@data[,"characteristics_ch1"]=="diagnosis: Nephrectomy"]
#' newd=phenoGSE36059@data[,"characteristics_ch1"][-removeGSE36059]
#' response<-ifelse(newd== "diagnosis: non-rejecting",0,1)
#' table(response)
#' rowN<-rowNames[-removeGSE36059]
#' GSE36059Res<-data.frame(rowN,response)
#'
#'
#'
#' phenoGSE48581<-phenoData(gsetGSE48581)
#' rowNames<-rownames(phenoGSE48581@data)
#' table(phenoGSE48581@data[,"characteristics_ch1.1"])
#' #These patients are not eligible so we take them out
#' removeGSE48581<-c(1:length(rowNames))[phenoGSE48581@data[,"characteristics_ch1.1"]
#'                                       =="diagnosis (tcmr, abmr, mixed, non-rejecting, nephrectomy): nephrectomy"]
#' newd=phenoGSE48581@data[,"characteristics_ch1.1"][-removeGSE48581]
#' response<-ifelse(newd== "diagnosis (tcmr, abmr, mixed, non-rejecting, nephrectomy): non-rejecting",0,1)
#' table(response)
#' rowN<-rowNames[-removeGSE48581]
#' GSE48581Res<-data.frame(rowN,response)
#' sum(table(GSE48581Res[,2]))
#'
#'
#' phenoGSE50058<-phenoData(gsetGSE50058)
#' rowNames<-rownames(phenoGSE50058@data)
#' table(phenoGSE50058@data[,"characteristics_ch1"])
#' response<-ifelse(phenoGSE50058@data[,"characteristics_ch1"]=="patient group: stable patient (STA)",0,1)
#' GSE50058Res<-data.frame(rowNames,response)
#' table(GSE50058Res[,2])
#'
#'
#' #### Get expression data
#' X_GSE21374<-exprs(gsetGSE21374)
#' X_GSE36059<-exprs(gsetGSE36059)
#' X_GSE48581<-exprs(gsetGSE48581)
#' X_GSE50058<-exprs(gsetGSE50058)
#'
#' ### log2 transformations
#' X_GSE21374=t(apply(X_GSE21374,1, transFunc))
#' X_GSE36059=t(apply(X_GSE36059,1, transFunc))
#' X_GSE48581=t(apply(X_GSE48581,1, transFunc))
#' X_GSE50058=t(apply(X_GSE50058,1, transFunc))
#'
#' #Removing ineligible patients
#' X_GSE36059<-X_GSE36059[,-removeGSE36059]
#' X_GSE48581<-X_GSE48581[,-removeGSE48581]
#'
#'
#'
#' z_GSE50058<-apply(X_GSE50058,1,zfunc,GSE50058Res[,2])
#' nas<-c(1:length(z_GSE50058))[is.na(z_GSE50058)]
#' ### Removing bad genes
#'
#' X_GSE21374<-X_GSE21374[-nas,]
#' X_GSE36059<-X_GSE36059[-nas,]
#' X_GSE48581<-X_GSE48581[-nas,]
#' X_GSE50058<-X_GSE50058[-nas,]
#'
#' ### list of datasets and their responses
#' X<-list(X_GSE21374,
#'         X_GSE36059,
#'         X_GSE48581,
#'        X_GSE50058)
#'
#' Y<-list(GSE21374Res[,2],
#'         GSE36059Res[,2],
#'         GSE48581Res[,2],
#'         GSE50058Res[,2])
#' ### row and column scaling, two iterations
#' for(i in 1:2)
#' {
#'   X<-lapply(X,function(x)
#'   {
#'     x<-scale(x)
#'     x<-scale(t(x))
#'     t(x)
#'   })
#' }
#' ### z-values
#' Z<-mapply(function(x,y){apply(x,1,zfunc,y)},X,Y)
#' res<-EBayes(Z, Y, df = 7, breaks = 120)
#' alpha(res$EB_Delta,res$EB_tauSq,Y,X,170)
#' }
#' @return
#' \item{alpha_cor_pos}{correlation correted estimate of alpha_+ for each predictor
#' added to the linear classifier }
#' \item{alpha_cor_neg}{correlation correcd estimate of alpha_- for each predictor
#' added to the linear classifier }
#' \item{alpha_naive_pos}{naive estimate of alpha_+ for each predictor added to the
#' linear classifier }
#' \item{alpha_naive_neg}{naive estimate of alpha_- for each predictor added to the
#' linear classifier }
#' @details
#' Starting from the linear classifier, \code{D}.
#' The binary outcome Y is coded as +1 and 1.
#' The prediction rule is if D>1 classify subject to
#' the group +1 otherwise classify subject to the group -1.
#' The misclassification error rates \code{alpha} is
#' made up of two parts \code{alpha_+} and \code{alpha_-}
#' for which \code{alpha}=\code{alpha_+Pr(Y=+1) + alpha_-Pr(Y=-1)}.
#' This function computes estimates of alpha_+ and alpha_- ;
#' both the estimated correlated corrected versions of
#' \code{alpha_+} and \code{alpha_-} are returned as
#' \code{alpha_cor_pos} and \code{alpha_cor_neg}, respectively.
#' Their naive versions are also returned \code{alpha_niave_pos}
#' and \code{alpha_naive_neg}.
#' @export

alpha<-function(beta,tauSq,Y,X,num_Pred)
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
