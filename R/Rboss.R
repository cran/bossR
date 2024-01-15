#'  Select Optimal cutoff for a biomarker
#' @description Given a set of data including survival time ,censor status and Biomarker values,
#'    return the Optimal cutoff for the biomarker.
#' @param data A data frame which contains 3 columns for cox regression : survival time, censor status, Biomarker values.
#'    2 columns for linear regression :  Y, X.
#' @param cutoff Numeric vector of candidate cutoffs.
#' @param type A number; if =1, will perform linear regression;if =2(default) will perform cox regerssion.
#'
#' @return Optimal cutoff for the biomarker, the FWER of the model
#' @export
#' @references BOSS - Biomarker Optimal Segmentation System
#' @importFrom survival coxph Surv
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats lm
#' @examples cutoff=c(56,112,167,223,278,334,389,445)
#'     data(myGene)
#'     getOC(data=myGene,cutoff)
getOC=function(data,cutoff,type=2){

  if(type==1){

    colnames(data)=c('y','X')
    data=data[order(-data$X),]
    fmod=lm(data$y~data$X)
    sig=summary(fmod)$var[1]$coefficients[2,2]
    n=nrow(data)
    len=length(cutoff)
    beta_hat=c()
    beta_hatsd=c()
    m=c()
    se=c()
    for(i in 1:length(cutoff)){
      point=cutoff[i]
      result=getbeta(data,point,type)
      beta_hat[i]=result[1]
      m[i]=result[2]
      se[i]=result[3]
    }

    sigma<-matrix(0,nrow=length(cutoff),ncol=length(cutoff))
    for(i in 1:length(cutoff)){
      for(j in 1:length(cutoff)){
        mmax=max(m[i],m[j])
        mmin=min(m[j],m[i])
        sigma[i,j]=n/(mmax*(n-mmin))*sqrt((n-m[i])*m[i]/n)*sqrt((n-m[j])*m[j]/n)
      }
    }

    for(i in 1:length(cutoff)){
      beta_hatsd[i]=beta_hat[i]/se[i]
    }


    Sigma=sigma
    beta=rep(0,length(cutoff))
    idk=which(abs(beta_hatsd)==max(abs(beta_hatsd)))
    mbetaksd=beta_hatsd[idk]
    mbetak=beta_hat[idk]
    pk=getpvalue(mbetaksd,beta,length(cutoff),Sigma)
    c=0


    return(data.frame('optimal cutoff'=cutoff[idk],'FWER'=pk[1]))

  }else if(type==2){

    colnames(data)=c('time','status','X')
    data=data[order(-data$X),]
    fmod=coxph(Surv(data$time,data$status)~data$X)
    sig=sqrt(fmod$var[1])
    n=nrow(data)
    len=length(cutoff)
    beta_hat=c()
    beta_hatsd=c()
    m=c()
    se=c()
    for(i in 1:length(cutoff)){
      point=cutoff[i]
      result=getbeta(data,point,type)
      beta_hat[i]=result[1]
      m[i]=result[2]
      se[i]=result[3]
    }

    sigma<-matrix(0,nrow=length(cutoff),ncol=length(cutoff))
    for(i in 1:length(cutoff)){
      for(j in 1:length(cutoff)){
        mmax=max(m[i],m[j])
        mmin=min(m[j],m[i])
        sigma[i,j]=n/(mmax*(n-mmin))*sqrt((n-m[i])*m[i]/n)*sqrt((n-m[j])*m[j]/n)
      }
    }

    for(i in 1:length(cutoff)){
      beta_hatsd[i]=beta_hat[i]/se[i]#sqrt(n*sig/(m[i]*(n-m[i])))
    }

    Sigma=sigma
    beta=rep(0,length(cutoff))
    idk=which(abs(beta_hatsd)==max(abs(beta_hatsd)))
    mbetaksd=beta_hatsd[idk]
    mbetak=beta_hat[idk]
    pk=getpvalue(mbetaksd,beta,length(cutoff),Sigma)



    return(data.frame('optimal cutoff'=cutoff[idk],'FWER'=pk[1]))
  }else{

    stop('error in input data type')
    }

}

#' Get regression coefficient
#' @description Computes the regression coefficient of certain regression based on certain cutoff.
#' @param data A data frame which contains 3 columns for cox regression : survival time, censor status, Biomarker values.
#'    2 columns for linear regression :  Y, X.
#' @param point A number to cut biomarker or X.
#' @param type A number; if =1, will perform linear regression; if =2(default) will perform cox regression.
#'
#' @return An object with 3 class: Coefficient beta, number of samples of which the biomarker is greater than the point, standard error of coefficient estimation.
#' @importFrom survival coxph Surv
#' @importFrom stats lm
getbeta=function(data,point,type=2){
  if(type==1){#1-lin 2-cox
    xi=c()
    xi[0:point]=rep(1,point)
    xi[(point+1):nrow(data)]=rep(0,nrow(data)-point)
    sdata=cbind(data,xi)
    mod=lm(sdata$y~xi)
    betai=mod$coefficients[2]
    sigma=summary(mod)$coefficients[2,2]
    mi=sum(xi)
    return(c(betai,mi,sigma))
  }else if(type==2){
    xi=c()
    xi[0:point]=rep(1,point)
    xi[(point+1):nrow(data)]=rep(0,nrow(data)-point)
    sdata=cbind(data,xi)
    mod=coxph(Surv(sdata$time,sdata$status)~xi)
    betai=mod$coefficients[1]
    sigma=sqrt(mod$var[1])
    mi=sum(xi)
    return(c(betai,mi,sigma))
  }else{
    stop('error in input data type ')
  }

}
#' Computes the distribution function of the multivariate normal distribution
#' @description Computes the distribution function of the multivariate normal distribution.
#' @param threshold A number.
#' @param mu The mean vector of length n.
#' @param n A number indicates dimension.
#' @param Sigma The correlation matrix of dimension n.
#' @return The evaluated distribution function
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats lm
getpvalue=function(threshold,mu,n,Sigma){
  if(threshold>0){
    p=1-pmvnorm(lower=-Inf,upper=rep(threshold,n),mean=mu, sigma=Sigma)
  }else{
    p=1-pmvnorm(lower=rep(threshold,n),upper=Inf,mean=mu, sigma=Sigma)
  }
  return(p)
}



