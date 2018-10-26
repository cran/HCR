L<-function(X,YPrime,Y,score_type="bic",is_anm=FALSE,is_cyclic=FALSE,...){

  freqx=as.data.frame(table(X))$Freq
  freqx<-freqx[freqx!=0]
  px= freqx/sum(freqx)
  nx<-uniqueN(X)
  nyp<-uniqueN(YPrime)

  if (is_anm) {
    # data preprocessing for ANM
    if(is_cyclic){
      Y=(Y-YPrime)%%(max(Y)-min(Y)+1)
    }else{
      Y=Y-YPrime
    }
    tab=table(Y) # |Y'|=1
    nyp=1 # |Y'|=1
    freq=as.data.frame(tab)$Freq
    B=as.data.frame(tab/sum(tab))$Freq # |Y'|=1
  }else{
    tab=table(YPrime,Y)
    freq=as.data.frame(tab)$Freq
    B=as.data.frame(tab/rowSums(tab))$Freq
  }

  # remove the categories that have zero frequency.
  sel<-(freq!=0)
  B=B[sel]
  freq=freq[sel]

  #d=nyp*(uniqueN(Y)-1)+(nx-1)+nx
  d=nyp*(uniqueN(Y)-1)+nx-1
  if(score_type=="log"){
    return(sum(freq*log(B))+sum(freqx*log(px))) #log-likelihood
  }else if(score_type=="bic"){
    return(sum(freq*log(B))+sum(freqx*log(px))-d/2*log(length(X))) #bic
  }else if(score_type=="aic"){
    return(sum(freq*log(B))+sum(freqx*log(px))-d) #aic
  }else if(score_type=="aicc"){
    return(sum(freq*log(B))+sum(freqx*log(px))-d-(d*(d+1)/(length(X)-d-1))) #aicc
  }

}


#' @title Hidden Compact Representation Model
#' @description Causal Discovery from Discrete Data using Hidden Compact Representation.
#' @param X The data of cause.
#' @param Y The data of effect.
#' @param score_type You can choose "bic","aic","aicc","log" as the type of score to fit the HCR model. Default: bic
#' @param is_anm If is_anm=TRUE, it will enable a data preprocessing to adjust for the additive noise model.
#' @param is_cyclic If is_anm=TRUE and is_cyclic=TRUE, it will enable a data preprocessing to adjust the cyclic additive noise model.
#' @param verbose Show the score at each iteration.
#' @param max_iteration The maximum iteration.
#' @param ... Other arguments passed on to methods. Not currently used.
#' @return The fitted HCR model and its score.
#' @export
#' @examples
#' library(data.table)
#' set.seed(10)
#' data=simuXY(sample_size=200)
#' r1<-HCR(data$X,data$Y)
#' r2<-HCR(data$Y,data$X)
#' # The canonical hidden representation
#' unique(r1$data[,c("X","Yp")])
#' # The recovery of hidden representation
#' unique(data.frame(data$X,data$Yp))
#'
HCR<-function(X,Y,score_type="bic",is_anm=FALSE,is_cyclic=FALSE,verbose=FALSE,max_iteration=1000,...){

  Yp=NULL
  .=NULL
  setx=unique(X)
  sety=unique(Y)
  nx=length(setx)
  ny=length(sety)

  dt=data.table(X,Y,key="X")



  if(is_anm){
    if(is.character(Y)||is.factor(Y)||is.numeric(Y)){
      Y=as.integer(factor(Y,labels=1:ny))
      dt=data.table(X,Y,key="X")
      setx=unique(dt$Y)
    }
    if(is.character(X)||is.factor(X)||is.numeric(X)){
      X=as.integer(factor(X,labels=1:nx))
      dt=data.table(X,Y,key="X")
      setx=unique(dt$X)
    }
    setyp=min(sety):max(sety)
  }else{
    setyp=1:nx #random assign :default setting
  }
  for(i in setx){
    temp=dt[.(i),.N,by=Y]
    dt[.(i),Yp:=temp$Y[which.max(temp$N)]] #set Yp to the most freq Y
  }

  if(!is_anm){
    # rearrange Y'
    dt$Yp<-factor(dt$Yp,labels=1:uniqueN(dt$Yp))
    dt$Yp<-as.integer(dt$Yp)
  }


  newScore=L(X=dt$X,YPrime=dt$Yp,Y=dt$Y,score_type=score_type,is_anm=is_anm,is_cyclic=is_cyclic,...)


  bestdt=dt
  iteration=0L
  oldScore=-Inf
  while(newScore>oldScore){
    if(iteration>max_iteration){
      break
    }

    if(verbose){
      #if(iteration%%10==0){
        print(sprintf("iteration:%d  score:%f",iteration,newScore))
      #}
    }
    iteration=iteration+1L
    oldScore=newScore

    for(i in setx){
      for(j in setyp){
        temp<-copy(dt)
        temp[.(i),Yp:=j]

        score<-L(X=temp$X,YPrime=temp$Yp,Y=temp$Y,score_type=score_type,is_anm=is_anm,is_cyclic=is_cyclic,...)
        if(score>newScore){
          newScore=score
          bestdt=copy(temp)
        }
      }
    }
    dt=bestdt
  }
  newScore<-oldScore

  if(verbose){
      print(sprintf("iteration:%d  score:%f",iteration,newScore))
  }
  return(list(data=bestdt,score=newScore))
}




#' @title The Fast Version for Fitting Hidden Compact Representation Model
#' @description A fast implementation for fitting the HCR model.
#' This implementation caches all intermediate results to speed up the greedy search.
#' The basic idea is that if there are two categories need to be combined, for instance, X=1 and X=2 mapping to the same Y'=1, then the change of the score only depend on the frequency of the data where X=1 and X=2.
#' Therefore, after combination, if the increment of the likelihood is greater than the penalty, then we will admit such combination.
#' @param X The data of cause.
#' @param Y The data of effect.
#' @param score_type You can choose "bic","aic","aicc","log" as the type of score to fit the HCR model. Default: bic
#' @param ... Other arguments passed on to methods. Not currently used.
#' @return The fitted HCR model and its score.
#' @export
#' @examples
#' library(data.table)
#' set.seed(1)
#' data=simuXY(sample_size=2000)
#' r1=HCR.fast(data$X,data$Y)
#' r2=HCR.fast(data$Y,data$X)
#' # The canonical hidden representation
#' unique(r1$data[,c("X","Yp")])
#' # The recovery of hidden representation
#' unique(data.frame(data$X,data$Yp))
#'
HCR.fast<-function(X,Y,score_type="bic",...){
  Yp=NULL
  .=NULL
  X=as.integer(factor(X))
  Y=as.integer(factor(Y))
  setx=unique(X)
  sety=unique(Y)
  nx=length(setx)
  ny=length(sety)
  m=length(X)
  dt=data.table(X,Y,key="X")
  dt[,Yp:=X] #set Yp to the most freq Y
  nyp=nx
  N=table(X,Y)
  Nx=rowSums(N)
  pxy=N/Nx

  pxy[pxy==0]<-1
  likxy<-rowSums(N*log(pxy)) #likelihood on P(Y|X=x)
  besteps=0
  repeat{
    # Combine best ith, jth categories at each iteration.
    iscontinue=F
    # d   =nyp*(ny-1)+(nx-1)+nx
    # dnew=(nyp-1)*(ny-1)+(nx-1)+nx
    d   =nyp*(ny-1)+nx-1
    dnew=(nyp-1)*(ny-1)+nx-1
    if (score_type=="bic") {
      penalty1=d/2*log(m)
      penalty2=dnew/2*log(m)
    }else if(score_type=="aic"){
      penalty1=d
      penalty2=dnew
    }else if(score_type=="log"){
      penalty1=0
      penalty2=0
    }else if (score_type=="aicc") {
      penalty1=d+(d*(d+1)/(m-d-1))
      penalty2=dnew+(dnew*(dnew+1)/(m-dnew-1))
    }
    #rowx=as.integer(rownames(N))
    for(i in 1:(nyp-1)){
      for (j in (i+1):nyp) {
        nij=N[i,]+N[j,]
        sel=nij!=0
        combine_lik=sum(nij[sel]*log(nij[sel]/(Nx[i]+Nx[j])))
        eps=combine_lik-penalty2-likxy[i]-likxy[j]+penalty1

        if (eps>=besteps) {
          besteps=eps
          besti=i
          bestj=j
          bestrownamesi=as.integer(rownames(N)[i])
          bestrownamesj=as.integer(rownames(N)[j])
          bestlik=combine_lik
          iscontinue=T
        }
      }
    }
    if(!iscontinue){
      break
    }
    besteps=0
    nyp=nyp-1
    N[besti,]=N[besti,]+N[bestj,]
    N=N[-bestj,]
    likxy[besti]=bestlik
    likxy=likxy[-bestj]
    dt[Yp==bestrownamesj,Yp:=dt[.(bestrownamesi),"Yp"][1]]


    if (is.null(nrow(N))) {
      Nx=sum(N)
      break
    }else{
      Nx=rowSums(N)
    }
  }

  newScore=L(X=dt$X,YPrime=dt$Yp,Y=dt$Y,score_type=score_type,...)
  return(list(data=dt,score=newScore))
}


