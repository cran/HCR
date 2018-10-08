#' @title Simulate the data of hidden compact representation model.
#' @description Generate the X->Y pair HCR data
#' @param sample_size Sample size
#' @param min_nx The minimum value of |X| (Default: 3)
#' @param max_nx The maximum value of |X| (Default: 15)
#' @param min_ny The minimum value of |Y| (Default: 3)
#' @param max_ny The maximum value of |Y| (Default: 15)
#' @param type type=0: standard version, type=1: |X|=|Y|, type=2: |Y'|=|Y|, type=3: |X|=|Y'|, type=4: |X|=|Y'|=|Y| (Default: type=0)
#' @param distribution The distribution of the cause X. The options are "multinomial","geom","hyper","nbinom","pois". Default: multinomial
#'
#' @return
#' return the synthetic data
#' @export
#'
#' @examples
#'df=simuXY(sample_size=100,type=0)
#'length(unique(df[,1]))
#'length(unique(df[,2]))
#'length(unique(df[,3]))
#'
#'df=simuXY(sample_size=100,type=1)
#'length(unique(df[,1]))
#'length(unique(df[,3]))
#'
#'df=simuXY(sample_size=100,type=2)
#'length(unique(df[,2]))
#'length(unique(df[,3]))
#'
#'df=simuXY(sample_size=100,type=3)
#'length(unique(df[,1]))
#'length(unique(df[,2]))
#'
#'df=simuXY(sample_size=100,type=4)
#'length(unique(df[,1]))
#'length(unique(df[,2]))
#'length(unique(df[,3]))
#'
simuXY<-function(sample_size=2000,min_nx=3,max_nx=15,min_ny=3,max_ny=15,type=0,distribution="multinomial"){
  Yp=NULL
  .=NULL
  if(is.null(distribution)||distribution=="multinomial"){
    setx<-1:sample(min_nx:max_nx,1)
    probx=abs(rnorm(length(setx)))
    probx=probx/sum(probx)
    X=sample(setx,sample_size,replace=T,prob = probx)
  }else if(distribution=="geom"){
    p=runif(1,0.4,0.7)
    X=rgeom(sample_size,p)
  }else if(distribution=="hyper"){
    M=sample(1:40,1)
    K=sample(1:M,1)
    N=sample(1:K,1)
    X=rhyper(sample_size,m = M,n = N,k = K)
  }else if(distribution=="nbinom"){
    R=sample.int(10,1)
    P=runif(1,0.8,0.9)
    X=rnbinom(sample_size,size=R,prob = P)
  }else if(distribution=="pois"){
    lambda=8*runif(1)
    X=rpois(sample_size,lambda)
  }else if(distribution=="uniform"){
    setx<-1:sample(min_nx:max_nx,1)
    X=sample(setx,sample_size,replace=T)
  }
  X=as.integer(factor(X))


  setx<-unique(X) # It is possible that the state of X will not be included.
  nx=length(setx)
  X=data.table(X,key="X")

  if(type==0){
    # type=0 standard version
    sety<-1:round(sample(min_ny:max_ny,1))
    ny=length(sety)
    setyp<-1:(sample(3:(min(nx,ny)),1))

    for(i in unique(X$X)){
      X[.(i),Yp:=sample(setyp,1)]
    }
    X$Yp<-as.integer(as.factor(X$Yp))
    setyp=unique(X$Yp)
    nyp=length(setyp)

  }else if(type==1){
    # type=1 nx=ny
    sety<-1:nx
    ny=length(sety)
    setyp<-1:round(sample(2:ny,1))

    for(i in unique(X$X)){
      X[.(i),Yp:=sample(setyp,1)]
    }
    X$Yp<-as.integer(as.factor(X$Yp))
    setyp=unique(X$Yp)
    nyp=length(setyp)

  }else if(type==2){
    #type=2 yp=y

    setyp<-1:nx
    for(i in unique(X$X)){
      X[.(i),Yp:=sample(setyp,1)]
    }
    X$Yp<-as.integer(as.factor(X$Yp))
    setyp=unique(X$Yp)
    nyp=length(setyp)
    sety<-1:nyp
    ny=nyp
  }else if(type==3){
    #type 3 nx=nyp
    X$Yp<-X$X
    setyp=unique(X$Yp)
    nyp=length(setyp)
    sety<-1:round(sample(nyp:max_ny,1))
    ny=length(sety)
  }else if(type==4){
    #type 4 nx=nyp=ny
    X$Yp<-X$X
    setyp=unique(X$Yp)
    nyp=length(setyp)
    sety<-1:nyp
    ny=nyp
  }
  A<-matrix(rnorm(nyp*ny),nrow=nyp,ncol=ny)
  A<-abs(A)
  A<-A/rowSums(A)
  rownames(A)=as.character(unique(X$Yp))
  colnames(A)=as.character(sety)
  Y=transfer(X$Yp,A)
  return(data.frame(X=X$X,Yp=X$Yp,Y=Y))
}

transfer<-function(data,A){
  # Transfer the data by a probaility table A
  if(uniqueN(data)>nrow(A)){
    stop("the row of matrix A shoule be smaller than the unique item of data")
  }
  data<-as.factor(data)
  data<-as.integer(data)
  term<-unique(data)
  data2<-data
  for(t in term){
    len=length(data[data==t])
    p=cumsum(A[as.integer(t),])
    r<-runif(len)
    tran<-c()
    for(j in 1:len){
      tran<-c(tran,sum(r[j]>p)+1)
    }
    data2[data==t]<-tran
  }
  return(data2)
}

transfer.fast<-function(X,A){
  Y=NULL
  r=NULL
  if(uniqueN(X)>nrow(A)){
    stop("the row of matrix A shoule be smaller than the unique item of data")
  }
  data<-data.table(X=X)
  term=unique(data$X)
  for(t in term){
    p=cumsum(A[as.character(t),])
    data[X==t,r:=runif(.N)]
    data[X==t,Y:=1]
    for(j in 1:(length(p)-1)){
      data[X==t&r>p[j],Y:=j+1]
    }

  }
  return(data$Y)
}

log_gamma<-function(x){
  return(sum(log(1:(x-1))))
}

Z_dirichlet<-function(n,alpha){
  if(length(n)!=length(alpha)){
    stop("n shoule be equal to alpha")
  }
  result=log_gamma(sum(n)+sum(alpha))
  for(i in 1:length(n)){
    result=result-log_gamma(n[i]+alpha[i])
  }
  return(result)
}




randomshuffle<-function(x){
  if(!is.factor(x)){
    c=class(x)
    x=as.factor(x)
    levels(x)<-sample(levels(x))
    return(as(as.character(x),c))
  }else{
    levels(x)<-sample(levels(x))
    return(x)
  }
}

roll <- function( x , n ){
  if( n == 0 )
    return( x )
  c( tail(x,n) , head(x,-n) )
}

scoretype<-function(score,d,m,score_type="bic"){
  if(score_type=="lik"){
    score=score  #lik
  }else if(score_type=="bic"){
    score=score-d/2*log(m)  #bic
  }else if(score_type=="aic"){
    score=score-d  #aic
  }else if(score_type=="aicc"){
    score=score-d-(d*(d+1)/(m-d-1))#aicc
  }
  return(score)
}

