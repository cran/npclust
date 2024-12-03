#' @title Nonparametric Clustered Data Analysis
#' @description Main function to calculate nonparametric effect sizes and conduct hypothesis tests.
#' @name ncda
#' @param formula An object of class "formula" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted. The details of model specification are
#' given under ‘Details’.
#' @param data a data frame in the long format.
#' @param period time indicator variable.
#' @param subject subject or cluster ID
#' @param indicator an optional vector of characters indicating the order of pre and post intervention period;
#' must match the levels of period argument if specified; if not specified, the pre and post intervention
#' period will be ordered in the alphabet order by default
#' @param Contrast an optional contrast matrix for effect sizes.
#'
#' @usage ncda(formula,data,period,subject,indicator=NULL,Contrast=NULL)
#'
#' @return An object with effect sizes and other test details.
#' @examples
#' \donttest{
#' ARTIS_analysis <- ncda(symptoms_pqol~tx, ARTIS, intervention, homeid,
#'                         indicator=c("0","1"),
#'                         Contrast=matrix(c(1,-1,1,-1,1,-1), nrow = 1))
#' names(ARTIS_analysis)
#' ARTIS_analysis$p.vector
#' }
#' skin_analysis <- ncda(score~tx, skin, intervention, subject,
#'                       indicator=c("control","treatment"),
#'                       Contrast=matrix(c(1,-1), nrow = 1))
#' skin_analysis$TotalSampleSize
#' skin_analysis$p.vector
#' @references Cui, Yue, Frank Konietschke, and Solomon W. Harrar. "The
#' nonparametric Behrens–Fisher problem in partially complete clustered data."
#' Biometrical Journal 63.1 (2021): 148-167.
#'
#' Harrar, Solomon W., and Yue Cui. "Nonparametric methods for clustered data
#' in pre-post intervention design." Journal of Statistical Planning and
#' Inference 222 (2023): 1-21.
#' @details The model has the form response ~ tx where response is the (numeric)
#' response variable and tx is the treatment variable.
#' @importFrom MASS ginv
#' @export

####repeated measurements within factorial design####

#modified - check if any of nc, n1 or n2 is 0
ncda <- function(formula,data,period,subject,indicator=NULL,Contrast=NULL){
  resp <- formula[[2]]
  tx <- formula[[3]]
  x <- eval(substitute(ProcessData(data,tx,period,subject,resp,indicator)))
  trt <- x$trt
  x1c <- x$x1c
  x2c <- x$x2c
  x1i <- x$x1i
  x2i <- x$x2i
  #store data into lists
  xc.list <-list(x1c,x2c)
  x1 <- lapply(1:trt, function(i){c(x1c[[i]],x1i[[i]])})
  x2 <- lapply(1:trt, function(i){c(x2c[[i]],x2i[[i]])})
  xi.list <- list(x1i,x2i)

  nc <- x$nc
  n1 <- x$n1
  n2 <- x$n2

  n.list <- list(n1,n2)

  #get cluster size
  m1c <- x$m1c
  m2c <- x$m2c
  m1i <- x$m1i
  m2i <- x$m2i

  mc.list <- list(m1c,m2c)
  mi.list <- list(m1i,m2i)

  #calculate sample size
  N1c <- sapply(1:trt, function(i){sum(m1c[[i]])})
  N2c <- sapply(1:trt, function(i){sum(m2c[[i]])})
  Nc.list <- list(N1c,N2c)

  N1i <- sapply(1:trt, function(i){sum(m1i[[i]])})
  N2i <- sapply(1:trt, function(i){sum(m2i[[i]])})
  Ni.list <- list(N1i,N2i)

  N1 <- sapply(1:trt, function(i){sum(m1c[[i]])+sum(m1i[[i]])})
  N2 <- sapply(1:trt, function(i){sum(m2c[[i]])+sum(m2i[[i]])})
  N <- N1+N2
  N.list <-lapply(1:trt, function(i){c(N1[i],N2[i])})
  #same n for each trt comb
  n <- nc+n1+n2

  Cal.w <- function(r,s,j,l){
    if((nc[r]==0)|(nc[j]==0)){sum1 <- 0}
    else{sum1 <- sum(sapply(1:nc[r], function(a){sum(sapply(1:nc[j],function(k){sum(sapply(1:mc.list[[l]][[j]][k],
                                     function(d){mc.list[[s]][[r]][a]*Fhat(xc.list[[l]][[j]][[k]][d],unlist(xc.list[[s]][[r]][[a]]))}))}))}))}
    if((nc[r]==0)|(n.list[[l]][j]==0)){sum2 <- 0}
    else{sum2 <- sum(sapply(1:nc[r], function(a){sum(sapply(1:n.list[[l]][j],function(k){sum(sapply(1:mi.list[[l]][[j]][k],
                                     function(d){mc.list[[s]][[r]][a]*Fhat(xi.list[[l]][[j]][[k]][d],unlist(xc.list[[s]][[r]][[a]]))}))}))}))}
    if((n.list[[s]][r]==0)|(nc[j]==0)){sum3 <- 0}
    else{sum3 <- sum(sapply(1:n.list[[s]][r], function(a){sum(sapply(1:nc[j],function(k){sum(sapply(1:mc.list[[l]][[j]][k],
                                              function(d){mi.list[[s]][[r]][a]*Fhat(xc.list[[l]][[j]][[k]][d],unlist(xi.list[[s]][[r]][[a]]))}))}))}))}
    if((n.list[[s]][r]==0)|(n.list[[l]][j]==0)){sum4 <- 0}
    else{sum4 <- sum(sapply(1:n.list[[s]][r], function(a){sum(sapply(1:n.list[[l]][j],function(k){sum(sapply(1:mi.list[[l]][[j]][k],
                                              function(d){mi.list[[s]][[r]][a]*Fhat(xi.list[[l]][[j]][[k]][d],unlist(xi.list[[s]][[r]][[a]]))}))}))}))}
    return(1/(N.list[[j]][l]*N.list[[r]][s])*(sum1+sum2+sum3+sum4))
  }


  EffSize <- function(j,l){
    return(sum(sapply(1:trt, function(r){sapply(1:2, function(s){Cal.w(r,s,j,l)})}))/(2*trt))
  }

  p.vector <- c(sapply(1:trt, function(j){c(EffSize(j,1), EffSize(j,2))}))

  w.all <- sapply(1:trt, function(i){sapply(1:2, function(s){sapply(1:trt, function(t){sapply(1:2,
                                                                                              function(l){Cal.w(i,s,t,l)})})})})

  W <- as.matrix(unlist(lapply(1:trt, function(i){w.all[1:2,i]})))

  for (a in seq(3,4*trt,2)) {
    W <- cbind(W, as.matrix(unlist(lapply(1:trt, function(i){w.all[a:(a+1),i]}))))
  }

  E <- t((2*trt)^(-1)*kronecker(rep(1,2*trt),diag(2*trt)))

  ######################################
  Fhat.estimate <- function(j,l,x){
    if(nc[j]==0){Fhatc <- 0} else{Fhatc <- sum(sapply(1:nc[j], function(k){mc.list[[l]][[j]][k]*Fhat(x,unlist(xc.list[[l]][[j]][[k]]))}))}
    if(n.list[[l]][j]==0){Fhati <- 0} else{Fhati <- sum(sapply(1:n.list[[l]][j], function(k){mi.list[[l]][[j]][k]*Fhat(x,unlist(xi.list[[l]][[j]][[k]]))}))}
    return((Fhatc+Fhati)/N.list[[j]][l])
  }

  Ybarc <- function(j,l,r,s,k){
    sum(sapply(1:mc.list[[l]][[j]][k], function(v){Fhat.estimate(r,s,xc.list[[l]][[j]][[k]][v])}))/mc.list[[l]][[j]][k]
  }

  Ybari <- function(j,l,r,s,k){
    sum(sapply(1:mi.list[[l]][[j]][k], function(v){Fhat.estimate(r,s,xi.list[[l]][[j]][[k]][v])}))/mi.list[[l]][[j]][k]
  }

  #########calculate S matrix##############
  #need to define C's under two different cases####
  tau <- function(r,s,l,p,q,p.p,q.q){
    if((nc[r]==0)|(nc[r]==1)){return(0)}
    else {covar <- 1/(nc[r]-1)*sum(sapply(1:nc[r], function(k){mc.list[[s]][[r]][k]*mc.list[[l]][[r]][k]*
        (Ybarc(r,s,p,q,k)-Cal.w(p,q,r,s))*(Ybarc(r,l,p.p,q.q,k)-Cal.w(p.p,q.q,r,l))}))
    return(nc[r]/(N.list[[r]][s]*N.list[[r]][l])*covar)}
  }

  eta <- function(r,s,l,p,q,p.p,q.q){
    if((n.list[[s]][r]==0)|(n.list[[s]][r]==1)){return(0)}
    else {covar <- 1/(n.list[[s]][r]-1)*sum(sapply(1:n.list[[s]][r], function(k){mi.list[[s]][[r]][k]^2*
        (Ybari(r,s,p,q,k)-Cal.w(p,q,r,s))*(Ybari(r,s,p.p,q.q,k)-Cal.w(p.p,q.q,r,s))}))
    return(n.list[[s]][r]/(N.list[[r]][s]^2)*covar)}
  }

  #####overall function to compute entries in S#########
  S.value.complete <- function(j,l,r,s,p,q,p.p,q.q){
    if((j==r)&(l==s)){
      C1 <- tau(p,q,q.q,r,s,r,s)
      C2 <- tau(p,q,s,r,s,p.p,q.q)
      C5 <- tau(r,s,q.q,p,q,r,s)
      C6 <- tau(r,s,s,p,q,p.p,q.q)

      if((p!=p.p)&(p!=r)&(p.p!=r)){return(C6)}
      if((p==p.p)&(p!=r)){return(C1+C6)}
      if((p!=p.p)&(p==r)&(q!=s)){return(-C2+C6)}
      if((p!=p.p)&(p.p==r)&(q.q!=s)){return(-C5+C6)}
      if((p==p.p)&(p==r)&(q!=s)&(q.q!=s)){return(C1-C2-C5+C6)}
      else {return(0)}
    }

    else{
      C1 <- tau(p,q,q.q,j,l,r,s)
      C2 <- tau(p,q,s,j,l,p.p,q.q)
      C5 <- tau(j,l,q.q,p,q,r,s)
      C6 <- tau(j,l,s,p,q,p.p,q.q)

      if((p!=p.p)&(p!=r)&(j!=p.p)&(j==r)){return(C6)}
      if((p==p.p)&(p!=r)&(j!=p.p)&(j!=r)){return(C1)}
      if((p!=p.p)&(p==r)&(j!=p.p)&(j!=r)){return(-C2)}
      if((p!=p.p)&(p!=r)&(j==p.p)&(j!=r)){return(-C5)}

      if((p==p.p)&(p==r)&(j!=r)&(q.q!=s)){return(C1-C2)}
      if((p==p.p)&(p==j)&(j!=r)&(q!=l)){return(C1-C5)}
      if((p==p.p)&(p.p!=j)&(j==r)&(l!=s)){return(C1+C6)}
      if((p!=p.p)&(p==r)&(j==r)&(q!=l)&(l!=s)){return(-C2+C6)}
      if((p!=p.p)&(p!=r)&(p.p==r)&(j==r)&(q.q!=s)&(l!=s)){return(-C5+C6)}
      if((p!=p.p)&(p==r)&(j==p.p)){return(-C2-C5)}

      if((p==p.p)&(p==r)&(p.p==j)&(q!=l)&(s!=q.q)){return(C1-C2-C5+C6)}#check main function q!=l and s!=q.q
      else {return(0)}
    }
  }

  S.value.incomplete <- function(j,l,r,s,p,q,p.p,q.q){
    if((j==r)&(l==s)){
      if(((p==p.p)&(q==q.q)&(r==p)&(s==q))|(q!=q.q)|(p!=p.p)){C11 <- 0}
      else{C11 <- eta(p,q,q.q,r,s,r,s)}

      if(((p==r)&(q==s))|((p.p==r)&(q.q==s))){C16 <- 0}
      else{C16 <- eta(r,s,s,p,q,p.p,q.q)}

      C12 <- C15 <- 0

      return(C11+C16-C12-C15)
    }

    if((j!=r)|(l!=s)){

      if((p!=p.p)|(q!=q.q)|((p==p.p)&(q==q.q)&(j==p)&(l==q))|((p==p.p)&(q==q.q)&(r==p.p)&(s==q.q))){C11 <- 0}
      else{C11 <- eta(p,q,q.q,j,l,r,s)}

      if((p!=r)|(q!=s)|((p==r)&(q==s)&(p==j)&(q==l))|((p==r)&(q==s)&(r==p.p)&(s==q.q))){C12 <- 0}
      else{C12 <- eta(p,q,s,j,l,p.p,q.q)}

      if((j!=p.p)|(l!=q.q)|((j==p.p)&(l==q.q)&(p==j)&(q==l))|((j==p.p)&(l==q.q)&(r==p.p)&(s==q.q))){C15 <- 0}
      else{C15 <- eta(j,l,q.q,p,q,r,s)}

      if((j!=r)|(l!=s)|((j==r)&(l==s)&(p==j)&(q==l))|((j==r)&(l==s)&(p.p==r)&(q.q==s))){C16 <- 0}
      else{C16 <- eta(j,l,s,p,q,p.p,q.q)}

      return(C11+C16-C12-C15)
    }
  }

  S.value <- function(j,l,r,s,p,q,p.p,q.q){
    return(S.value.complete(j,l,r,s,p,q,p.p,q.q) + S.value.incomplete(j,l,r,s,p,q,p.p,q.q))
  }

  #################Covariance Matrix of Z####################
  S.matrix <- function(j,l,r,s){
    return(matrix(c(sapply(1:trt, function(p){sapply(1:2, function(q){sapply(1:trt, function(p.p){sapply(1:2, function(q.q)
    {S.value(j,l,r,s,p,q,p.p,q.q)})})})})), nrow = trt*2))

    # matrix <- matrix(0,nrow = 2*trt, ncol = 2*trt)
    # for(a in 1:(2*trt)){
    #   for (b in 1:(2*trt)) {
    #     matrix[a,b] <- matrix[b,a] <- S.value(j,l,r,s,p,q,p.p,q.q)
    #   }
    # }
    # return(matrix)
  }

  S.column <- function(j,l){
    m.list <- lapply(1:trt, function(r){lapply(1:2, function(s){S.matrix(j,l,r,s)})})
    if(trt==1){
      a <- rbind(m.list[[1]][[1]], m.list[[1]][[2]])
    }
    if(trt>1){
      a <- rbind(m.list[[1]][[1]], m.list[[1]][[2]])
      for (i in 2:trt) {
        a <- rbind(a, m.list[[i]][[1]],m.list[[i]][[2]])
      }
    }
    return(a)
  }

  S.bycolumn <- lapply(1:trt, function(j){lapply(1:2, function(l){S.column(j,l)})})
  if(trt==1){
    Sigma <- cbind(S.bycolumn[[1]][[1]],S.bycolumn[[1]][[2]])
  }
  if(trt>1){
    Sigma <- cbind(S.bycolumn[[1]][[1]],S.bycolumn[[1]][[2]])
    for(i in 2:trt){
      Sigma <-  cbind(Sigma, S.bycolumn[[i]][[1]], S.bycolumn[[i]][[2]])
    }
  }

  #############check N here#############
  #Consistent Covariance Matrix Estimator
  V_N <- sum(N)*E%*%Sigma%*%t(E)
  #print(V_N)

  #Contrast for Different Tests (Default)
  T.trt <- kronecker(diag(trt)-1/trt*J(trt),1/2*J(2))
  T.time <- kronecker(1/trt*J(trt),diag(2)-1/2*J(2))
  T.inter <- kronecker(diag(trt)-1/trt*J(trt), diag(2)-1/2*J(2))

  T.list <- list(T.trt, T.time, T.inter)

  #Test Statistic Q_N(T)
  QNT.list <- sapply(1:3, function(i){sum(N)/(sum(diag(T.list[[i]]%*%V_N)))*t(p.vector)%*%T.list[[i]]%*%p.vector})

  ############Test Procedures############
  ##Box-type
  fhat.list <- lapply(1:3, function(i){sum(diag(T.list[[i]]%*%V_N))^2/sum(diag(T.list[[i]]%*%V_N%*%T.list[[i]]%*%V_N))})
  pvalue.box.list <- sapply(1:3, function(i){1-pchisq(fhat.list[[i]]*QNT.list[[i]],fhat.list[[i]])})

  z <- NULL
  #Formula
  z$formula <- formula
  #Number of Treatment
  z$trt <- trt
  #Sample Size - in the order of nc, n1, n2
  z$SampleSize <- lapply(1:trt, function(i){c(nc[i],n1[i],n2[i])})
  #Cluster Size - in the order of m1c, m2c, m1i, m2i (check)
  z$ClusterSize <- lapply(1:trt, function(i){list(mc.list[[1]][[i]],mc.list[[2]][[i]],mi.list[[1]][[i]],mi.list[[2]][[i]])})
  #complete and incomplete data by intervention
  z$x1c <- x1c
  z$x2c <- x2c
  z$x1i <- x1i
  z$x2i <- x2i
  #Total Sample Size - in the order of N1c, N2c, N1i, N2i (check)
  z$TotalSampleSize <- lapply(1:trt, function(i){c(Nc.list[[1]][[i]],Nc.list[[2]][[i]],Ni.list[[1]][[i]],Ni.list[[2]][[i]])})
  #Effect Size Vector - in the order of (p_11,p_12,p_21,p_22,...)
  z$p.vector <- p.vector
  #Contrast Matrix C and Contrast Matrix T
  if(is.null(Contrast)==FALSE) {
    #decide if Contrast is a vector or not
    if(is.vector(Contrast)==TRUE){
      ContrastT <- as.matrix(Contrast)%*%ginv(Contrast%*%as.matrix(Contrast))%*%Contrast
    }
    else{
      ContrastT <- t(Contrast)%*%ginv(Contrast%*%t(Contrast))%*%Contrast
    }
    ContrastTestStat <- round(sum(N)/(sum(diag(ContrastT%*%V_N)))*t(p.vector)%*%ContrastT%*%p.vector,4)
    fhat <- sum(diag(ContrastT%*%V_N))^2/sum(diag(ContrastT%*%V_N%*%ContrastT%*%V_N))
    pvalue.box <- round(1-pchisq(fhat*ContrastTestStat,fhat),4)
    #Input contrast matrix
    z$ContrastC <- Contrast
    #Projection of input contrast matrix
    z$ContrastT <- ContrastT
    #Test statistic of input contrast matrix
    z$ContrastTestStat <- ContrastTestStat
    #p-value of input contrast matrix
    z$pvalue.box <- pvalue.box
  }

  #Covariance Matrix
  z$CovMatrix <- round(V_N,4)
  #Test Statistic - Default (treatment,time,interaction)
  z$TestStat.default <- round(QNT.list,4)
  #p-value - Default (treatment,time,interaction)
  z$pvalue.default <- round(pvalue.box.list,4)
  z$indicator <- indicator
  class(z) <- "Main Results"
  z
}





