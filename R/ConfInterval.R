#' @title Confidence Interval
#'
#' @description  Construct confidence intervals for effect sizes.
#'
#' @param object a fitted model object from ncda().
#' @param level the confidence level required.
#' @param side a character string specifying the side of the confidence bound,
#' must be one of "two.sided" (default), "left" or "right".
#' @param adjust an optional character string specifying the multiple adjustment method,
#' by default there is no adjustment, if specified, must be one of "Bonferroni" or "Working-Hotelling".
#' You can specify just the initial letter.
#'
#' @usage ConfInterval(object, level, side="two.sided",
#' adjust=NULL)
#'
#' @return A list or a vector. If the confidence interval is two-sided, lower and
#' upper bounds are stored in lists for each nonparametric effect size estimate.
#' Otherwise, the lower/upper bounds are stored in vectors in the order of the effect size estimates.
#'
#' @examples
#' skin_analysis <- ncda(score~tx, skin, intervention, subject,
#'                       indicator=c("control","treatment"),
#'                       Contrast=matrix(c(1,-1), nrow = 1))
#' ConfInterval(skin_analysis,0.95)
#' \donttest{
#' ARTIS_analysis <- ncda(emot_pqol~tx, ARTIS, intervention, homeid,
#'                        indicator = c("0","1"))
#' ConfInterval(ARTIS_analysis,0.95)
#' ConfInterval(ARTIS_analysis,0.95,"two.sided","B")
#' ConfInterval(ARTIS_analysis,0.95,"left","W")
#' }
#' @importFrom stats pchisq qf qnorm
#' @export

ConfInterval <- function(object,level,side="two.sided",adjust=NULL){
  vars <- diag(object$CovMatrix)
  N <- sum(sapply(1:object$trt, function(i){sapply(object$ClusterSize[[i]],sum)}))
  #inverse function - logit(x)=log(x/(1-x))
  #default - two-sided CI
  if(is.null(side)==TRUE){side="two.sided"}
  if(side=="two.sided"){
    alpha <- 1/2*(1-level)
    #check if the CI's need to be adjusted
    #default - no adjustment
    if(is.null(adjust)==TRUE){
      lower.logit.origscale <- log(object$p.vector/(1-object$p.vector))-qnorm(1-alpha/2)/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      upper.logit.origscale <- log(object$p.vector/(1-object$p.vector))+qnorm(1-alpha/2)/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      lower.logit <- exp(lower.logit.origscale)/(1+exp(lower.logit.origscale))
      upper.logit <- exp(upper.logit.origscale)/(1+exp(upper.logit.origscale))
      #return(sapply(1:length(object$p.vector), function(i){paste0("(",round(lower.logit[i],3),",",round(upper.logit[i],3),")")}))
      #return(sapply(1:length(object$p.vector), function(i){paste0("(",round(lower.logit[i],3),",",round(upper.logit[i],3),")")}))
    }
    #Bonferroni adjustment
    else if(adjust=="Bonferroni" | adjust=="B"){
      lower.logit.origscale <- log(object$p.vector/(1-object$p.vector))-qnorm(1-alpha/(2*2*object$trt))/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      upper.logit.origscale <- log(object$p.vector/(1-object$p.vector))+qnorm(1-alpha/(2*2*object$trt))/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      lower.logit <- exp(lower.logit.origscale)/(1+exp(lower.logit.origscale))
      upper.logit <- exp(upper.logit.origscale)/(1+exp(upper.logit.origscale))
      #return(sapply(1:length(object$p.vector), function(i){paste0("(",round(lower.logit[i],3),",",round(upper.logit[i],3),")")}))
    }
    #Working-Hotelling adjustment
    else if(adjust=="Working-Hotelling" | adjust=="W"){
      lower.logit.origscale <- log(object$p.vector/(1-object$p.vector))-sqrt(2*qf(1-alpha,df1 = 2*object$trt,df2 = N-2*object$trt))/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      upper.logit.origscale <- log(object$p.vector/(1-object$p.vector))+sqrt(2*qf(1-alpha,df1 = 2*object$trt,df2 = N-2*object$trt))/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      lower.logit <- exp(lower.logit.origscale)/(1+exp(lower.logit.origscale))
      upper.logit <- exp(upper.logit.origscale)/(1+exp(upper.logit.origscale))
      #return(sapply(1:length(object$p.vector), function(i){paste0("(",round(lower.logit[i],3),",",round(upper.logit[i],3),")")}))
    }
    return(lapply(1:length(object$p.vector), function(i){c(round(lower.logit[i],3),round(upper.logit[i],3))}))
  }
  else if(side=="left"){
    alpha <- 1-level
    if(is.null(adjust)==TRUE){
      lower.logit.origscale <- log(object$p.vector/(1-object$p.vector))-qnorm(alpha)/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      lower.logit <- exp(lower.logit.origscale)/(1+exp(lower.logit.origscale))
      return(round(lower.logit,3))
    }
    else if(adjust=="Bonferroni" | adjust=="B"){
      lower.logit.origscale <- log(object$p.vector/(1-object$p.vector))-qnorm(1-alpha/(2*object$trt))/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      lower.logit <- exp(lower.logit.origscale)/(1+exp(lower.logit.origscale))
      return(round(lower.logit,3))
    }
    else if(adjust=="Working-Hotelling" | adjust=="W"){
      lower.logit.origscale <- log(object$p.vector/(1-object$p.vector))-sqrt(2*qf(1-alpha,df1 = 2*object$trt,df2 = N-2*object$trt))/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      lower.logit <- exp(lower.logit.origscale)/(1+exp(lower.logit.origscale))
      return(round(lower.logit,3))
    }
  }
  else if(side=="right"){
    alpha <- 1-level
    if(is.null(adjust)==TRUE){
      upper.logit.origscale <- log(object$p.vector/(1-object$p.vector))+qnorm(alpha)/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      upper.logit <- exp(upper.logit.origscale)/(1+exp(upper.logit.origscale))
      return(round(upper.logit,3))
    }
    else if(adjust=="Bonferroni" | adjust=="B"){
      upper.logit.origscale <- log(object$p.vector/(1-object$p.vector))+qnorm(1-alpha/(2*object$trt))/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      upper.logit <- exp(upper.logit.origscale)/(1+exp(upper.logit.origscale))
      return(round(upper.logit,3))
    }
    else if(adjust=="Working-Hotelling" | adjust=="W"){
      upper.logit.origscale <- log(object$p.vector/(1-object$p.vector))+sqrt(2*qf(1-alpha,df1 = 2*object$trt,df2 = N-2*object$trt))/sqrt(N)*sqrt(vars)*1/(object$p.vector*(1-object$p.vector))
      upper.logit <- exp(upper.logit.origscale)/(1+exp(upper.logit.origscale))
      return(round(upper.logit,3))
    }
  }
}
