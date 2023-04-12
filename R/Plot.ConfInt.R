#' @title Bar plots for two-sided confidence intervals
#' @name Plot.ConfInt
#' @param object a fitted model object from ncda().
#' @param level the confidence level required.
#' @param side a character string specifying the side of the confidence bound,
#' must be one of "two.sided" (default), "left" or "right".
#' @param adjust an optional character string specifying the multiple adjustment method,
#' by default there is no adjustment, if specified, must be one of "Bonferroni" or
#' "Working-Hotelling". You can specify just the initial letter.
#' @usage Plot.ConfInt(object, level, side="two.sided",
#' adjust=NULL)
#' @return Bar plots.
#' @examples
#' skin_analysis <- ncda(score~tx, skin, intervention, subject,
#'                       indicator=c("control","treatment"),
#'                       Contrast=matrix(c(1,-1), nrow = 1))
#' Plot.ConfInt(skin_analysis,0.95,"Two-Sided")
#' \donttest{
#' ARTIS_analysis <- ncda(emot_pqol~tx, ARTIS, intervention, homeid,
#'                        indicator = c("0","1"))
#' Plot.ConfInt(ARTIS_analysis,0.95,"Two-Sided")
#' Plot.ConfInt(ARTIS_analysis,0.95,"Two-Sided","Bonferroni")
#' Plot.ConfInt(ARTIS_analysis,0.95,"Two-Sided","Working-Hotelling")
#' }
#' @export

Plot.ConfInt <- function(object,level,side="two.sided",adjust=NULL){
  if(side=="left" | side=="right"){
    stop("Error Bar is only available for two-sided intervals.")
  }
  else{
    intervals <- ConfInterval(object,level,"two.sided",adjust)
    lower <- sapply(1:length(intervals), function(i){intervals[[i]][1]})
    upper <- sapply(1:length(intervals), function(i){intervals[[i]][2]})

    #par(mar = c(5, 4, 4, 2) + 0.1)

    plot(seq(1:(2*object$trt)),object$p.vector, ylim = c(0,1),
         xlim = c(0.5,(2*object$trt+0.5)), pch = 16,
         ylab = "Effect Size", xlab = "Arm",
         main = "Effect Size Estimates with Confidence Intervals")

    add.error.bars <- function(X,upper,lower,width,col=par( )$fg,lwd=1){
      segments(X,lower,X,upper,col=col,lwd=lwd,lend=1);
      segments(X-width/2,lower,X+width/2,lower,col=col,lwd=lwd,lend=1);
      segments(X-width/2,upper,X+width/2,upper,col=col,lwd=lwd,lend=1);
    }

    add.error.bars(seq(1:(2*object$trt)),lower,upper,width = 0.5)
  }
}

