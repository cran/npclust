#' @title Box plots.
#' @description Box plot of the input data set by treatments and time period.
#' @name Plot.box
#' @param object a fitted model object from ncda() or a processed data set from ProcessData()
#' @return Box plots.
#' @usage Plot.box(object)
#' @examples
#' #Plot from analysis object
#' \donttest{
#' ARTIS_analysis <- ncda(emot_pqol~tx, ARTIS, intervention, homeid,
#'                        indicator = c("0","1"))
#' Plot.box(ARTIS_analysis)
#' }
#' skin_analysis <- ncda(score~tx, skin, intervention, subject,
#'                       indicator = c("control","treatment"))
#' Plot.box(skin_analysis)
#' # Plot from processed data set
#' ARTIS_result <- ProcessData(ARTIS, tx, intervention, homeid, symptoms_pqol,
#'                             indicator = c("0","1"))
#' skin_result <- ProcessData(skin, tx, intervention, subject, score,
#'                            indicator = c("control","treatment"))
#' Plot.box(ARTIS_result)
#' Plot.box(skin_result)
#' @importFrom graphics par segments
#' @import ggplot2
#' @export



##############
##object is the output from DataProcess
Plot.box <- function(object){
  #if(is.null(unlist(object$x1c))==TRUE){object$x1c <- NA}
  #if(is.null(unlist(object$x1i))==TRUE){object$x1i <- NA}
  #if(is.null(unlist(object$x2c))==TRUE){object$x2c <- NA}
  #if(is.null(unlist(object$x2i))==TRUE){object$x2i <- NA}
  predata <- lapply(1:object$trt, function(i){unlist(c(object$x1c[[i]],object$x1i[[i]]))})
  postdata <- lapply(1:object$trt, function(i){unlist(c(object$x2c[[i]],object$x2i[[i]]))})

  Response <- Intervention <- Trt <- NULL

  df.long.pre <- lapply(1:object$trt, function(i){data.frame(Response=as.matrix(predata[[i]]),
                    Trt=as.matrix(rep(paste("Treatment ",i),length(predata[[i]]))),
                    Intervention=as.matrix(rep(object$indicator[1],2*length(predata[[i]]))))})

  df.long.post <- lapply(1:object$trt, function(i){data.frame(Response=as.matrix(postdata[[i]]),
                     Trt=as.matrix(rep(paste("Treatment ",i),length(postdata[[i]]))),
                     Intervention=as.matrix(rep(object$indicator[2],2*length(postdata[[i]]))))})

  df.stack <- rbind(Reduce(f = "rbind", df.long.pre,), Reduce(f = "rbind", df.long.post,))
  df.stack$Intervention <- as.factor(df.stack$Intervention)
  df.stack$Intervention <- factor(df.stack$Intervention, levels = object$indicator)
  #Boxplot
  box <- ggplot(df.stack, aes(x = Trt, y = Response, color = Intervention)) +
    geom_boxplot()
  return(box)
}


