#' @title Process data set.
#' @description Sample size and cluster size calculation for the imported data set.
#' @param data a data frame in the long format.
#' @param tx treatment variable.
#' @param period time indicator variable.
#' @param subject subject or cluster ID
#' @param resp response variable to be analyzed.
#' @param indicator an optional vector of characters indicating the order of pre and post intervention period;
#' must match the levels of period argument if specified; if not specified, the pre and post intervention
#' period will be ordered in the alphabet order by default
#' @return a list containing the following components:
#' \describe{
#'   \item{\code{trt }}{number of treatments}
#'   \item{\code{nc }}{complete cluster sample size within each treatment group}
#'   \item{\code{n1 }}{incomplete cluster sample size pre intervention within each treatment group}
#'   \item{\code{n2 }}{incomplete cluster sample size post intervention within each treatment group}
#'   \item{\code{m1c }}{complete cluster sizes pre-intervention within each treatment group}
#'   \item{\code{m2c }}{complete cluster size post-intervention within each treatment group}
#'   \item{\code{m1i }}{incomplete cluster sizes pre-intervention within each treatment group}
#'   \item{\code{m2i }}{incomplete cluster sizes post-intervention within each treatment group}
#'   \item{\code{x1c }}{complete data pre-intervention within each treatment group}
#'   \item{\code{x2c }}{complete data post-intervention within each treatment group}
#'   \item{\code{x1i }}{incomplete data pre-intervention within each treatment group}
#'   \item{\code{x2i }}{incomplete data post-intervention within each treatment group}
#' }
#'
#' @usage ProcessData(data, tx, period, subject, resp, indicator=NULL)
#'
#' @examples
#' ARTIS_result <- ProcessData(ARTIS, tx, intervention, homeid, symptoms_pqol,
#'                             c("0","1"))
#' names(ARTIS_result)
#' skin_result <- ProcessData(skin, tx, intervention, subject, score,
#'                            c("control","treatment"))
#' skin_result$nc
#' skin_result$n1
#' skin_result$n2
#' @importFrom utils head tail
#' @export

ProcessData <- function(data,tx,period,subject,resp,indicator=NULL){
  #check if input data is data frame
  if(is.data.frame(data)==FALSE) stop("data must be a data frame")
  x <- data.frame(tx=eval(substitute(tx),data),intervention=eval(substitute(period),data),
                  clusterid=eval(substitute(subject),data),resp=eval(substitute(resp),data))
  if(length(unique(x$intervention))!=2) stop("only comparison between two intervention periods is allowed")
  x$intervention <- as.character(x$intervention)
  #count number of treatments
  trt <- sum(table(x$tx)!=0)
  #delete obs with NA
  df.trt <- lapply(1:trt, function(i){x[which((x$tx==levels(factor(x$tx))[i])&(is.na(x$resp)==FALSE)),]})
  #count complete and incomplete clusters
  pattern <- lapply(1:trt, function(i){table(df.trt[[i]]$clusterid,df.trt[[i]]$intervention)})
  #rearrange dataset by the input pre/post indicator
  if(is.null(indicator)==FALSE){
    if(identical(sort(unique(x$intervention)),sort(indicator))==FALSE) stop("indicator argument must match levels of period argument in the data set")
    else pattern <- lapply(1:trt, function(i){pattern[[i]][,indicator]})
  }
  #identify complete/incomplete clusters
  cluster.complete <- lapply(1:trt, function(j){as.numeric(rownames(pattern[[j]])[sapply(1:dim(pattern[[j]])[1], function(i){(pattern[[j]][i,1]!=0)&(pattern[[j]][i,2]!=0)})])})
  cluster.incomplete.pre <- lapply(1:trt, function(j){as.numeric(rownames(pattern[[j]])[sapply(1:dim(pattern[[j]])[1], function(i){(pattern[[j]][i,1]!=0)&(pattern[[j]][i,2]==0)})])})
  cluster.incomplete.post <- lapply(1:trt, function(j){as.numeric(rownames(pattern[[j]])[sapply(1:dim(pattern[[j]])[1], function(i){(pattern[[j]][i,1]==0)&(pattern[[j]][i,2]!=0)})])})

  #calculate sample size
  nc <- sapply(1:trt, function(i){length(cluster.complete[[i]])})
  n1 <- sapply(1:trt, function(i){length(cluster.incomplete.pre[[i]])})
  n2 <- sapply(1:trt, function(i){length(cluster.incomplete.post[[i]])})

  #get complete and incomplete data (check zero sample size)
  data.complete <- lapply(1:trt, function(i){if(nc[i]==0){NULL} else{df.trt[[i]][which(df.trt[[i]]$clusterid %in% cluster.complete[[i]]),]}})
  data.incomplete.pre <- lapply(1:trt, function(i){if(n1[i]==0){NULL} else {df.trt[[i]][which(df.trt[[i]]$clusterid %in% cluster.incomplete.pre[[i]]),]}})
  data.incomplete.post <- lapply(1:trt, function(i){if(n2[i]==0) {NULL} else {df.trt[[i]][which(df.trt[[i]]$clusterid %in% cluster.incomplete.post[[i]]),]}})

  #generate xc x1i x2i
  if(trt==1){
    xc <- if(nc[1]!=0){list(lapply(1:nc[[1]], function(j){data.complete[[1]][which(data.complete[[1]]$clusterid==cluster.complete[[1]][j]),"resp"]}))} else{NULL}
    x1i <- if(n1[1]!=0){list(lapply(1:n1[[1]], function(j){data.incomplete.pre[[1]][which(data.incomplete.pre[[1]]$clusterid==cluster.incomplete.pre[[1]][j]),"resp"]}))} else{NULL}
    x2i <- if(n2[1]!=0){list(lapply(1:n2[[1]], function(j){data.incomplete.post[[1]][which(data.incomplete.post[[1]]$clusterid==cluster.incomplete.post[[1]][j]),"resp"]}))} else{NULL}
  }else{
    xc <- sapply(1:trt, function(i){if(nc[i]!=0){lapply(1:nc[[i]], function(j){data.complete[[i]][which(data.complete[[i]]$clusterid==cluster.complete[[i]][j]),"resp"]})} else{NULL}})
    x1i <- sapply(1:trt, function(i){if(n1[i]!=0){lapply(1:n1[[i]], function(j){data.incomplete.pre[[i]][which(data.incomplete.pre[[i]]$clusterid==cluster.incomplete.pre[[i]][j]),"resp"]})} else{NULL}})
    x2i <- sapply(1:trt, function(i){if(n2[i]!=0){lapply(1:n2[[i]], function(j){data.incomplete.post[[i]][which(data.incomplete.post[[i]]$clusterid==cluster.incomplete.post[[i]][j]),"resp"]})} else{NULL}})
  }

  #calculate cluster size
  m1c <- lapply(1:trt, function(i){pattern[[i]][which(as.numeric(rownames(pattern[[i]])) %in% cluster.complete[[i]]),1]})
  m2c <- lapply(1:trt, function(i){pattern[[i]][which(as.numeric(rownames(pattern[[i]])) %in% cluster.complete[[i]]),2]})
  m1i <- lapply(1:trt, function(i){pattern[[i]][which(as.numeric(rownames(pattern[[i]])) %in% cluster.incomplete.pre[[i]]),1]})
  m2i <- lapply(1:trt, function(i){pattern[[i]][which(as.numeric(rownames(pattern[[i]])) %in% cluster.incomplete.post[[i]]),2]})

  #separate complete data into lists
  x1c <- lapply(1:trt, function(i){if(nc[i]==0){NULL} else{lapply(1:nc[i],function(arg){head(unlist(xc[[i]][arg]),m1c[[i]][arg])})}})
  x2c <- lapply(1:trt, function(i){if(nc[i]==0){NULL} else{lapply(1:nc[i],function(arg){tail(unlist(xc[[i]][arg]),-m1c[[i]][arg])})}})

  #output
  z <- NULL
  z$trt <- trt
  z$nc <- nc
  z$n1 <- n1
  z$n2 <- n2
  z$m1c <- m1c
  z$m2c <- m2c
  z$m1i <- m1i
  z$m2i <- m2i
  z$x1c <- x1c
  z$x2c <- x2c
  z$x1i <- x1i
  z$x2i <- x2i
  z$indicator <- indicator
  class(z) <- "Transformed Data"
  z
}








