Fhat <- function(x,x11){
  sum((x>x11)+1/2*(x==x11))/length(x11)
}
