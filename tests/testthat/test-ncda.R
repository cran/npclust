require(npclust)
ncda(score~tx,skin,intervention,subject,Contrast=matrix(c(1,-1), nrow = 1))
