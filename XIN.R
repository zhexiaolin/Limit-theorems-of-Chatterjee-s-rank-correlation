XIncalculate = function(xvec, yvec){
  n = length(xvec)
  xrank = rank(xvec, ties.method = "random")
  yrank = rank(yvec, ties.method = "random")
  ord = order(xrank)
  yrank = yrank[ord]
  coef.sum = sum(abs(yrank[2:n] - yrank[1:n-1]))
  coef.value = 1-3*coef.sum/(n^2-1)
  return(coef.value)
}

XInvarcalculate = function(xvec, yvec){
  n = length(xvec)
  xrank = rank(xvec, ties.method = "random")
  yrank = rank(yvec, ties.method = "random")
  ord = order(xrank)
  yrank = yrank[ord]
  yrank1 = yrank[c(2:n,n)]
  yrank2 = yrank[c(3:n,n,n)]
  yrank3 = yrank[c(4:n,n,n,n)]
  term1 = pmin(yrank,yrank1)
  term2 = pmin(yrank,yrank2)
  term3 = pmin(yrank2,yrank3)
  term4 = sapply(1:n, function(i){sum(yrank[i]<=(term1[-i]))})
  term5 = pmin(yrank1,yrank2)
  sum1 = mean((term1/n)^2)
  sum2 = mean(term1*term2/n^2)
  sum3 = mean(term1*term3/n^2)
  sum4 = mean(term4*term1/(n*(n-1)))
  sum5 = mean(term4*term5/(n*(n-1)))
  sum6 = mean(sapply(1:n, function(i){sum(pmin(term1[i],term1[-i]))})/(n*(n-1)))
  sum7 = (mean(term1/n))^2
  return(max(0,36*(sum1+2*sum2-2*sum3+4*sum4-2*sum5+sum6-4*sum7)))
}


