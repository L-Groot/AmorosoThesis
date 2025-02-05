#############################################
### The Amoroso PDF defined from scratch! ###
#############################################

# --> just for fun :)
# --> equivalent to dgg4()! from the 'AmoRosoDistrib' package

dAmoroso <- function(x, a, l, c, mu) {
  c1 <- 1/(gamma(l))
  c2 <- abs(c/a)
  c3 <- ((x-mu)/a)^(l*c-1)
  c4 <- exp(-((x-mu)/a)^c)
  return(c1*c2*c3*c4)
}