decay_titer <- function (
    t,GMT0,k1,k2,f
  ){
  GMT0+log10(f*exp(-k1*t)+(1-f)*exp(-k2*t))
}