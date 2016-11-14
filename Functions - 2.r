
# psf(w)
psf <- function(w)pe*w^(-b)

# wf(ps)
wf <- function(ps)(ps/pe)^(-1/b)

# PLC(px)
PLCf <- function(px)1-exp(-(-px/d)^c)

# 50 % PLC loss
# c
Psi50fc <- function(c){
  f1 <- function(px)exp(-(-px/d)^c)-0.5
  res <- uniroot(f1, c(-100, 0), tol=.Machine$double.eps)$root
  return(res)
}

# d
Psi50fd <- function(d){
  f1 <- function(px)exp(-(-px/d)^c)-0.5
  res <- uniroot(f1, c(-100, 0), tol=.Machine$double.eps)$root
  return(res)
}

# Inverse d
InvPsi50fd <- function(px){
  f1 <- function(d)exp(-(-px/d)^c)-0.5
  res <- uniroot(f1, c(0, 100), tol=.Machine$double.eps)$root
  return(res)
}

# xylem conductance function
kxf <- function(px)kxmax*exp(-(-px/d)^c)

# minimum xylem water potential function at given w
pxminf <- function(w){
  ps <- psf(w)
  f1 <- function(px)(ps-px)*h2*kxf(px)
  res <- optimize(f1, c(-20, 0), tol=.Machine$double.eps, maximum=T)$maximum
  return(res)
}

# gsmaxf(w)
gsmaxf <- Vectorize(function(w){
  ps <- psf(w)
  pxmin <- pxminf(w)
  res <- (ps-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
})

# xylem water potential function
pxf <- function(w, gs){
  ps <- psf(w)
  pxmin <- pxminf(w)
  f1 <- function(px)((ps-px)*h2*kxf(px)-h*VPD*gs)^2
  res <- ifelse(pxmin<ps, optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum, ps)
  return(res)
}

# Af(gs)
Af <- function(gs)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# dAdgs
dAdgsf <- function(gs)(1/2)*LAI*(ca+Km+((-ca^2)*gs-gs*Km^2-Km*Rd-2*cp*Vcmax-Km*Vcmax+ca*(-2*gs*Km-Rd+Vcmax))/sqrt((ca*gs-gs*Km+Rd-Vcmax)^2+4*gs*(ca*gs*Km+Km*Rd+cp*Vcmax)))

# PLCwgsf(w, gs)
PLCwgsf <- function(w, gs){
  px <- pxf(w, gs)
  res <- PLCf(px)
  return(res)
}

# mf(w, gs)
mf <- function(w, gs)h3*PLCwgsf(w, gs)

# B(w, gs)
Bf <- function(w, gs)Af(gs)-mf(w, gs)

# ESS gs(w)
ESSf <- Vectorize(function(w){
  f1 <- function(gs)Bf(w, gs)
  res <- ifelse(gsmaxf(w)>0, optimize(f1, c(0, gsmaxf(w)), tol=.Machine$double.eps, maximum=T)$maximum, 0)
  return(res)
})

# ESS A(w)
ESSAf <- function(w)Af(ESSf(w))

# ESS m(w)
ESSmf <- Vectorize(function(w)mf(w, ESSf(w)))

# ESS B(w)
ESSBf <- function(w)ESSAf(w)-ESSmf(w)

# ESS Ev(w)
ESSEvf <- function(w)h*VPD*ESSf(w)

# ESS Ev(ps)
ESSEvpsf <- function(ps)h*VPD*ESSf(wf(ps))

# ESS g1(ps)
ESSg1psf <- Vectorize(function(ps){
  f1 <- function(w)psf(w)-ps
  w <- uniroot(f1, c(0.001, 1), tol=.Machine$double.eps)$root
  res <- sqrt(VPD*100)*(ca*ESSf(w)/(a*ESSAf(w))-1)
  return(res)
})

# ESS gs(ps)
ESSpsf <- function(ps){
  w <- wf(ps)
  res <- ESSf(w)
  return(res)
}

# ESS px(w)
ESSpxpsf <- Vectorize(function(ps){
  w <- wf(ps)
  res <- pxf(w, ESSf(w))
  return(res)
})

# ESS PLC(w)
ESSPLCf <- function(w){
  px <- pxf(w, ESSf(w))
  res <- PLCf(px)
  return(res)
}

# ESS PLC(ps)
ESSPLCpsf <- function(ps){
  w <- wf(ps)
  px <- pxf(w, ESSf(w))
  res <- PLCf(px)
  return(res)
}

# ESS cica(w)
ESScicaf <- function(w)1-ESSAf(w)/ESSf(w)/ca
