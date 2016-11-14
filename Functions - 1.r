
# PLC(px)
PLCf <- function(px, kxmax=5, c=2.64, d=3.54){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  kx <- kxf(px)
  res <- 1-kx/kxmax
  return(res)
}

# ps(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# pxmin(w)
pxminf <- function(w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54){
  
  f1 <- function(x)-((psf(w)-x)*h2*kxmax*exp(-(-x/d)^c))
  res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
  return(res)
}

# intact gsmax
gsmaxf <- function(w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    f1 <- function(x)(psf(w)-x)*h2*kxf(x)
    res <- optimize(f1, c(-20, pe), tol=.Machine$double.eps, maximum=T)$maximum
    return(res)
  }
  
  pxmin <- pxminf(w)
  res <- (psf(w)-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
}

# family ESS
gswLf <- function(w, wL,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54){
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  pxmin <- psf(wL)
  res <- (psf(w)-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
}

# family ESS as px(w)
pxwLf <- function(w, wL,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                  h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  pxmin <- psf(wL)
  res <- psf(w)-h*VPD*gswLf(w, wL)/kxf(pxmin)/h2
  return(res)
}

# family ESS as px(ps)
pxpswLf <- function(ps, wL, pe=-1.58*10^-3, b=4.38){
  
  w <- (ps/pe)^(-1/b)
  res <- pxwLf(w, wL)
  return(res)
}

# A
Af <- function(gs, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1){LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))}

# averBi
averBif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  gswLfr <- function(w)ifelse(w<wLr+1e-10, 0, gswLf(w, wLr))
  gswLfi <- function(w)ifelse(w<wLi+1e-10, 0, gswLf(w, wLi))
  
  Ef <- function(w){h*VPD*gswLfr(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.5)$value})
  
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wLr-1e-10)/(1-wLr-1e-10)-k*integralrEf(w)*1/(1-wLr-1e-10))*1/(1-wLr-1e-10)
  fLnoc <- 1/k*exp(-k*integrate(rEf, wLr+1e-10, 1, rel.tol=.Machine$double.eps^0.5)$value*1/(1-wLr-1e-10))
  cPDF <- 1/(fLnoc+integrate(fnoc, wLr+1e-10, 1, rel.tol=.Machine$double.eps^0.4)$value)
  
  f1 <- function(w)cPDF*fnoc(w)
  f2 <- function(w)Af(gswLfi(w))*f1(w)
  fL <- cPDF*fLnoc
  res1 <- integrate(f2, wLr+1e-10, 1, rel.tol=.Machine$double.eps^0.3)
  res2 <- fL*Af(gswLfi(wLr+1e-10))
  res <- res1$value+res2
  return(res)
}
