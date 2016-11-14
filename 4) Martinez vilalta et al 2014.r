
source("Functions - 2.r")
data <- read.csv("Data/Martinez-vilalta 2014.csv")

# Parameterization
ca <- 400
k <- 0.05
MAP <- 1000
LAI <- 1
Vcmax <- 50
cp <- 30
Km <- 703
Rd <- 1
a <- 1.6
nZ <- 0.5
p <- 43200
l <- 1.8e-5
VPD <- 0.02
pe <- -1.58*10^-3
b <- 4.38
kxmax <- 5
c <- 2.64
d <- 3.54
h <- l*a*LAI/nZ*p
h2 <- l*LAI/nZ*p/1000
h3 <- 10
gamma <- 1/((MAP/365/k)/1000)*nZ

# Martinez-vilalta 2014
rawdata <- read.csv("Data/Martinez-vilalta 2014.csv")
data <- data.frame(P50=rawdata$P50.Mpa, sigma=rawdata$sigma)
data <- data[order(data$P50), ]

fit <- nls(sigma ~ a*(-P50)^b+c, data=data, start=list(a=0.8459, b=0.1261, c=-0.1320),
           control=c(minFactor=1e-5))

# Sensitivity analysis
SA <- seq(1, 10, by=0.5)
Psi50 <- sapply(SA, Psi50fd)
Slope <- numeric(length=length(SA))
for(i in 1:length(SA)){
  d <- SA[i]
  wL <- uniroot(ESSBf, c(0.1, 1), tol=.Machine$double.eps)$root
  psL <- psf(wL)
  Slope[i] <- (ESSpxpsf(psL)-ESSpxpsf(pe))/(psL-pe)
}

# Figure
windows(18, 12)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", mar=c(8, 8, 1.2, 1), mfrow=c(1, 1))
plot(0, 0,
     type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA,
     xlim=c(-10, 0), ylim=c(0.2, 1.4), lwd=8)

points(data$P50, data$sigma, type="p", lwd=3, cex=3)
lines(data$P50, predict(fit), lty=3, lwd=8)
points(Psi50, Slope, type="l", cex=3, lwd=8)

axis(1, xlim=c(-10, 0), pos=0.2, cex.axis=2, lwd=4)
mtext(expression(psi[x50]~(MPa)),side=1,line=6.5, cex=5)
axis(2, ylim=c(0.2, 1.4), pos=-10, cex.axis=2, lwd=4)
mtext(expression(Slope~of~psi[x]*(psi[s])), side=2, line=3, cex=5)

legend("bottomleft",
       legend=expression(Martinez-vilalta~italic(et~al.)~2014, "Our model prediction"),
       lty=c(3, 1), pch=c(1, NA), cex=3, lwd=6, box.lwd=8)
box(lwd=8)

dev.copy2pdf(file = "Figures/Martinez vilalta et al 2014.pdf")
