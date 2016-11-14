
source("Functions - 2.r")
data <- read.csv("data/Choat2012.csv")
dataAng <- subset(data, Type=="Angiosperm", select=c("Psi50", "Psimin"))
dataAng <- dataAng[order(dataAng$Psi50), ]
dataGym <- subset(data, Type=="Gymnosperm", select=c("Psi50", "Psimin"))
dataGym <- dataGym[order(dataGym$Psi50), ]

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
gamma <- 1/((MAP/365/k)/1000)*nZ

# Regression
fitAng <- nls(Psimin ~ -a*Psi50+b, data=dataAng, start=list(a=-0.4, b=-1),
              control=c(minFactor=1e-5))
fitGym <- nls(Psimin ~ -a*Psi50+b, data=dataGym, start=list(a=-0.4, b=-1),
              control=c(minFactor=1e-5))

# Sensitivity Analysis
SA1 <- seq(1/4.5, 4.5, by=0.1)*d
SA2 <- c(1, 25, 100)
res <- data.frame(PLC50=numeric(length(SA1)), PLCmin=numeric(length(SA1)),
                  PLC50=numeric(length(SA1)), PLCmin=numeric(length(SA1)),
                  PLC50=numeric(length(SA1)), PLCmin=numeric(length(SA1)))

for(i in 1:length(SA2)){
  h3 <- SA2[i]
  for(j in 1:length(SA1)){
    d <- SA1[j]
    res[j, 2*i-1] <- Psi50fd(d)
    wL <- uniroot(ESSBf, c(0.05, 1), tol=.Machine$double.eps)$root
    res[j, 2*i] <- pxf(wL, ESSf(wL))
  }
}

# Figures
Cols <- c("lightblue", "lightpink", "purple", "forestgreen", "orange")
windows(18, 12)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", mar=c(8, 8, 1, 1), mfrow=c(1, 1))
plot(0, 0, type="n",
     xaxt="n", yaxt="n", xlab=NA, ylab=NA,
     xlim=c(-15, 0), ylim=c(-15, 0), lwd=8)

axis(1, xlim=c(-15, 0), pos=-15, at=c(-15, -10, -5, 0), cex.axis=2, lwd=4)
mtext(expression(P[50]~(MPa)),side=1,line=6.5, cex=5)
axis(2, ylim=c(-15, 0), pos=-15, at=c(-15, -10, -5, 0), cex.axis=2, lwd=4)
mtext(expression(P[min]~(MPa)),side=2,line=3, cex=5)
abline(a=0, b=1, lwd=3, lty=3)

points(dataAng, type="p", col=Cols[1], pch=1, cex=3, lwd=3)
lines(dataAng$Psi50, predict(fitAng), col=Cols[1], lwd=8, lty=2)
points(dataGym, type="p", col=Cols[2], pch=2, cex=3, lwd=3)
lines(dataGym$Psi50, predict(fitGym), col=Cols[2], lwd=8, lty=2)

points(res[1:2], type="l", col=Cols[3], lwd=8)
points(res[3:4], type="l", col=Cols[4], lwd=8)
points(res[5:6], type="l", col=Cols[5], lwd=8)

legend("bottomright", legend=SA2, title=expression(beta), lty=c(1), col=Cols[3:5], cex=3, lwd=8, box.lwd=8)
legend("bottomleft", title=expression(Choat~italic(et~al.)~2012), c("Angiosperm", "Gymnosperm"), pch=15, col=Cols[1:2], cex=3, box.lwd=8, bg="white")
box(lwd=8)

dev.copy2pdf(file = "Figures/Choat et al 2012.pdf")
