
source("Functions - 2.r")

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
#d <- 3.54
h <- l*a*LAI/nZ*p
h2 <- l*LAI/nZ*p/1000
gamma <- 1/((MAP/365/k)/1000)*nZ

# Sensitivity Analysis
SA1 <- seq(0.5, 10, by=0.5)
SA2 <- c(1, 25, 100)
data <- data.frame(PLC50=numeric(length(SA1)), gs50=numeric(length(SA1)),
                   PLC50=numeric(length(SA1)), gs50=numeric(length(SA1)),
                   PLC50=numeric(length(SA1)), gs50=numeric(length(SA1)))

for(i in 1:length(SA2)){
  h3 <- SA2[i]
  for(j in 1:length(SA1)){
    d <- SA1[j]
    g1 <- ESSf(1)
    f1 <- function(w)ESSf(w)-g1*0.5
    wL <- uniroot(ESSBf, c(0.12, 1), tol=.Machine$double.eps)$root
    w50 <- uniroot(f1, c(wL, 1), tol=.Machine$double.eps)$root
    data[j, 2*i-1] <- Psi50fd(d)
    data[j, 2*i] <- pxf(w50, ESSf(w50))
  }
}

# Figures
Cols <- c("black", "purple", "forestgreen", "orange")
windows(18, 12)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", mar=c(8, 8, 1, 1), mfrow=c(1, 1))
plot(0, 0, type="n",
     xaxt="n", yaxt="n", xlab=NA, ylab=NA,
     xlim=c(-10, 0), ylim=c(-10, 0), cex.lab=1.3, lwd=8)

axis(1, xlim=c(-10, 0), pos=-10, cex.axis=2, lwd=4)
mtext(expression(P[50]~(MPa)), side=1, line=6.5, cex=5)
axis(2, ylim=c(-10, 0), pos=-10, cex.axis=2, lwd=4)
mtext(expression(P[50*", "*italic(g[s])]~(MPa)), side=2, line=2.8, cex=5)
abline(a=0, b=1, lwd=3, lty=3)

points(data[1:2], type="l", col=Cols[2], lwd=8, lty=2)
points(data[3:4], type="l", col=Cols[3], lwd=8, lty=2)
points(data[5:6], type="l", col=Cols[4], lwd=8, lty=2)

curve(0.49*x-0.42, -7, -1, lty=1, add=T, lwd=8)

legend("bottomright", legend=SA2, title=expression(With~beta*"="), lty=2, col=Cols[2:4], cex=3, lwd=8, box.lwd=8)
legend("bottomleft", c("Klein 2014"), lty=1, col=Cols[1], cex=3, lwd=8, box.lwd=8, bg="white")
box(lwd=8)

dev.copy2pdf(file = "Figures/Klein 2014.pdf")
