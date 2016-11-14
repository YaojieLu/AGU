
source("Functions - 1.r")

# Initialize
pe <- -1.58*10^-3
psL <- -3.5
pkx <- 0.5

PLCf1 <- function(px)PLCf(px)*100
PLCmax <- PLCf1(psL)
PLCmf <- function(px)PLCmax-(PLCmax-PLCf1(px))*pkx

# Figure
windows(18, 12)
par(mgp=c(2.2, 1, 0), xaxs="i", lwd=8, mar=c(8, 8, 1, 1), mfrow=c(1, 1))
curve(PLCf1, -6, pe,
      xlim=c(-6, 0), ylim=c(0, 100),
      xlab=NA, ylab=NA, xaxt="n", yaxt="n")

#curve(PLCf1, psL, pe, col="red", add=T, lty=2)
#segments(psL, PLCmax, pe, PLCmax, col="orange", lty=2)
curve(PLCmf, psL, pe, col="blue", add=T, lty=2)

axis(1, xlim=c(-6, 0), pos=-100*0.04, cex.axis=2, lwd=4)
mtext(expression(psi[x]~(MPa)),side=1,line=6.5, cex=5)
axis(2, ylim=c(0, 100), pos=-6, cex.axis=2, lwd=4)
mtext("PLC (%)",side=2,line=4, cex=5)

arrows(-1, PLCmf(-1)*1.1, -1, PLCmax, lwd=4, col="blue")
arrows(-1, PLCmf(-1)*0.9, -1, PLCf1(-1), lwd=4, col="blue")
arrows(psL, PLCmax, psL, -4, lty=2, lwd=4)
text(psL-0.59, PLCmax/2, expression(psi[xmin]), cex=5)
arrows(0, PLCmax, -6, PLCmax, lty=2, lwd=4)
text((-6-psL)/2+psL, PLCmax+8, expression(PLC[max]), cex=5)

legend("topright", c("Original curve", "Scenario 2"),
       lty=c(1, 2, 2, 2), col=c("black", "blue"), cex=3, bty="n")
box()

dev.copy2pdf(file = "Figures/Scenario figure.pdf")
