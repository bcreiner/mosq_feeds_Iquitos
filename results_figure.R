
V_vec <- seq(1,10,by = .5)
D_vec <- seq(1,5, by = .1)
I_vec <- seq(6,16,by = 1)


Pred <- array(NA,dim=c(length(V_vec), length(D_vec), length(I_vec), 3))

for (i in 1:length(V_vec)){
  for (j in 1:length(D_vec)){
    for (k in 1:length(I_vec)){
      if (I_vec[k] <= 8.5 | I_vec[k] >= 11.5){
        tmpval <- ifelse(I_vec[k] < 10, 1,2)
        tmpdat <- data.frame(viremia = V_vec[i], DoI = D_vec[j], bin_EIP = tmpval)
        tmp <- predict(mod_to_plot, newdata = tmpdat, se = TRUE)
        Pred[i,j,k,1] <- ilogit(tmp$fit - 1.96 * tmp$se.fit)
        Pred[i,j,k,2] <- ilogit(tmp$fit)
        Pred[i,j,k,3] <- ilogit(tmp$fit + 1.96 * tmp$se.fit)
      }
    }
  }
}

wholeV <- which(V_vec %% 1 == 0)
wholeD <- which(D_vec %% 1 == 0)
wholeI <- which(I_vec %% 1 == 0)

tmpdata <- which(df_D2$DoI == 1)

VBINS <- 1:10
VColor <- brewer.pal(length(VBINS) - 1, "Reds")
DColor <- rev(brewer.pal(6, "BuGn"))
IColor <- brewer.pal(length(I_vec) - 1, "Purples")
VColLoc <- sapply(df_D2$viremia, function(x) findInterval(x,VBINS, all.inside = TRUE))


YLIM <- c(0, max(df_D2$fracPos, max(Pred, na.rm=TRUE)))

pdf(file="output/results_v2.pdf",height=14,width=12)

par(mfrow=c(2,2), mar=c(5.1,4.1,1.1,1.1))
plot(jitter(df_D2$viremia), df_D2$fracPos, type='n',ylim= YLIM, 
     xlim=range(VBINS), pch=19, cex = df_D2$num_tested / 5, col=VColor[VColLoc], axes=FALSE, ann=FALSE)
box()
axis(1, VBINS)
axis(2,pretty(YLIM), paste0(100*pretty(YLIM),"%"))
mtext("Log viremia", 1, line=2.6)
mtext("Transmission probability", 2, line=2.6)

LegBins2 <- seq(.3, YLIM[2], length = length(wholeD) + 2)+.025
XGap2 <- .25
XStart2 <- 3
YGap2 <- .005


LegBins1 <- seq(.3, YLIM[2] - .1, length = 3)+.025
XStart1 <- 1.5
XGap1 <- .45
tmpVal <- c(5,10,25)
points(rep(XStart1,3), LegBins1, cex = tmpVal / 5, pch=21)
text(rep(XStart1+XGap1,3), LegBins1, tmpVal)
text(XStart1,  YLIM[2] - .025 , "Sample
size")


for (i in 1:length(wholeD)){
  lines(V_vec,Pred[,wholeD[i],7,2],col=DColor[i], lwd=2)
  text(XStart2 + XGap2, LegBins2[i], pos=4, i)
  rect(XStart2,LegBins2[i] - YGap2, XStart2 + XGap2, LegBins2[i] + YGap2, col=DColor[i])
}
text(XStart2 + 3.75*XGap2, LegBins2[3], "Day of infection", srt = 90)
points(df_D2$viremia, df_D2$fracPos, pch=21, cex = df_D2$num_tested / 5, col=1, bg = DColor[df_D2$DoI])

###
###
###

plot(jitter(df_D2$DoI), df_D2$fracPos, type='n',ylim= YLIM, xlim=c(.75,6), pch=19, cex = df_D2$num_tested / 5, col=VColor[VColLoc], axes=FALSE, ann=FALSE)
box()
axis(1,1:5)
axis(2,pretty(YLIM), paste0(100*pretty(YLIM),"%"))
mtext("Day of illness", 1, line=2.6)
mtext("Transmission probability", 2, line=2.6)



LegBins1 <- seq(.3, YLIM[2] - .1, length = 3)+.025
XStart1 <- 5
XGap1 <- .3
tmpVal <- c(5,10,25)
points(rep(XStart1,3), LegBins1, cex = tmpVal / 5, pch=21)
text(rep(XStart1+XGap1,3), LegBins1, tmpVal)
text(XStart1,  YLIM[2] - .025 , "Sample
size")

LegBins2 <- seq(.1, YLIM[2], length = length(wholeV) + 1)+.025
XGap2 <- .2
XStart2 <- 5.5


for (i in 1:length(wholeV)){
  lines(D_vec,Pred[wholeV[i],,7,2],col=VColor[i], lwd=2)
  text(XStart2 + XGap2, LegBins2[i], pos=4, VBINS[i])
  if (i < length(wholeV)){
    rect(XStart2,LegBins2[i], XStart2 + XGap2, LegBins2[i+1], col=VColor[i])
  }
}
text(XStart2 + 2.3*XGap2, mean(LegBins2), "Log viremia", srt = 90)
points(jitter(df_D2$DoI), df_D2$fracPos, pch=21, cex = df_D2$num_tested / 5, col=1, bg = VColor[VColLoc])

###
###
###


plot(jitter(df_D2$EIP), df_D2$fracPos, type='n',ylim= YLIM, xlim=c(5,20), pch=19, cex = df_D2$num_tested / 5, col=VColor[VColLoc], axes=FALSE, ann=FALSE)
box()
axis(1,6:16)
axis(2,pretty(YLIM), paste0(100*pretty(YLIM),"%"))
mtext("Extrinsic incubation period", 1, line=2.6)
mtext("Transmission probability", 2, line=2.6)



LegBins1 <- seq(.3, YLIM[2] - .1, length = 3)+.025
XStart1 <- 16.5
XGap1 <- .75
tmpVal <- c(5,10,25)
points(rep(XStart1,3), LegBins1, cex = tmpVal / 5, pch=21)
text(rep(XStart1+XGap1,3), LegBins1, tmpVal)
text(XStart1,  YLIM[2] - .025 , "Sample
size")

LegBins2 <- seq(.1, YLIM[2], length = length(wholeV) + 1)+.025
XGap2 <- .5
XStart2 <- 18


for (i in 1:length(wholeV)){
  lines(I_vec,Pred[wholeV[i],wholeD[2],,2],col=VColor[i], lwd=2)
  text(XStart2 + XGap2, LegBins2[i], pos=4, VBINS[i])
  if (i < length(wholeV)){
    rect(XStart2,LegBins2[i], XStart2 + XGap2, LegBins2[i+1], col=VColor[i])
  }
}
text(XStart2 + 3*XGap2, mean(LegBins2), "Log viremia", srt = 90)
points(jitter(df_D2$EIP), df_D2$fracPos, pch=21, cex = df_D2$num_tested / 5, col=1, bg = VColor[VColLoc])

###
###
###

tmp_YLIM <- range(VBINS)
plot(jitter(df_D2$DoI), df_D2$viremia, type='n',ylim= tmp_YLIM, xlim=c(.75,6), pch=19, cex = df_D2$num_tested / 5, col=VColor[VColLoc], axes=FALSE, ann=FALSE)
box()
axis(1,1:5)
axis(2,pretty(tmp_YLIM))
mtext("Day of illness", 1, line=2.6)
mtext("Log viremia", 2, line=2.6)

LegBins1 <- seq(tmp_YLIM[2] - 3, tmp_YLIM[2] - 1, length = 3)+.025
XStart1 <- 5
XGap1 <- .3
tmpVal <- c(5,10,25)
points(rep(XStart1,3), LegBins1, cex = tmpVal / 5, pch=21)
text(rep(XStart1+XGap1,3), LegBins1, tmpVal)
text(XStart1,  tmp_YLIM[2] - .075 , "Sample
size")

LegBins2 <- seq(1, tmp_YLIM[2], length = length(wholeV) + 1)+.025
XGap2 <- .25
XStart2 <- 5.5


for (i in 1:length(wholeV)){
  lines(D_vec,Pred[wholeV[i],,7,2],col=VColor[i], lwd=2)
  text(XStart2 + XGap2, LegBins2[i], pos=4, VBINS[i])
  if (i < length(wholeV)){
    rect(XStart2,LegBins2[i], XStart2 + XGap2, LegBins2[i+1], col=VColor[i])
  }
}
text(XStart2 + 2.25*XGap2, mean(LegBins2), "Log viremia", srt = 90)
points(jitter(df_D2$DoI), df_D2$viremia, pch=21, cex = df_D2$num_tested / 5, col=1, bg = VColor[VColLoc])



dev.off()
