  library(plotrix)
  m1 <- read.csv("SM24_m1.csv", header = F, sep = "\t", dec = ",")
  m2 <- read.csv("SM24_m2.csv", header = F, sep = "\t", dec = ",")
  m3 <- read.csv("SM24_m3.csv", header = F, sep = "\t", dec = ",")
  m4 <- read.csv("SM24_m4.csv", header = F, sep = "\t", dec = ",")
  m5 <- read.csv("SM24_m5.csv", header = F, sep = "\t", dec = ",")
  #m6 <- read.csv("SM24_m6.csv", header = F, sep = "\t", dec = ",") # Cells are no longer growing subconfluently (plateau has been reached)
  
  mutant_2_raw <- cbind(m1[3:6, 3],m2[3:6, 3],m3[3:6, 3],m4[3:6, 3],m5[3:6, 3])
  mutant_6_raw<- cbind(m1[3:6, 5],m2[3:6, 5],m3[3:6, 5],m4[3:6, 5],m5[3:6, 5])
  mutant_11_raw<- cbind(m1[3:6, 7],m2[3:6, 7],m3[3:6, 7],m4[3:6, 7],m5[3:6, 7])
  mutant_12_raw<- cbind(m1[3:6, 9],m2[3:6, 9],m3[3:6, 9],m4[3:6, 9],m5[3:6, 9])
  wt_3_raw<- cbind(m1[3:6, 4],m2[3:6, 4],m3[3:6, 4],m4[3:6, 4],m5[3:6, 4])
  wt_4_raw<- cbind(m1[3:6,6],m2[3:6, 6],m3[3:6, 6],m4[3:6, 6],m5[3:6, 6])
  wt_5_raw<- cbind(m1[3:6, 8],m2[3:6, 8],m3[3:6, 8],m4[3:6, 8],m5[3:6, 8])
  wt_14_raw<- cbind(m1[3:6, 10],m2[3:6, 10],m3[3:6, 10],m4[3:6, 10],m5[3:6, 10])
  
  mutant_2_raw[4, 5] <- NA # Experimental mistake (too much Cell Titer Glo reagent)
  
  time <- c(0, 9, 24,33,47.5)
  
  mutant_2_m_raw <- apply(mutant_2_raw, 2, mean, na.rm = TRUE)
  mutant_6_m_raw<- apply(mutant_6_raw, 2, mean, na.rm = TRUE)
  mutant_11_m_raw<- apply(mutant_11_raw, 2, mean, na.rm = TRUE)
  mutant_12_m_raw<- apply(mutant_12_raw, 2, mean, na.rm = TRUE)
  wt_3_m_raw<- apply(wt_3_raw, 2, mean, na.rm = TRUE)
  wt_4_m_raw<- apply(wt_4_raw, 2, mean, na.rm = TRUE)
  wt_5_m_raw<- apply(wt_5_raw, 2, mean, na.rm = TRUE)
  wt_14_m_raw<- apply(wt_14_raw, 2, mean, na.rm = TRUE)
  
  
  #mutant_2 <- mutant_2_raw
  #mutant_6<- mutant_6_raw
  #mutant_11<- mutant_11_raw
  #mutant_12<- mutant_12_raw
  #wt_3<- wt_3_raw
  #wt_4<- wt_4_raw
  #wt_5<- wt_5_raw
  #wt_14<- wt_14_raw
  mutant_2 <- mutant_2_raw / mutant_2_m_raw[1]
  mutant_6<- mutant_6_raw / mutant_6_m_raw[1]
  mutant_11<- mutant_11_raw / mutant_11_m_raw[1]
  mutant_12<- mutant_12_raw / mutant_12_m_raw[1]
  wt_3<- wt_3_raw / wt_3_m_raw[1]
  wt_4<- wt_4_raw / wt_4_m_raw[1]
  wt_5<- wt_5_raw / wt_5_m_raw[1]
  wt_14<- wt_14_raw / wt_14_m_raw[1]
  
  mutant_2_m <- apply(log(mutant_2) / log(2), 2, mean, na.rm = TRUE)
  mutant_6_m<- apply(log(mutant_6) / log(2), 2, mean, na.rm = TRUE)
  mutant_11_m<- apply(log(mutant_11) / log(2), 2, mean, na.rm = TRUE)
  mutant_12_m<- apply(log(mutant_12) / log(2), 2, mean, na.rm = TRUE)
  wt_3_m<- apply(log(wt_3) / log(2), 2, mean, na.rm = TRUE)
  wt_4_m<- apply(log(wt_4) / log(2), 2, mean, na.rm = TRUE)
  wt_5_m<-apply(log(wt_5) / log(2), 2, mean, na.rm = TRUE)
  wt_14_m<- apply(log(wt_14) / log(2), 2, mean, na.rm = TRUE)
  
  mutant_2_se <- apply(log(mutant_2) / log(2), 2, std.error, na.rm = TRUE)
  mutant_6_se<- apply(log(mutant_6) / log(2), 2, std.error, na.rm = TRUE)
  mutant_11_se<- apply(log(mutant_11) / log(2), 2, std.error, na.rm = TRUE)
  mutant_12_se<- apply(log(mutant_12) / log(2), 2, std.error, na.rm = TRUE)
  wt_3_se<- apply(log(wt_3) / log(2), 2, std.error, na.rm = TRUE)
  wt_4_se<- apply(log(wt_4) / log(2), 2, std.error, na.rm = TRUE)
  wt_5_se<-apply(log(wt_5) / log(2), 2, std.error, na.rm = TRUE)
  wt_14_se<- apply(log(wt_14) / log(2), 2, std.error, na.rm = TRUE)
  
  #y_range <- range(pretty(range(c(wt_5_m + wt_5_se, wt_5_m - wt_5_se,mutant_12_m + mutant_12_se, mutant_12_m -mutant_12_se), na.rm = T)))
  y_range=c(0,3.5)
 x_range <- c(0, 50)
  pdf('Proliferation_HCT116_KO_WT.pdf',width=6,height=6)
  plotCI(time, wt_3_m,wt_3_se, pch = 20, col = "#636363", ylim = y_range, xlim = x_range, xlab = "Time (h)", ylab = "Log-transformed number of cells (relative to initial value)")
  par(new = T)
  plotCI(time, wt_4_m, wt_4_se, pch = 20, col = "#969696", ylim = y_range, xlim = x_range, xlab = "", ylab = "", axes = F)
  par(new = T)
  plotCI(time, wt_5_m, wt_5_se, pch = 20, col = "#bdbdbd", ylim = y_range, xlim = x_range, xlab = "", ylab = "", axes = F)
  par(new = T)
  plotCI(time, wt_14_m, wt_14_se, pch = 20, col = "#252525", ylim = y_range, xlim = x_range, xlab = "", ylab = "", axes = F)
  par(new = T)
  plotCI(time, mutant_2_m, mutant_2_se, pch = 20, col = "#de2d26", ylim = y_range, xlim = x_range, xlab = "", ylab = "", axes = F)
  par(new = T)
  plotCI(time, mutant_6_m, mutant_6_se, pch = 20, col = "#fb6a4a", ylim = y_range, xlim = x_range, xlab = "", ylab = "", axes = F)
  par(new = T)
  plotCI(time, mutant_11_m, mutant_11_se, pch = 20, col = "#fc9272", ylim = y_range, xlim = x_range, xlab = "", ylab = "", axes = F)
  par(new = T)
  plotCI(time, mutant_12_m, mutant_12_se, pch = 20, col = "#a50f15", ylim = y_range, xlim = x_range, xlab = "", ylab = "", axes = F)
  
  mutant_2_for_lm <- array(log(mutant_2) / log(2), dim = c(1, length(mutant_2)))[1, ]
  mutant_6_for_lm<- array(log(mutant_6) / log(2), dim = c(1, length(mutant_6)))[1, ]
  mutant_11_for_lm<- array(log(mutant_11) / log(2), dim = c(1, length(mutant_11)))[1, ]
  mutant_12_for_lm<- array(log(mutant_12) / log(2), dim = c(1, length(mutant_12)))[1, ]
  wt_3_for_lm<- array(log(wt_3) / log(2), dim = c(1, length(wt_3)))[1, ]
  wt_4_for_lm<- array(log(wt_4) / log(2), dim = c(1, length(wt_4)))[1, ]
  wt_5_for_lm<-array(log(wt_5) / log(2), dim = c(1, length(wt_5)))[1, ]
  wt_14_for_lm<- array(log(wt_14) / log(2), dim = c(1, length(wt_14)))[1, ]
  
  time_for_lm <- c()
  for (t in time) {
    time_for_lm <- append(time_for_lm, rep(t, 4))
  }
  
  #mutant_2_lm <- lm(mutant_2_for_lm ~ time_for_lm)
  #mutant_6_lm<-lm(mutant_6_for_lm ~ time_for_lm)
  #mutant_11_lm<- lm(mutant_11_for_lm ~ time_for_lm)
  #mutant_12_lm<- lm(mutant_12_for_lm ~ time_for_lm)
  #wt_3_lm<- lm(wt_3_for_lm ~ time_for_lm)
  #wt_4_lm<- lm(wt_4_for_lm ~ time_for_lm)
  #wt_5_lm<-lm(wt_5_for_lm ~ time_for_lm)
  #wt_14_lm<-lm(wt_14_for_lm ~ time_for_lm)
  #pred.mutant_2 <- data.frame(time_for_lm)
  #pp.mutant_2 <- predict(mutant_2_lm, int = "c", newdata = pred.mutant_2)
  #matlines(time_for_lm, pp.mutant_2, col = "red", lty = c(2, 3, 3))
  #pred.mutant_6 <- data.frame(time_for_lm)
  #pp.mutant_6 <- predict(mutant_6_lm, int = "c", newdata = pred.mutant_6)
  #matlines(time_for_lm, pp.mutant_6, col = "red", lty = c(2, 3, 3))
  #pred.mutant_11 <- data.frame(time_for_lm)
  #pp.mutant_11 <- predict(mutant_11_lm, int = "c", newdata = pred.mutant_11)
  #matlines(time_for_lm, pp.mutant_11, col = "red", lty = c(2, 3, 3))
  #pred.mutant_12 <- data.frame(time_for_lm)
  #pp.mutant_12 <- predict(mutant_12_lm, int = "c", newdata = pred.mutant_12)
  #matlines(time_for_lm, pp.mutant_12, col = "red", lty = c(2, 3, 3))
  #pred.wt_3 <- data.frame(time_for_lm)
  #pp.wt_3 <- predict(wt_3_lm, int = "c", newdata = pred.wt_3)
  #matlines(time_for_lm, pp.wt_3, col = "black", lty = c(2, 3, 3))
  #pred.wt_4 <- data.frame(time_for_lm)
  #pp.wt_4 <- predict(wt_4_lm, int = "c", newdata = pred.wt_4)
  #matlines(time_for_lm, pp.wt_4, col = "black", lty = c(2, 3, 3))
  #pred.wt_5 <- data.frame(time_for_lm)
  #pp.wt_5 <- predict(wt_5_lm, int = "c", newdata = pred.wt_5)
  #matlines(time_for_lm, pp.wt_5, col = "black", lty = c(2, 3, 3))
  #pred.wt_14 <- data.frame(time_for_lm)
  #pp.wt_14 <- predict(wt_14_lm, int = "c", newdata = pred.wt_14)
  #matlines(time_for_lm, pp.wt_14, col = "black", lty = c(2, 3, 3))
  
  wt_for_lm=c(wt_3_for_lm,wt_4_for_lm,wt_5_for_lm,wt_14_for_lm)
  mutant_for_lm=c(mutant_2_for_lm,mutant_6_for_lm,mutant_11_for_lm,mutant_12_for_lm)
  time_for_popul_lm=rep(time_for_lm,4)
  
  y <- c(wt_for_lm, mutant_for_lm)
  x <- c(time_for_popul_lm, time_for_popul_lm)
  genotype <- as.factor(c(rep("wt", times = length(wt_for_lm)), rep("mutant", times = length(mutant_for_lm))))
  summary(lm(y ~ x * genotype))
  
  wt_lm=lm(wt_for_lm ~ time_for_popul_lm)
  mutant_lm=lm(mutant_for_lm ~ time_for_popul_lm)
  pred.wt <- data.frame(time_for_popul_lm)
  pp.wt <- predict(wt_lm, int = "c", newdata = pred.wt)
#matlines(unique(time_for_popul_lm[seq(1,length(time_for_popul_lm),by=4)]),unique(pp.wt[seq(1,nrow(pp.wt),by=4),]), col = "black", lty = c(0, 2, 2))
  polygon(c(unique(time_for_popul_lm[seq(1,length(time_for_popul_lm),by=4)]),rev(unique(time_for_popul_lm[seq(1,length(time_for_popul_lm),by=4)]))),c(unique(pp.wt[seq(1,nrow(pp.wt),by=4),2]),rev(unique(pp.wt[seq(1,nrow(pp.wt),by=4),3]))), col = rgb(0, 0, 0, 0.15),border=NA)
  abline(wt_lm,col='black',lty=2,lwd=2)
  pred.mutant <- data.frame(time_for_popul_lm)
  pp.mutant <- predict(mutant_lm, int = "c", newdata = pred.mutant)
 # matlines(unique(time_for_popul_lm[seq(1,length(time_for_popul_lm),by=4)]), unique(pp.mutant[seq(1,nrow(pp.wt),by=4),]), col = "red", lty = c(0, 2, 2))
  polygon(c(unique(time_for_popul_lm[seq(1,length(time_for_popul_lm),by=4)]),rev(unique(time_for_popul_lm[seq(1,length(time_for_popul_lm),by=4)]))),c(unique(pp.mutant[seq(1,nrow(pp.mutant),by=4),2]),rev(unique(pp.mutant[seq(1,nrow(pp.mutant),by=4),3]))), col = rgb(1, 0, 0, 0.15),border=NA)
  abline(mutant_lm,col='red',lty=2,lwd=2)
  
  double_time_wt <- as.numeric(1 / wt_lm$coefficients[2])
  double_time_wt # 37.19992 +/- 0.1791
  double_time_mutant <- as.numeric(1 / mutant_lm$coefficients[2])
  double_time_mutant # 42.86441 +/- 0.1822
  
  #legend("bottomright", c(paste("wt clones",round(double_time_wt,2),'h',sep=' '), paste("miR-34a KO clones",round(double_time_mutant,2),'h',sep=' ')), col = c("black", "red"), pch = 20, cex = 1)
  
  dev.off()
