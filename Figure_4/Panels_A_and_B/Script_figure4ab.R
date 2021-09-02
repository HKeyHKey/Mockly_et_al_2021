library(plotrix)
#setwd("~/Documents/HCT_116_miR34a_duplex_Assay_6/")
## Import dataset
m1 <- read.csv("SM43_m1.csv", header = F, sep = ",")
m1[7, 4:5] <- NA # Cell seeding mistake
m2 <- read.csv("SM43_m2.csv", header = F, sep = ",")
m3 <- read.csv("SM43_m3.csv", header = F, sep = ",")
m4 <- read.csv("SM43_m4.csv", header = F, sep = ",")
m5 <- read.csv("SM43_m5.csv", header = F, sep = ",")
m6 <- read.csv("SM43_m6.csv", header = F, sep = ",")

time <- c(6, 22, 31, 46.5, 55, 70.25)

## 100 nM transfection of duplexes
# Please note that in the paper, only 10 and 1 nM transfections results are shown.

ctrl_100_raw <- cbind(as.numeric(c(m1[2:7, 6])), as.numeric(c(m2[2:7, 6])), as.numeric(c(m3[2:7, 6])), as.numeric(c(m4[2:7, 6])), as.numeric(c(m5[2:7, 6])), as.numeric(c(m6[2:7, 6])))
miR34a_100_raw <- cbind(as.numeric(c(m1[2:7, 9])), as.numeric(c(m2[2:7, 9])), as.numeric(c(m3[2:7, 9])), as.numeric(c(m4[2:7, 9])), as.numeric(c(m5[2:7, 9])), as.numeric(c(m6[2:7, 9])))

ctrl_100_m_raw <- apply(ctrl_100_raw, 2, mean, na.rm = TRUE)
miR34a_100_m_raw <- apply(miR34a_100_raw, 2, mean, na.rm = TRUE)

# Normalized by t= 0 h
ctrl_100 <- ctrl_100_raw / ctrl_100_m_raw[1]
miR34a_100 <- miR34a_100_raw / miR34a_100_m_raw[1]

## log2-transformation
ctrl_100_m <- apply(log(ctrl_100) / log(2), 2, mean, na.rm = TRUE)
miR34a_100_m <- apply(log(miR34a_100) / log(2), 2, mean, na.rm = TRUE)

ctrl_100_se <- apply(log(ctrl_100) / log(2), 2, std.error, na.rm = TRUE)
miR34a_100_se <- apply(log(miR34a_100) / log(2), 2, std.error, na.rm = TRUE)

## Cell proliferation plot
# y_range <- range(pretty(range(c(ctrl_100_m + ctrl_100_se, ctrl_100_m - ctrl_100_se, miR34a_100_m + miR34a_100_se, miR34a_100_m - miR34a_100_se), na.rm = T)))
y_range <- c(0, 2.5)
x_range <- c(5, 72)
# pdf("Proliferation_HCT_100_mimic.pdf", width = 6, height = 6)
plotCI(time, ctrl_100_m, ctrl_100_se, pch = 20, col = "black", ylim = y_range, xlim = x_range, xlab = "Time (h)", ylab = "Log-transformed number of cells (relative to initial value)")
par(new = T)
plotCI(time, miR34a_100_m, miR34a_100_se, pch = 20, col = "forestgreen", ylim = y_range, xlim = x_range, xlab = "", ylab = "", axes = F)

## Representation of the 95% confidence interval
ctrl_100_for_lm <- array(log(ctrl_100) / log(2), dim = c(1, length(ctrl_100)))[1, ]
miR34a_100_for_lm <- array(log(miR34a_100) / log(2), dim = c(1, length(miR34a_100)))[1, ]
time_for_lm <- c()

for (t in time) {
  time_for_lm <- append(time_for_lm, rep(t, 6))
}

ctrl_100_lm <- lm(ctrl_100_for_lm ~ time_for_lm)
miR34a_100_lm <- lm(miR34a_100_for_lm ~ time_for_lm)

pred.ctrl_100 <- data.frame(time_for_lm)
pp.ctrl_100 <- predict(ctrl_100_lm, int = "c", newdata = pred.ctrl_100)
polygon(c(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]), rev(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]))), c(unique(pp.ctrl_100[seq(1, nrow(pp.ctrl_100), by = 6), 2]), rev(unique(pp.ctrl_100[seq(1, nrow(pp.ctrl_100), by = 6), 3]))), col = rgb(0, 0, 0, 0.15), border = NA)
abline(ctrl_100_lm, col = "black", lty = 2, lwd = 2)
# matlines(time_for_lm,pp.ctrl_100,col='black',lty=c(2,3,3))
pred.miR34a_100 <- data.frame(time_for_lm)
pp.miR34a_100 <- predict(miR34a_100_lm, int = "c", newdata = pred.miR34a_100)
polygon(c(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]), rev(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]))), c(unique(pp.miR34a_100[seq(1, nrow(pp.miR34a_100), by = 6), 2]), rev(unique(pp.miR34a_100[seq(1, nrow(pp.miR34a_100), by = 6), 3]))), col = rgb(0.6, 1, 0.6, 0.3), border = NA)
abline(miR34a_100_lm, col = "forestgreen", lty = 2, lwd = 2)
# matlines(time_for_lm,pp.miR34a_100,col='red',lty=c(2,3,3))

## ANOVA 100 nM control duplex vs 100 nM miR-34a duplex
y <- c(ctrl_100_for_lm, miR34a_100_for_lm)
x <- c(time_for_lm, time_for_lm)
genotype <- as.factor(c(rep("ctrl_100", times = length(ctrl_100_for_lm)), rep("miR34a_100", times = length(miR34a_100_for_lm))))
model <- lm(y ~ x * genotype)
summary.aov(model)
# Doubling time
double_time_ctrl_100 <- as.numeric(1 / ctrl_100_lm$coefficients[2])
double_time_miR34a_100 <- as.numeric(1 / miR34a_100_lm$coefficients[2])

legend("topleft", c(paste("Control duplex (100 nM) - Doubling time = ", round(double_time_ctrl_100, 1), " h", sep = ""), paste("miR-34a duplex (100 nM) - Doubling time = ", round(double_time_miR34a_100, 1), " h", sep = "")), col = c("black", "forestgreen"), pch = 20, cex = 0.8)

# dev.off()


## 10 nM transfection of duplexes
# Please note that in the paper, this result is shown on figure 3A.

ctrl_10_raw <- cbind(as.numeric(c(m1[2:7, 5])), as.numeric(c(m2[2:7, 5])), as.numeric(c(m3[2:7, 5])), as.numeric(c(m4[2:7, 5])), as.numeric(c(m5[2:7, 5])), as.numeric(c(m6[2:7, 5])))
miR34a_10_raw <- cbind(as.numeric(c(m1[2:7, 8])), as.numeric(c(m2[2:7, 8])), as.numeric(c(m3[2:7, 8])), as.numeric(c(m4[2:7, 8])), as.numeric(c(m5[2:7, 8])), as.numeric(c(m6[2:7, 8])))

ctrl_10_m_raw <- apply(ctrl_10_raw, 2, mean, na.rm = TRUE)
miR34a_10_m_raw <- apply(miR34a_10_raw, 2, mean, na.rm = TRUE)

ctrl_10 <- ctrl_10_raw / ctrl_10_m_raw[1]
miR34a_10 <- miR34a_10_raw / miR34a_10_m_raw[1]

ctrl_10_m <- apply(log(ctrl_10) / log(2), 2, mean, na.rm = TRUE)
miR34a_10_m <- apply(log(miR34a_10) / log(2), 2, mean, na.rm = TRUE)

ctrl_10_se <- apply(log(ctrl_10) / log(2), 2, std.error, na.rm = TRUE)
miR34a_10_se <- apply(log(miR34a_10) / log(2), 2, std.error, na.rm = TRUE)

 pdf("Proliferation_HCT_10_mimic.pdf", width = 6, height = 6)
# y_range <- range(pretty(range(c(ctrl_10_m + ctrl_10_se, ctrl_10_m - ctrl_10_se, miR34a_10_m + miR34a_10_se, miR34a_10_m - miR34a_10_se), na.rm = T)))
y_range <- c(0, 2.5)
x_range <- c(5, 72)
plotCI(time, ctrl_10_m, ctrl_10_se, pch = 20, col = "black", ylim = y_range, xlim = x_range, xlab = "Time (h)", ylab = "Log-transformed number of cells (relative to initial value)")
par(new = T)
plotCI(time, miR34a_10_m, miR34a_10_se, pch = 20, col = "forestgreen", ylim = y_range, xlim = x_range, xlab = "", ylab = "", axes = F)

ctrl_10_for_lm <- array(log(ctrl_10) / log(2), dim = c(1, length(ctrl_10)))[1, ]
miR34a_10_for_lm <- array(log(miR34a_10) / log(2), dim = c(1, length(miR34a_10)))[1, ]
time_for_lm <- c()

for (t in time) {
  time_for_lm <- append(time_for_lm, rep(t, 6))
}

ctrl_10_lm <- lm(ctrl_10_for_lm ~ time_for_lm)
miR34a_10_lm <- lm(miR34a_10_for_lm ~ time_for_lm)

pred.ctrl_10 <- data.frame(time_for_lm)
pp.ctrl_10 <- predict(ctrl_10_lm, int = "c", newdata = pred.ctrl_10)
polygon(c(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]), rev(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]))), c(unique(pp.ctrl_10[seq(1, nrow(pp.ctrl_10), by = 6), 2]), rev(unique(pp.ctrl_10[seq(1, nrow(pp.ctrl_10), by = 6), 3]))), col = rgb(0, 0, 0, 0.15), border = NA)
abline(ctrl_10_lm, col = "black", lty = 2, lwd = 2)
# matlines(time_for_lm,pp.ctrl_10,col='black',lty=c(2,3,3))
pred.miR34a_10 <- data.frame(time_for_lm)
pp.miR34a_10 <- predict(miR34a_10_lm, int = "c", newdata = pred.miR34a_10)
polygon(c(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]), rev(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]))), c(unique(pp.miR34a_10[seq(1, nrow(pp.miR34a_10), by = 6), 2]), rev(unique(pp.miR34a_10[seq(1, nrow(pp.miR34a_10), by = 6), 3]))), col = rgb(0.6, 1, 0.6, 0.3), border = NA)
abline(miR34a_10_lm, col = "forestgreen", lty = 2, lwd = 2)
# matlines(time_for_lm,pp.miR34a_10,col='red',lty=c(2,3,3))

## ANOVA 10 nM control duplex vs 10 nM miR-34a duplex
y <- c(ctrl_10_for_lm, miR34a_10_for_lm)
x <- c(time_for_lm, time_for_lm)
genotype <- as.factor(c(rep("ctrl_10", times = length(ctrl_10_for_lm)), rep("miR34a_10", times = length(miR34a_10_for_lm))))
model <- lm(y ~ x * genotype)
summary.aov(model)
# Doubling time
double_time_ctrl_10 <- as.numeric(1 / ctrl_10_lm$coefficients[2])
double_time_miR34a_10 <- as.numeric(1 / miR34a_10_lm$coefficients[2])

legend("bottomright", c(paste("Control duplex (10 nM) - Doubling time = ", round(double_time_ctrl_10, 1), " h", sep = ""), paste("miR-34a duplex (10 nM) - Doubling time = ", round(double_time_miR34a_10, 1), " h", sep = "")), col = c("black", "forestgreen"), pch = 20, cex = 0.8)

 dev.off()



## 1 nM transfection of duplexes
# Please note that in the paper, this result is shown on figure 3B.

ctrl_1_raw <- cbind(as.numeric(c(m1[2:7, 4])), as.numeric(c(m2[2:7, 4])), as.numeric(c(m3[2:7, 4])), as.numeric(c(m4[2:7, 4])), as.numeric(c(m5[2:7, 4])), as.numeric(c(m6[2:7, 4])))
miR34a_1_raw <- cbind(as.numeric(c(m1[2:7, 7])), as.numeric(c(m2[2:7, 7])), as.numeric(c(m3[2:7, 7])), as.numeric(c(m4[2:7, 7])), as.numeric(c(m5[2:7, 7])), as.numeric(c(m6[2:7, 7])))

ctrl_1_m_raw <- apply(ctrl_1_raw, 2, mean, na.rm = TRUE)
miR34a_1_m_raw <- apply(miR34a_1_raw, 2, mean, na.rm = TRUE)

ctrl_1 <- ctrl_1_raw / ctrl_1_m_raw[1]
miR34a_1 <- miR34a_1_raw / miR34a_1_m_raw[1]

ctrl_1_m <- apply(log(ctrl_1) / log(2), 2, mean, na.rm = TRUE)
miR34a_1_m <- apply(log(miR34a_1) / log(2), 2, mean, na.rm = TRUE)

ctrl_1_se <- apply(log(ctrl_1) / log(2), 2, std.error, na.rm = TRUE)
miR34a_1_se <- apply(log(miR34a_1) / log(2), 2, std.error, na.rm = TRUE)

 pdf("Proliferation_HCT_1_mimic.pdf", width = 6, height = 6)
# y_range <- range(pretty(range(c(ctrl_1_m + ctrl_1_se, ctrl_1_m - ctrl_1_se, miR34a_1_m + miR34a_1_se, miR34a_1_m - miR34a_1_se), na.rm = T)))
y_range <- c(0, 2.5)
x_range <- c(5, 72)
plotCI(time, ctrl_1_m, ctrl_1_se, pch = 20, col = "black", ylim = y_range, xlim = x_range, xlab = "Time (h)", ylab = "Log-transformed number of cells (relative to initial value)")
par(new = T)
plotCI(time, miR34a_1_m, miR34a_1_se, pch = 20, col = "forestgreen", ylim = y_range, xlim = x_range, xlab = "", ylab = "", axes = F)

ctrl_1_for_lm <- array(log(ctrl_1) / log(2), dim = c(1, length(ctrl_1)))[1, ]
miR34a_1_for_lm <- array(log(miR34a_1) / log(2), dim = c(1, length(miR34a_1)))[1, ]
time_for_lm <- c()

for (t in time) {
  time_for_lm <- append(time_for_lm, rep(t, 6))
}

ctrl_1_lm <- lm(ctrl_1_for_lm ~ time_for_lm)
miR34a_1_lm <- lm(miR34a_1_for_lm ~ time_for_lm)

pred.ctrl_1 <- data.frame(time_for_lm)
pp.ctrl_1 <- predict(ctrl_1_lm, int = "c", newdata = pred.ctrl_1)
polygon(c(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]), rev(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]))), c(unique(pp.ctrl_1[seq(1, nrow(pp.ctrl_1), by = 6), 2]), rev(unique(pp.ctrl_1[seq(1, nrow(pp.ctrl_1), by = 6), 3]))), col = rgb(0, 0, 0, 0.15), border = NA)
abline(ctrl_1_lm, col = "black", lty = 2, lwd = 2)
# matlines(time_for_lm,pp.ctrl_1,col='black',lty=c(2,3,3))
pred.miR34a_1 <- data.frame(time_for_lm)
pp.miR34a_1 <- predict(miR34a_1_lm, int = "c", newdata = pred.miR34a_1)
polygon(c(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]), rev(unique(time_for_lm[seq(1, length(time_for_lm), by = 6)]))), c(unique(pp.miR34a_1[seq(1, nrow(pp.miR34a_1), by = 6), 2]), rev(unique(pp.miR34a_1[seq(1, nrow(pp.miR34a_1), by = 6), 3]))), col = rgb(0.6, 1, 0.6, 0.3), border = NA)
abline(miR34a_1_lm, col = "forestgreen", lty = 2, lwd = 2)
# matlines(time_for_lm,pp.miR34a_1,col='red',lty=c(2,3,3))

## ANOVA 1 nM control duplex vs 1 nM miR-34a duplex
y <- c(ctrl_1_for_lm, miR34a_1_for_lm)
x <- c(time_for_lm, time_for_lm)
genotype <- as.factor(c(rep("ctrl_1", times = length(ctrl_1_for_lm)), rep("miR34a_1", times = length(miR34a_1_for_lm))))
model <- lm(y ~ x * genotype)
summary.aov(model)
# Doubling time
double_time_ctrl_1 <- as.numeric(1 / ctrl_1_lm$coefficients[2])
double_time_miR34a_1 <- as.numeric(1 / miR34a_1_lm$coefficients[2])

legend("bottomright", c(paste("Control duplex (1 nM) - Doubling time = ", round(double_time_ctrl_1, 1), "h", sep = ""), paste("miR-34a duplex (1 nM) - Doubling time = ", round(double_time_miR34a_1, 1), "h", sep = "")), col = c("black", "forestgreen"), pch = 20, cex = 0.8)

 dev.off()
