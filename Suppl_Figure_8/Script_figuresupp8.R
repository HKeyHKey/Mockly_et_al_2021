library("plotrix")
library("readxl")
# Assay 1 (20220225)
X20220225_flowjo <- read_excel("20220225_flowjo.xls", range = "B1:D81")
X20220225_flowjo_data <- X20220225_flowjo[-c(1:20), ]
# Assay 2 (20220228)
X20220228_flowjo_data <- read_excel("20220228_flowjo.xls", range = "B1:D121")

table <- data.frame(sample = character(), alive = integer(), dead = integer(), apoptosis = integer(), other = integer(), stringsAsFactors = FALSE)
for (i in 0:11) {
  table[i + 1, 1] <- gsub("2022-..-.._(.*)_.000..mqd", "\\1", as.character(X20220225_flowjo_data[1 + (i * 5), 1]))
  table[i + 1, 2] <- as.numeric(X20220225_flowjo_data[5 + (i * 5), 2])
  table[i + 1, 3] <- as.numeric(X20220225_flowjo_data[3 + (i * 5), 2])
  table[i + 1, 4] <- as.numeric(X20220225_flowjo_data[4 + (i * 5), 2])
  table[i + 1, 5] <- as.numeric(X20220225_flowjo_data[2 + (i * 5), 2])
}
for (i in 0:23) {
  table[i + 13, 1] <- gsub("2022-..-.._(.*)_.000..mqd", "\\1", as.character(X20220228_flowjo_data[1 + (i * 5), 1]))
  table[i + 13, 2] <- as.numeric(X20220228_flowjo_data[5 + (i * 5), 2])
  table[i + 13, 3] <- as.numeric(X20220228_flowjo_data[3 + (i * 5), 2])
  table[i + 13, 4] <- as.numeric(X20220228_flowjo_data[4 + (i * 5), 2])
  table[i + 13, 5] <- as.numeric(X20220228_flowjo_data[2 + (i * 5), 2])
}

gamme_concentration <- c(1e-12, 10e-12, 100e-12, 1e-9, 10e-9, 100e-9)

# pdf("FACS_apoptosis_plots.pdf", width = 9, height = 3)
par(mfrow = c(1, 3))

# DEAD CELLS#
table_dead <- rbind(
  table$dead[which(table$sample == "1pm_ctrl" | table$sample == "1pM_ctrl")],
  table$dead[which(table$sample == "10pm_ctrl" | table$sample == "10pM_ctrl")],
  table$dead[which(table$sample == "100pm_ctrl" | table$sample == "100pM_ctrl")],
  table$dead[which(table$sample == "1nm_ctrl" | table$sample == "1nM_ctrl")],
  table$dead[which(table$sample == "10nm_ctrl" | table$sample == "10nM_ctrl")],
  table$dead[which(table$sample == "100nm_ctrl" | table$sample == "100nM_ctrl")],
  table$dead[which(table$sample == "1pm_miR34" | table$sample == "1pM_miR34")],
  table$dead[which(table$sample == "10pm_miR34" | table$sample == "10pM_miR34")],
  table$dead[which(table$sample == "100pm_miR34" | table$sample == "100pM_miR34")],
  table$dead[which(table$sample == "1nm_miR34" | table$sample == "1nM_miR34")],
  table$dead[which(table$sample == "10nm_miR34" | table$sample == "10nM_miR34")],
  table$dead[which(table$sample == "100nm_miR34" | table$sample == "100nM_miR34")]
)

anova_dead <- data.frame(
  var = c(table_dead[1, ], table_dead[2, ], table_dead[3, ], table_dead[4, ], table_dead[5, ], table_dead[6, ], table_dead[7, ], table_dead[8, ], table_dead[9, ], table_dead[10, ], table_dead[11, ], table_dead[12, ]),
  treatment = factor(c(rep("ctrl", 18), rep("miR34", 18))),
  concentration = factor(rep(c(rep("1pM", 3), rep("10pM", 3), rep("100pM", 3), rep("1nM", 3), rep("10nM", 3), rep("100nM", 3)), 2))
)
fit_dead <- aov(var ~ treatment + concentration, anova_dead)
summary(fit_dead)
# treatment p-value=0.365 Not significant
# concentration p-value=0.465 Not significant

plotCI(gamme_concentration,
  apply(table_dead[1:6, ], 1, mean), apply(table_dead[1:6, ], 1, std.error), apply(table_dead[1:6, ], 1, std.error),
  ylim = c(0, 100),
  col = "black", scol = "black", main = "Dead cells", xlab = "duplex concentration (M)", ylab = "Percentage of dead cells", pch = 16, log = "x"
)

axis(1, labels = T, c(1e-11, 1e-9, 1e-7))
axis(1, labels = F, lwd.ticks = 0.5, c(
  9e-12, 8e-12, 7e-12, 6e-12, 5e-12, 4e-12, 3e-12, 2e-12,
  9e-11, 8e-11, 7e-11, 6e-11, 5e-11, 4e-11, 3e-11, 2e-11,
  9e-10, 8e-10, 7e-10, 6e-10, 5e-10, 4e-10, 3e-10, 2e-10,
  9e-9, 8e-9, 7e-9, 6e-9, 5e-9, 4e-9, 3e-9, 2e-9,
  9e-8, 8e-8, 7e-8, 6e-8, 5e-8, 4e-8, 3e-8, 2e-8
))
par(new = TRUE)
plotCI(gamme_concentration,
  apply(table_dead[7:12, ], 1, mean), apply(table_dead[7:12, ], 1, std.error), apply(table_dead[7:12, ], 1, std.error),
  ylim = c(0, 100),
  xlab = "", ylab = "", axes = FALSE, col = "forestgreen", scol = "forestgreen", pch = 16, log = "x"
)
legend("topright",
  legend = c("mimic control", "mimic miR-34a"),
  col = c("black", "forestgreen"), lty = 1, cex = 0.8
)

# APOTOTIC CELLS#
table_apoptotic <- rbind(
  table$apoptosis[which(table$sample == "1pm_ctrl" | table$sample == "1pM_ctrl")],
  table$apoptosis[which(table$sample == "10pm_ctrl" | table$sample == "10pM_ctrl")],
  table$apoptosis[which(table$sample == "100pm_ctrl" | table$sample == "100pM_ctrl")],
  table$apoptosis[which(table$sample == "1nm_ctrl" | table$sample == "1nM_ctrl")],
  table$apoptosis[which(table$sample == "10nm_ctrl" | table$sample == "10nM_ctrl")],
  table$apoptosis[which(table$sample == "100nm_ctrl" | table$sample == "100nM_ctrl")],
  table$apoptosis[which(table$sample == "1pm_miR34" | table$sample == "1pM_miR34")],
  table$apoptosis[which(table$sample == "10pm_miR34" | table$sample == "10pM_miR34")],
  table$apoptosis[which(table$sample == "100pm_miR34" | table$sample == "100pM_miR34")],
  table$apoptosis[which(table$sample == "1nm_miR34" | table$sample == "1nM_miR34")],
  table$apoptosis[which(table$sample == "10nm_miR34" | table$sample == "10nM_miR34")],
  table$apoptosis[which(table$sample == "100nm_miR34" | table$sample == "100nM_miR34")]
)

anova_apoptosis <- data.frame(
  var = c(table_apoptotic[1, ], table_apoptotic[2, ], table_apoptotic[3, ], table_apoptotic[4, ], table_apoptotic[5, ], table_apoptotic[6, ], table_apoptotic[7, ], table_apoptotic[8, ], table_apoptotic[9, ], table_apoptotic[10, ], table_apoptotic[11, ], table_apoptotic[12, ]),
  treatment = factor(c(rep("ctrl", 18), rep("miR34", 18))),
  concentration = factor(rep(c(rep("1pM", 3), rep("10pM", 3), rep("100pM", 3), rep("1nM", 3), rep("10nM", 3), rep("100nM", 3)), 2))
)
fit_apoptosis <- aov(var ~ treatment + concentration, anova_apoptosis)
summary(fit_apoptosis)
# treatment p-value=0.104 Not significant
# concentration p-value=1.02e-07 Significant
fit_apoptosis <- aov(var ~ treatment * concentration, anova_apoptosis)
summary(fit_apoptosis)
# treatment p-value=0.106 Not significant
# concentration p-value=4.17e-07 Significant
# treatment:concentration p-value=0.423 Not significant, so the additive model is indicated in the paper.

plotCI(gamme_concentration,
  apply(table_apoptotic[1:6, ], 1, mean), apply(table_apoptotic[1:6, ], 1, std.error), apply(table_apoptotic[1:6, ], 1, std.error),
  ylim = c(0, 100),
  col = "black", scol = "black", main = "Apoptotic cells", xlab = "duplex concentration (M)", ylab = "Percentage of apoptotic cells", pch = 16, log = "x"
)
axis(1, labels = T, c(1e-11, 1e-9, 1e-7))
axis(1, labels = F, lwd.ticks = 0.5, c(
  9e-12, 8e-12, 7e-12, 6e-12, 5e-12, 4e-12, 3e-12, 2e-12,
  9e-11, 8e-11, 7e-11, 6e-11, 5e-11, 4e-11, 3e-11, 2e-11,
  9e-10, 8e-10, 7e-10, 6e-10, 5e-10, 4e-10, 3e-10, 2e-10,
  9e-9, 8e-9, 7e-9, 6e-9, 5e-9, 4e-9, 3e-9, 2e-9,
  9e-8, 8e-8, 7e-8, 6e-8, 5e-8, 4e-8, 3e-8, 2e-8
))
par(new = TRUE)
plotCI(gamme_concentration,
  apply(table_apoptotic[7:12, ], 1, mean), apply(table_apoptotic[7:12, ], 1, std.error), apply(table_apoptotic[7:12, ], 1, std.error),
  ylim = c(0, 100),
  xlab = "", ylab = "", axes = FALSE, col = "forestgreen", scol = "forestgreen", pch = 16, log = "x"
)
legend("topright",
  legend = c("mimic control", "mimic miR-34a"),
  col = c("black", "forestgreen"), lty = 1, cex = 0.8
)

# ALIVE CELLS#
table_alive <- rbind(
  table$alive[which(table$sample == "1pm_ctrl" | table$sample == "1pM_ctrl")],
  table$alive[which(table$sample == "10pm_ctrl" | table$sample == "10pM_ctrl")],
  table$alive[which(table$sample == "100pm_ctrl" | table$sample == "100pM_ctrl")],
  table$alive[which(table$sample == "1nm_ctrl" | table$sample == "1nM_ctrl")],
  table$alive[which(table$sample == "10nm_ctrl" | table$sample == "10nM_ctrl")],
  table$alive[which(table$sample == "100nm_ctrl" | table$sample == "100nM_ctrl")],
  table$alive[which(table$sample == "1pm_miR34" | table$sample == "1pM_miR34")],
  table$alive[which(table$sample == "10pm_miR34" | table$sample == "10pM_miR34")],
  table$alive[which(table$sample == "100pm_miR34" | table$sample == "100pM_miR34")],
  table$alive[which(table$sample == "1nm_miR34" | table$sample == "1nM_miR34")],
  table$alive[which(table$sample == "10nm_miR34" | table$sample == "10nM_miR34")],
  table$alive[which(table$sample == "100nm_miR34" | table$sample == "100nM_miR34")]
)

anova_alive <- data.frame(
  var = c(table_alive[1, ], table_alive[2, ], table_alive[3, ], table_alive[4, ], table_alive[5, ], table_alive[6, ], table_alive[7, ], table_alive[8, ], table_alive[9, ], table_alive[10, ], table_alive[11, ], table_alive[12, ]),
  treatment = factor(c(rep("ctrl", 18), rep("miR34", 18))),
  concentration = factor(rep(c(rep("1pM", 3), rep("10pM", 3), rep("100pM", 3), rep("1nM", 3), rep("10nM", 3), rep("100nM", 3)), 2))
)
fit_alive <- aov(var ~ treatment + concentration, anova_alive)
summary(fit_alive)
# treatment p-value= 0.128706 Not significant
# concentration p-value=0.000584 Significant
fit_alive <- aov(var ~ treatment * concentration, anova_alive)
summary(fit_alive)
# treatment p-value=0.101185 Not significant
# concentration p-value=0.000303 Significant
# treatment:concentration p-value=0.101627 Not significant, so the additive model is indicated in the paper.

plotCI(gamme_concentration,
  apply(table_alive[1:6, ], 1, mean), apply(table_alive[1:6, ], 1, std.error), apply(table_alive[1:6, ], 1, std.error),
  ylim = c(0, 100),
  col = "black", scol = "black", main = "Alive cells", xlab = "duplex concentration (M)", ylab = "Percentage of alive cells", pch = 16, log = "x"
)
axis(1, labels = T, c(1e-11, 1e-9, 1e-7))
axis(1, labels = F, lwd.ticks = 0.5, c(
  9e-12, 8e-12, 7e-12, 6e-12, 5e-12, 4e-12, 3e-12, 2e-12,
  9e-11, 8e-11, 7e-11, 6e-11, 5e-11, 4e-11, 3e-11, 2e-11,
  9e-10, 8e-10, 7e-10, 6e-10, 5e-10, 4e-10, 3e-10, 2e-10,
  9e-9, 8e-9, 7e-9, 6e-9, 5e-9, 4e-9, 3e-9, 2e-9,
  9e-8, 8e-8, 7e-8, 6e-8, 5e-8, 4e-8, 3e-8, 2e-8
))
par(new = TRUE)
plotCI(gamme_concentration,
  apply(table_alive[7:12, ], 1, mean), apply(table_alive[7:12, ], 1, std.error), apply(table_alive[7:12, ], 1, std.error),
  ylim = c(0, 100),
  xlab = "", ylab = "", axes = FALSE, col = "forestgreen", scol = "forestgreen", pch = 16, log = "x"
)
legend("topright",
  legend = c("mimic control", "mimic miR-34a"),
  col = c("black", "forestgreen"), lty = 1, cex = 0.8
)
# dev.off()