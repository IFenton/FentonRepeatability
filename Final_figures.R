# Repeatability analysis figures
# Project: Repeatability
# Author: Isabel Fenton
# Date created: 3/8/2018
# Date last edited: 3/8/2018
# 
# Code for figures of the repeatability paper
# Fenton et al on foram repeatability
# 
# Previous file: Repeatability.R

library(sjPlot) # plotting glms

# Figure 1
png("Docs/Final/fig1a.png")
with(people.df, boxplot(fracCorrSp ~ Experienced, ylim = c(0, 1), xlim = c(-0.5, 2.5), xaxt = "n", yaxt = "n", boxwex = 0.4, col = c("blue3", "purple2"), cex.axis = 1.5, las = 1))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSp, add = TRUE, at = 0, boxwex = 0.8, col = "red2", xaxt = "n", yaxt = "n"))
axis(1, c(0,1,2), labels = FALSE, cex.axis = 1.5)
axis(2, seq(0,1,.2), labels = FALSE, cex.axis = 1.5)
dev.off()

png("Docs/Final/fig1b.png")
plot(NULL, ylim = c(0, 1), xlim = c(0, 4.8), yaxt = "n", xaxt = "n", cex.axis = 1.5, xlab = "", ylab = "", las = 2)
with(people.df, boxplot(fracCorrSN ~ Experienced, boxwex = 0.2, col = c("steelblue3", "plum2"), at = c(2.4, 3.9), add = TRUE, xaxt = "n", yaxt = "n", xlim = c(0, 4.8)))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSN, add = TRUE, at = 0.9, boxwex = 0.4, col = "indianred2", xaxt = "n", yaxt = "n"))
axis(1, c(0.6,2.1,3.6), labels = FALSE, cex.axis = 1.5)
axis(2, seq(0,1,.2), labels = FALSE, cex.axis = 1.5)

with(people.df, boxplot(fracCorrSM ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("royalblue3", "darkorchid2"), cex.axis = 1.5, at = c(2.1, 3.6), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSM, add = TRUE, at = 0.6, boxwex = 0.4, col = "firebrick2", xaxt = "n", yaxt = "n"))

with(people.df, boxplot(fracCorrSY ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("blue4", "purple4"), cex.axis = 1.5, at = c(1.8, 3.3), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSY, add = TRUE, at = 0.3, boxwex = 0.4, col = "red4", xaxt = "n", yaxt = "n"))

rect(4.65, .964, 4.8, .992, col = "purple4")
rect(4.65, .914, 4.8, .942, col = "darkorchid2")
rect(4.65, .864, 4.8, .892, col = "plum2")
rect(4.4, .964, 4.55, .992, col = "blue4")
rect(4.4, .914, 4.55, .942, col = "royalblue3")
rect(4.4, .864, 4.55, .892, col = "steelblue3")
rect(4.15, .964, 4.3, .992, col = "red4")
rect(4.15, .914, 4.3, .942, col = "firebrick2")
rect(4.15, .864, 4.3, .892, col = "indianred2")
dev.off()

png("Docs/Final/fig1c.png")
plot(NULL, ylim = c(0, 1), xlim = c(0, 4.8), yaxt = "n", xaxt = "n", cex.axis = 1.5, xlab = "", ylab = "", las = 1)
with(people.df, boxplot(fracCorrY ~ Experienced, boxwex = 0.2, col = c("blue4", "purple4"), at = c(1.8, 3.3), add = TRUE, xaxt = "n", yaxt = "n", xlim = c(0, 4.8)))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrY, add = TRUE, at = 0.3, boxwex = 0.4, col = "red4", xaxt = "n", yaxt = "n"))
axis(1, c(0.6,2.1,3.6), labels = FALSE, cex.axis = 1.5)
axis(2, seq(0,1,.2), labels = FALSE, cex.axis = 1.5)

with(people.df, boxplot(fracCorrM ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("royalblue3", "darkorchid2"), cex.axis = 1.5, at = c(2.1, 3.6), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrM, add = TRUE, at = 0.6, boxwex = 0.4, col = "firebrick2", xaxt = "n", yaxt = "n"))

with(people.df, boxplot(fracCorrNwoUID ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("steelblue3", "plum2"), cex.axis = 1.5, at = c(2.4, 3.9), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrNwoUID, add = TRUE, at = 0.9, boxwex = 0.4, col = "indianred2", xaxt = "n", yaxt = "n"))

rect(4.65, .964, 4.8, .992, col = "purple4")
rect(4.65, .914, 4.8, .942, col = "darkorchid2")
rect(4.65, .864, 4.8, .892, col = "plum2")
rect(4.4, .964, 4.55, .992, col = "blue4")
rect(4.4, .914, 4.55, .942, col = "royalblue3")
rect(4.4, .864, 4.55, .892, col = "steelblue3")
rect(4.15, .964, 4.3, .992, col = "red4")
rect(4.15, .914, 4.3, .942, col = "firebrick2")
rect(4.15, .864, 4.3, .892, col = "indianred2")
dev.off()

png("Docs/Final/fig1d.png")
plot(NULL, bty = "l", ylim = c(0, 1), xlim = c(0, 4.8), yaxt = "n", xaxt = "n", cex.axis = 1.5, xlab = "", ylab = "", las = 2)
with(people.df, boxplot(fracCorr125 ~ Experienced, boxwex = 0.2, col = c("steelblue3", "plum2"), at = c(1.8, 3.3), add = TRUE, xaxt = "n", yaxt = "n", xlim = c(0, 4.8)))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorr125, add = TRUE, at = 0.3, boxwex = 0.4, col = "indianred2", xaxt = "n", yaxt = "n"))
axis(1, c(0.6,2.1,3.6), labels = FALSE, cex.axis = 1.5)
axis(2, seq(0,1,.2), labels = FALSE, cex.axis = 1.5)

with(people.df, boxplot(fracCorr200 ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("royalblue3", "darkorchid2"), cex.axis = 1.5, at = c(2.1, 3.6), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorr200, add = TRUE, at = 0.6, boxwex = 0.4, col = "firebrick2", xaxt = "n", yaxt = "n"))

with(people.df, boxplot(fracCorr400 ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("blue4", "purple4"), cex.axis = 1.5, at = c(2.4, 3.9), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorr400, add = TRUE, at = 0.9, boxwex = 0.4, col = "red4", xaxt = "n", yaxt = "n"))

rect(4.65, .089, 4.8, .117, col = "purple4")
rect(4.65, .039, 4.8, .067, col = "darkorchid2")
rect(4.65, -.011, 4.8, .017, col = "plum2")
rect(4.4, .089, 4.55, .117, col = "blue4")
rect(4.4, .039, 4.55, .067, col = "royalblue3")
rect(4.4, -.011, 4.55, .017, col = "steelblue3")
rect(4.15, .089, 4.3, .117, col = "red4")
rect(4.15, .039, 4.3, .067, col = "firebrick2")
rect(4.15, -.011, 4.3, .017, col = "indianred2")

dev.off()


# Figure 2
png("Docs/Final/fig2.png", 264.6, 185.2, units = "mm", res = 200)
conf_mat(IDs.long, "ID", "DefinitiveID", axes.same = FALSE, abb.end = c("juvenile", "nonmacro", "unIDd"), sp.exc = "lost", xlab = "Individual ID", grid = TRUE)
dev.off()


# Figure 3
png("Docs/Final/fig3a.png", 160, 160, units = "mm", res = 600)
set_theme(theme_classic(), axis.title.size = 1.4, axis.textsize = 1.2, panel.gridcol = "white")
p <- sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "Taught"), title = "", show.ci = TRUE, facet.grid = FALSE, axis.title = c(expression(paste("Mean Diameter / ", mu, "m")), "Percentage Correct"), geom.size = 1.2, point.alpha = 0.5, prnt.plot = FALSE)
p$plot + scale_x_continuous(breaks = c(-0.9975478, 0.3995266, 1.216763, 1.796601), labels = c(200, 400, 600, 800))
dev.off()
rm(p)

png("Docs/Final/fig3b.png", 160, 160, units = "mm", res = 600)
set_theme(theme_classic(), axis.title.size = 1.4, axis.textsize = 1.2, panel.gridcol = "white")
p2 <- sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "HowLong"), title = "", show.ci = TRUE, facet.grid = FALSE, axis.title = c(expression(paste("Mean Diameter / ", mu, "m")), "Percentage Correct"), geom.size = 1.2, point.alpha = 0.5, prnt.plot = FALSE)
p2$plot + scale_x_continuous(breaks = c(-0.9975478, 0.3995266, 1.216763, 1.796601), labels = c(200, 400, 600, 800))
dev.off()
rm(p2)


png("Docs/Final/fig3c.png", 160, 160, units = "mm", res = 600)
set_theme(theme_classic(), axis.title.size = 1.4, axis.textsize = 1.2, panel.gridcol = "white")
p1 <- sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "Experience"), title = "", show.ci = TRUE, facet.grid = FALSE, axis.title = c(expression(paste("Mean Diameter / ", mu, "m")), "Percentage Correct"), geom.size = 1.2, point.alpha = 0.5, prnt.plot = FALSE)
p1$plot + scale_x_continuous(breaks = c(-0.9975478, 0.3995266, 1.216763, 1.796601), labels = c(200, 400, 600, 800))
dev.off()
rm(p1)
