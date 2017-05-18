## repeatability

x <- read.csv("../../Box Sync/Forams - adaptive zones/IFwork/Project/Repeatability/Data/sizes_unscaled.csv")
head(x)

png("../../Box Sync/Forams - adaptive zones/IFwork/Project/Repeatability/Figures/size_hist.png", 400, 400)
hist(x$X1000, main = "", xlab = expression(paste("Max. diameter / ", mu, "m")), col = "maroon3")
dev.off()