# Repeatability analysis
# Project: Repeatability
# Author: Isabel Fenton
# Date created: 9/3/2017
# Date last edited: 6/9/2017
# 
# Code for analysis of the repeatability results. For more details see Lab book Repeatability.docx)
# 
# Previous file: size histogram.R
# Next file:

# Inputs ------------------------------------------------------------------
# CompleteIDs.csv
# SpeciesData.csv
# PersonData.csv
# sizes_unscaled.csv - the initial measurements
# measurements2.csv - the final (accurate) measurements

# Outputs -----------------------------------------------------------------

# Source files / libraries ------------------------------------------------
library(cluster) # for the dendrogram
library(caret) # for the confusion matrix
library(colorRamps) # colours
library(ade4) # for distinctiveness
library(sjPlot) # plotting glms
library(lme4)
library(stringr) # for confusion matrix axes
library(broom) # alternative plotting
library(ggplot2) # alternative plotting
library(car) # p-values in the anova
library(MuMIn) # r2 values for the model

rm(list = ls())

source("../../../../../Dropbox/Documents/AdrianaDePalma/plotLmerMeans.R")
source("Code/plotLmerMeansRepeat.R")

dev.off()
par.def <- par()


# 1. Load in the data -----------------------------------------------------
SpeciesData <- read.csv("Data/SpeciesData.csv")
str(SpeciesData)

completeIDs <- read.csv("Data/CompleteIDs.csv")
str(completeIDs)

personData <- read.csv("Data/PersonData.csv", check.names = FALSE) # to keep the spaces etc. in the column names
str(personData)

sp.size <- read.csv("Data/MeasurementsFull.csv")
str(sp.size)

# set the levels of each ID column to the SpeciesData
for (i in grep("ID", names(completeIDs), value = TRUE)) {
  completeIDs[,i] <- factor(completeIDs[,i], levels = unique(c(levels(SpeciesData$Species), levels(SpeciesData$ChecklistGenus))))
}
rm(i)
str(completeIDs)

# 2. Species information  ---------------------------------------------

# 2a. Species / genus list ------------------------------------------------
species <- grep(" ", SpeciesData$Species, value = TRUE)
genera <- unique(grep("^[[:upper:]]", SpeciesData$ChecklistGenus, value = TRUE))

# 2b. Functional distinctiveness ------------------------------------------
str(SpeciesData)

# remove the NAs (those that aren't extant species)
summary(SpeciesData)
SpeciesData.nna <- SpeciesData[complete.cases(SpeciesData), 6:10]
rownames(SpeciesData.nna) <- SpeciesData[complete.cases(SpeciesData), 4]
rownames(SpeciesData.nna)[rownames(SpeciesData.nna) == "ruber (pink)"] <- "ruber_pink"
summary(SpeciesData.nna)
head(SpeciesData.nna)

# create a functional dendrogram
sp.dist <- daisy(SpeciesData.nna)
sp.dend <- hclust(sp.dist)
png("Figures/FD_dend.png", 800, 800)
plot(sp.dend) # n.b. this looks convincing checking a few - the polytomies seem to have identical traits
dev.off()

# calculate the distinctiveness
sp.fd <- originality(hclust2phylog(sp.dend), c(1:4, 6:7))
summary(sp.fd)
str(sp.fd)

SpeciesData.nna <- merge(SpeciesData.nna, sp.fd, by = 0)
str(SpeciesData.nna)

# compare one metric with the dendrogram
par(mar = c(8.1, 4.1, 2.1, 2.1), mfrow = c(2, 1))
plot(sp.dend) 
with(SpeciesData.nna, barplot(VW, names.arg = Row.names, las = 2, col = "blue4"))
par(mar = c(5.1, 4.1, 2.1, 2.1), mfrow = c(1,1))

# compare all the different metrics
par(mar = c(8.1, 4.1, 2.1, 2.1), mfrow = c(6, 1))
with(SpeciesData.nna, barplot(VW, names.arg = Row.names, las = 2, col = "blue4"))
with(SpeciesData.nna, barplot(M, names.arg = Row.names, las = 2, col = "blue4"))
with(SpeciesData.nna, barplot(SpeciesData.nna[, "NWU*"], names.arg = Row.names, las = 2, col = "blue4"))
with(SpeciesData.nna, barplot(NWW, names.arg = Row.names, las = 2, col = "blue4"))
with(SpeciesData.nna, barplot(ED, names.arg = Row.names, las = 2, col = "blue4"))
with(SpeciesData.nna, barplot(eqsplit, names.arg = Row.names, las = 2, col = "blue4"))
par(mar = c(5.1, 4.1, 2.1, 2.1), mfrow = c(1,1))

# ED seems to have the best properties, so add that to the dataframe
SpeciesData <- merge(SpeciesData, SpeciesData.nna[, names(SpeciesData.nna) %in% c("Row.names", "ED")], by.x = "ChecklistSpecies", by.y = "Row.names", all = TRUE)
SpeciesData$ED[SpeciesData$ChecklistSpecies == "ruber (pink)"] <- SpeciesData.nna$ED[(SpeciesData.nna$Row.names == "ruber_pink")]
SpeciesData <- SpeciesData[SpeciesData$ChecklistSpecies != "ruber_pink", ]

# 3. Person level data -------------------------------

# 3a. Create dataframe ----------------------------------------------
# list of people in the analysis
people <- gsub("ID", "", grep("ID", names(completeIDs), value = TRUE))
people.col <- grep("ID", names(completeIDs), value = TRUE)
p.conf.col <- grep("C$", names(completeIDs), value = TRUE)

people.df <- data.frame(Name = people)
people.df$Experienced <- factor(c(rep(NA, 2), rep("Experienced", 4), rep("Student", nrow(people.df) - 6)))
str(people.df)

# 3b. Person level statistics -----------------------------
## All individuals 
# Number of species identified
people.df$NumSpID <- apply(completeIDs[, people.col], 2, function(x) sum(unique(x) %in% species))

# Number of genera identified
unique(grep("^[[:upper:]]", unlist(strsplit(as.character(completeIDs$DefinitiveID), " ")), value = TRUE))
people.df$NumGenID <- apply(completeIDs[, people.col], 2, function(x) length(unique(grep("^[[:upper:]]", unlist(strsplit(as.character(x), " ")), value = TRUE))))

# Number of specimens lost
sum(completeIDs$PooleID != "lost")
sum(completeIDs$PooleID != "lost"& completeIDs$PooleID %in% SpeciesData$Species)
people.df$numLost <- apply(completeIDs[, people.col], 2, function(x) sum(x == "lost"))
people.df$TotalIndSp <- apply(completeIDs[, people.col], 2, function(x) sum(x != "lost" & x %in% SpeciesData$Species))
people.df$TotalIndGe <- apply(completeIDs[, people.col], 2, function(x) sum(x != "lost"))

# Fraction of species identified
people.df$numUnIDd <- apply(completeIDs[, people.col], 2, function(x) sum(x == "unIDd"))

people.df$fracSpecIDd <- people.df$NumSpID / people.df$NumSpID[people.df$Name == "Definitive"]

# Fraction of genera identified
people.df$fracGenIDd <- people.df$NumGenID / people.df$NumGenID[people.df$Name == "Definitive"]


## Percentage accuracy per person
sum(completeIDs[,grep("Fenton", names(completeIDs))[1]] == completeIDs$DefinitiveID) / sum(completeIDs[,grep("Fenton", names(completeIDs))[1]] != "lost")

# species
people.df$fracCorrSp <- apply(completeIDs[, people.col], 2, function(x) sum(x == completeIDs$DefinitiveID) / sum(x != "lost"))
mean(people.df$fracCorrSp)
median(people.df$fracCorrSp[!is.na(people.df$Experienced)])
median(people.df$fracCorrSp[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrSp[people.df$Experienced == "Student"], na.rm = TRUE)

# genus
gsub(" .*", "", completeIDs$DefinitiveID)

sum(gsub(" .*", "", completeIDs$FentonID) == gsub(" .*", "", completeIDs$DefinitiveID)) / sum(completeIDs[,grep("Fenton", names(completeIDs))[1]] != "lost")

people.df$fracCorrGen <- apply(completeIDs[, people.col], 2, function(x) sum(gsub(" .*", "", x) == gsub(" .*", "", completeIDs$DefinitiveID)) / sum(x != "lost"))
mean(people.df$fracCorrGen)

people.df$numUnID <- NA
people.df$numUnID <- apply(completeIDs[, people.col], 2, function(x) sum(x == "unIDd", na.rm = TRUE))

# box plot of percentage accuracy split by students and Experienced
png("Figures/Percent accuracy.png")
with(people.df, boxplot(fracCorrSp ~ Experienced, ylim = c(0, 1), xlim = c(-0.5, 2.5), xaxt = "n", boxwex = 0.4, col = c("blue3", "purple2"), cex.axis = 1.5, las = 1))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSp, add = TRUE, at = 0, boxwex = 0.8, col = "red2", xaxt = "n", yaxt = "n"))
axis(1, c(0,1,2), c("Full", levels(people.df$Experienced)), cex.axis = 1.5)
text(0:2, 0.95, paste("(", c(sum(!is.na(people.df$Experienced)), sum(people.df$Experienced == "Experienced", na.rm = TRUE), sum(people.df$Experienced == "Student", na.rm = TRUE)), ")", sep = ""))
dev.off()

median(people.df$fracCorrSp[!is.na(people.df$Experienced)])
median(people.df$fracCorrSp[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrSp[people.df$Experienced == "Student"], na.rm = TRUE)

## Split by confidence 
# number of identifications by confidence (working at specimen level)
sum(completeIDs$FentonC == "y", na.rm = TRUE)

people.df$numCy <- NA
people.df$numCy[3:nrow(people.df)] <- apply(completeIDs[, p.conf.col], 2, function(x) sum(x == "y", na.rm = TRUE))

people.df$numCm <- NA
people.df$numCm[3:nrow(people.df)] <- apply(completeIDs[, p.conf.col], 2, function(x) sum(x == "m", na.rm = TRUE))

people.df$numCn <- NA
people.df$numCn[3:nrow(people.df)] <- apply(completeIDs[, p.conf.col], 2, function(x) sum(x == "n", na.rm = TRUE))

people.df

# fraction of confident identifications
sum(completeIDs[, people.col[3]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[1]] == "y", na.rm = TRUE) / sum(completeIDs[, p.conf.col[1]] == "y", na.rm = TRUE)

sum(completeIDs[, people.col[3]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[1]] == "m", na.rm = TRUE) / sum(completeIDs[, p.conf.col[1]] == "m", na.rm = TRUE)

sum(completeIDs[, people.col[3]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[1]] == "n", na.rm = TRUE) / sum(completeIDs[, p.conf.col[1]] == "n", na.rm = TRUE)

people.df$fracCorrY <- NA
people.df$fracCorrM <- NA
people.df$fracCorrN <- NA
people.df$fracCorrNwoUID <- NA
for (i in 1:length(p.conf.col)){
  people.df$fracCorrY[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[i]] == "y", na.rm = TRUE) / sum(completeIDs[, p.conf.col[i]] == "y", na.rm = TRUE)
  people.df$fracCorrM[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[i]] == "m", na.rm = TRUE) / sum(completeIDs[, p.conf.col[i]] == "m", na.rm = TRUE)
  people.df$fracCorrN[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[i]] == "n", na.rm = TRUE) / sum(completeIDs[, p.conf.col[i]] == "n", na.rm = TRUE)
  people.df$fracCorrNwoUID[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[i]] == "n" & completeIDs[, people.col[2+i]] != "unIDd", na.rm = TRUE) / sum(completeIDs[, p.conf.col[i]] == "n" & completeIDs[, people.col[2+i]] != "unIDd", na.rm = TRUE, na.rm = TRUE)
}
# where there are no specimens classified as 'n', this is currently giving NaN
people.df[is.na(people.df)] <- NA

hist(people.df$fracCorrY)
hist(people.df$fracCorrM)
hist(people.df$fracCorrN)
hist(people.df$fracCorrNwoUID)

# percentage accuracy by confidence
png("Figures/Percent accuracy conf.png")
plot(NULL, ylim = c(0, 1), xlim = c(0, 4.8), xaxt = "n", cex.axis = 1.5, xlab = "", ylab = "", las = 1)
with(people.df, boxplot(fracCorrY ~ Experienced, boxwex = 0.2, col = c("blue3", "purple2"), at = c(1.8, 3.3), add = TRUE, xaxt = "n", yaxt = "n", xlim = c(0, 4.8)))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrY, add = TRUE, at = 0.3, boxwex = 0.4, col = "red2", xaxt = "n", yaxt = "n"))
axis(1, c(0.6,2.1,3.6), c("Full", levels(people.df$Experienced)), cex.axis = 1.5)

with(people.df, boxplot(fracCorrM ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("royalblue3", "darkorchid2"), cex.axis = 1.5, at = c(2.1, 3.6), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrM, add = TRUE, at = 0.6, boxwex = 0.4, col = "firebrick2", xaxt = "n", yaxt = "n"))

with(people.df, boxplot(fracCorrNwoUID ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("steelblue3", "plum2"), cex.axis = 1.5, at = c(2.4, 3.9), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrNwoUID, add = TRUE, at = 0.9, boxwex = 0.4, col = "indianred2", xaxt = "n", yaxt = "n"))

legend("topright", c("yes", "maybe", "no"), fill = c("black", "grey50", "grey80"), bty = "n", cex = 1.2)
dev.off()

# percent correct if yes
median(people.df$fracCorrY[!is.na(people.df$Experienced)])
median(people.df$fracCorrY[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrY[people.df$Experienced == "Student"], na.rm = TRUE)

# percent correct if maybe
median(people.df$fracCorrM[!is.na(people.df$Experienced)])
median(people.df$fracCorrM[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrM[people.df$Experienced == "Student"], na.rm = TRUE)

# percent correct if no
median(people.df$fracCorrNwoUID[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorrNwoUID[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrNwoUID[people.df$Experienced == "Student"], na.rm = TRUE)


## accuracy based on species level confidence
# generate columns
people.df$numcSY <- NA
people.df$numcSM <- NA
people.df$numcSN <- NA
people.df$fracCorrSY <- NA
people.df$fracCorrSM <- NA
people.df$fracCorrSN <- NA

sum(completeIDs[, people.col[2+5]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[5+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y", na.rm = TRUE) / sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[5+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y" & completeIDs[, people.col[2+5]] != "lost", na.rm = TRUE)

for (i in 1:(length(people.col) - 2)){
  # sum (correct ID and confident in species) / sum(confident in species and not lost)
  people.df$numcSY[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y" & completeIDs[, people.col[2+i]] != "lost", na.rm = TRUE)
  people.df$fracCorrSY[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y", na.rm = TRUE) / people.df$numcSY[2 + i]
  people.df$numcSM[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "m" & completeIDs[, people.col[2+i]] != "lost", na.rm = TRUE)
  people.df$fracCorrSM[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "m", na.rm = TRUE) / people.df$numcSM[2 + i] 
  people.df$numcSN[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "n" & completeIDs[, people.col[2+i]] != "lost", na.rm = TRUE) 
  people.df$fracCorrSN[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "n", na.rm = TRUE) / people.df$numcSN[2 + i]
}

png("Figures/Percent accuracy species level.png")
plot(NULL, ylim = c(0, 1), xlim = c(0, 4.8), xaxt = "n", cex.axis = 1.5, xlab = "", ylab = "", las = 2)
with(people.df, boxplot(fracCorrSN ~ Experienced, boxwex = 0.2, col = c("steelblue3", "plum2"), at = c(2.4, 3.9), add = TRUE, xaxt = "n", yaxt = "n", xlim = c(0, 4.8)))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSN, add = TRUE, at = 0.9, boxwex = 0.4, col = "indianred2", xaxt = "n", yaxt = "n"))
axis(1, c(0.6,2.1,3.6), c("Full", levels(people.df$Experienced)), cex.axis = 1.5)

with(people.df, boxplot(fracCorrSM ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("royalblue3", "darkorchid2"), cex.axis = 1.5, at = c(2.1, 3.6), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSM, add = TRUE, at = 0.6, boxwex = 0.4, col = "firebrick2", xaxt = "n", yaxt = "n"))

with(people.df, boxplot(fracCorrSY ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("blue3", "purple2"), cex.axis = 1.5, at = c(1.8, 3.3), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSY, add = TRUE, at = 0.3, boxwex = 0.4, col = "red2", xaxt = "n", yaxt = "n"))

legend("topright", c("yes", "maybe", "no"), fill = c("black", "grey50", "grey80"), bty = "n", cex = 1.2)
dev.off()

# percent correct if yes
median(people.df$fracCorrSY[!is.na(people.df$Experienced)])
median(people.df$fracCorrSY[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrSY[people.df$Experienced == "Student"], na.rm = TRUE)

# percent correct if maybe
median(people.df$fracCorrSM[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorrSM[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrSM[people.df$Experienced == "Student"], na.rm = TRUE)

# percent correct if no
median(people.df$fracCorrSN[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorrSN[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrSN[people.df$Experienced == "Student"], na.rm = TRUE)


## accuracy based on species and specimen level confidence
# percentage accuracy by size (based on mean diameter)
people.df$numcYy <- NA
people.df$numcYm <- NA
people.df$numcYn <- NA
people.df$numcMy <- NA
people.df$numcMm <- NA
people.df$numcMn <- NA
people.df$numcNy <- NA
people.df$numcNm <- NA
people.df$numcNn <- NA
people.df$fracCorrYy <- NA
people.df$fracCorrYm <- NA
people.df$fracCorrYn <- NA
people.df$fracCorrMy <- NA
people.df$fracCorrMm <- NA
people.df$fracCorrMn <- NA
people.df$fracCorrNy <- NA
people.df$fracCorrNm <- NA
people.df$fracCorrNn <- NA

sum(completeIDs[, people.col[2+5]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[5+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y" & completeIDs[, p.conf.col[5]] == "y", na.rm = TRUE) / sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[5+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y" & completeIDs[, people.col[2+5]] != "lost" & completeIDs[, p.conf.col[5]] == "y", na.rm = TRUE)

sum(completeIDs[, people.col[2+5]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[5+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "n" & completeIDs[, p.conf.col[5]] == "y", na.rm = TRUE) / sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[5+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "n" & completeIDs[, people.col[2+5]] != "lost" & completeIDs[, p.conf.col[5]] == "y", na.rm = TRUE)

for (i in 1:(length(people.col) - 2)){
  people.df$numcYy[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y" & completeIDs[, people.col[2+i]] != "lost" & completeIDs[, p.conf.col[i]] == "y", na.rm = TRUE)
  people.df$fracCorrYy[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y" & completeIDs[, p.conf.col[i]] == "y", na.rm = TRUE) / people.df$numcYy[2 + i]
  people.df$numcYm[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y" & completeIDs[, people.col[2+i]] != "lost" & completeIDs[, p.conf.col[i]] == "m", na.rm = TRUE)
  people.df$fracCorrYm[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y" & completeIDs[, p.conf.col[i]] == "m", na.rm = TRUE) / people.df$numcYm[2 + i]
  people.df$numcYn[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y" & completeIDs[, people.col[2+i]] != "lost" & completeIDs[, p.conf.col[i]] == "n", na.rm = TRUE)
  people.df$fracCorrYn[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "y" & completeIDs[, p.conf.col[i]] == "n", na.rm = TRUE) / people.df$numcYn[2 + i]
  people.df$numcMy[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "m" & completeIDs[, people.col[2+i]] != "lost" & completeIDs[, p.conf.col[i]] == "y", na.rm = TRUE)
  people.df$fracCorrMy[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "m" & completeIDs[, p.conf.col[i]] == "y", na.rm = TRUE) / people.df$numcMy[2 + i] 
  people.df$numcMm[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "m" & completeIDs[, people.col[2+i]] != "lost" & completeIDs[, p.conf.col[i]] == "m", na.rm = TRUE)
  people.df$fracCorrMm[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "m" & completeIDs[, p.conf.col[i]] == "m", na.rm = TRUE) / people.df$numcMm[2 + i]
  people.df$numcMn[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "m" & completeIDs[, people.col[2+i]] != "lost" & completeIDs[, p.conf.col[i]] == "n", na.rm = TRUE)
  people.df$fracCorrMn[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "m" & completeIDs[, p.conf.col[i]] == "n", na.rm = TRUE) / people.df$numcMn[2 + i]
  people.df$numcNy[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "n" & completeIDs[, people.col[2+i]] != "lost" & completeIDs[, p.conf.col[i]] == "y", na.rm = TRUE)
  people.df$fracCorrNy[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "n" & completeIDs[, p.conf.col[i]] == "y", na.rm = TRUE) / people.df$numcNy[2 + i]
  people.df$numcNm[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "n" & completeIDs[, people.col[2+i]] != "lost" & completeIDs[, p.conf.col[i]] == "m", na.rm = TRUE)
  people.df$fracCorrNm[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "n" & completeIDs[, p.conf.col[i]] == "m", na.rm = TRUE) / people.df$numcNm[2 + i] 
  people.df$numcNn[2 + i] <- sum(sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "n" & completeIDs[, people.col[2+i]] != "lost" & completeIDs[, p.conf.col[i]] == "n", na.rm = TRUE)
  people.df$fracCorrNn[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sapply(1:100, function (x) ifelse(completeIDs$DefinitiveID[x] %in% species, as.character(personData[personData$Person == people[i+2], as.character(completeIDs$DefinitiveID[x])]), NA)) == "n" & completeIDs[, p.conf.col[i]] == "n", na.rm = TRUE) / people.df$numcNn[2 + i] 
}

png("Figures/Percent accuracy species level.png")
plot(NULL, ylim = c(0, 1), xlim = c(0, 4.8), xaxt = "n", cex.axis = 1.5, xlab = "", ylab = "", las = 2)
with(people.df, boxplot(fracCorrSN ~ Experienced, boxwex = 0.2, col = c("steelblue3", "plum2"), at = c(2.4, 3.9), add = TRUE, xaxt = "n", yaxt = "n", xlim = c(0, 4.8)))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSN, add = TRUE, at = 0.9, boxwex = 0.4, col = "indianred2", xaxt = "n", yaxt = "n"))
axis(1, c(0.6,2.1,3.6), c("Full", levels(people.df$Experienced)), cex.axis = 1.5)

with(people.df, boxplot(fracCorrSM ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("royalblue3", "darkorchid2"), cex.axis = 1.5, at = c(2.1, 3.6), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSM, add = TRUE, at = 0.6, boxwex = 0.4, col = "firebrick2", xaxt = "n", yaxt = "n"))

with(people.df, boxplot(fracCorrSY ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("blue3", "purple2"), cex.axis = 1.5, at = c(1.8, 3.3), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSY, add = TRUE, at = 0.3, boxwex = 0.4, col = "red2", xaxt = "n", yaxt = "n"))

legend("topright", c("yes", "maybe", "no"), fill = c("black", "grey50", "grey80"), bty = "n", cex = 1.2)
dev.off()

# percent correct if Yes yes
median(people.df$fracCorrYy[!is.na(people.df$Experienced)])
median(people.df$fracCorrYy[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrYy[people.df$Experienced == "Student"], na.rm = TRUE)
# percent correct if Yes maybe
median(people.df$fracCorrYm[!is.na(people.df$Experienced)])
median(people.df$fracCorrYm[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrYm[people.df$Experienced == "Student"], na.rm = TRUE)
# percent correct if Yes no
median(people.df$fracCorrYn[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorrYn[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrYn[people.df$Experienced == "Student"], na.rm = TRUE)

# percent correct if Maybe yes
median(people.df$fracCorrMy[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorrMy[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrMy[people.df$Experienced == "Student"], na.rm = TRUE)
# percent correct if Maybe maybe
median(people.df$fracCorrMm[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorrMm[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrMm[people.df$Experienced == "Student"], na.rm = TRUE)
# percent correct if Maybe no
median(people.df$fracCorrMn[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorrMn[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrMn[people.df$Experienced == "Student"], na.rm = TRUE)

# percent correct if No yes
median(people.df$fracCorrNy[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorrNy[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrNy[people.df$Experienced == "Student"], na.rm = TRUE)
# percent correct if No maybe
median(people.df$fracCorrNm[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorrNm[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrNm[people.df$Experienced == "Student"], na.rm = TRUE)
# percent correct if No no
median(people.df$fracCorrNn[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorrNn[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorrNn[people.df$Experienced == "Student"], na.rm = TRUE)

# how many people does this apply to
sum(!is.na(people.df$fracCorrYy))
sum(!is.na(people.df$fracCorrYm))
sum(!is.na(people.df$fracCorrYn))
sum(!is.na(people.df$fracCorrMy))
sum(!is.na(people.df$fracCorrMm))
sum(!is.na(people.df$fracCorrMn))
sum(!is.na(people.df$fracCorrNy))
sum(!is.na(people.df$fracCorrNm))
sum(!is.na(people.df$fracCorrNn))
sum(!is.na(people.df$fracCorrYy[people.df$Experienced == "Experienced"]))
sum(!is.na(people.df$fracCorrYm[people.df$Experienced == "Experienced"]))
sum(!is.na(people.df$fracCorrYn[people.df$Experienced == "Experienced"]))
sum(!is.na(people.df$fracCorrMy[people.df$Experienced == "Experienced"]))
sum(!is.na(people.df$fracCorrMm[people.df$Experienced == "Experienced"]))
sum(!is.na(people.df$fracCorrMn[people.df$Experienced == "Experienced"]))
sum(!is.na(people.df$fracCorrNy[people.df$Experienced == "Experienced"]))
sum(!is.na(people.df$fracCorrNm[people.df$Experienced == "Experienced"]))
sum(!is.na(people.df$fracCorrNn[people.df$Experienced == "Experienced"]))
sum(!is.na(people.df$fracCorrYy[people.df$Experienced == "Student"]))
sum(!is.na(people.df$fracCorrYm[people.df$Experienced == "Student"]))
sum(!is.na(people.df$fracCorrYn[people.df$Experienced == "Student"]))
sum(!is.na(people.df$fracCorrMy[people.df$Experienced == "Student"]))
sum(!is.na(people.df$fracCorrMm[people.df$Experienced == "Student"]))
sum(!is.na(people.df$fracCorrMn[people.df$Experienced == "Student"]))
sum(!is.na(people.df$fracCorrNy[people.df$Experienced == "Student"]))
sum(!is.na(people.df$fracCorrNm[people.df$Experienced == "Student"]))
sum(!is.na(people.df$fracCorrNn[people.df$Experienced == "Student"]))

# number of Yes yes
median(people.df$numcYy[!is.na(people.df$Experienced) & people.df$numcYy > 0])
median(people.df$numcYy[people.df$Experienced == "Experienced" & people.df$numcYy > 0], na.rm = TRUE)
median(people.df$numcYy[people.df$Experienced == "Student" & people.df$numcYy > 0], na.rm = TRUE)
# number of Yes maybe
median(people.df$numcYm[!is.na(people.df$Experienced) & people.df$numcYm > 0])
median(people.df$numcYm[people.df$Experienced == "Experienced" & people.df$numcYm > 0], na.rm = TRUE)
median(people.df$numcYm[people.df$Experienced == "Student" & people.df$numcYm > 0], na.rm = TRUE)
# number of Yes no
median(people.df$numcYn[!is.na(people.df$Experienced) & people.df$numcYn > 0], na.rm = TRUE)
median(people.df$numcYn[people.df$Experienced == "Experienced" & people.df$numcYn > 0], na.rm = TRUE)
median(people.df$numcYn[people.df$Experienced == "Student" & people.df$numcYn > 0], na.rm = TRUE)

# number of Maybe yes
median(people.df$numcMy[!is.na(people.df$Experienced) & people.df$numcMy > 0], na.rm = TRUE)
median(people.df$numcMy[people.df$Experienced == "Experienced" & people.df$numcMy > 0], na.rm = TRUE)
median(people.df$numcMy[people.df$Experienced == "Student" & people.df$numcMy > 0], na.rm = TRUE)
# number of Maybe maybe
median(people.df$numcMm[!is.na(people.df$Experienced) & people.df$numcMm > 0], na.rm = TRUE)
median(people.df$numcMm[people.df$Experienced == "Experienced" & people.df$numcMm > 0], na.rm = TRUE)
median(people.df$numcMm[people.df$Experienced == "Student" & people.df$numcMm > 0], na.rm = TRUE)
# number of Maybe no
median(people.df$numcMn[!is.na(people.df$Experienced) & people.df$numcMn > 0], na.rm = TRUE)
median(people.df$numcMn[people.df$Experienced == "Experienced" & people.df$numcMn > 0], na.rm = TRUE)
median(people.df$numcMn[people.df$Experienced == "Student" & people.df$numcMn > 0], na.rm = TRUE)

# number of No yes
median(people.df$numcNy[!is.na(people.df$Experienced) & people.df$numcNy > 0], na.rm = TRUE)
median(people.df$numcNy[people.df$Experienced == "Experienced" & people.df$numcNy > 0], na.rm = TRUE)
median(people.df$numcNy[people.df$Experienced == "Student" & people.df$numcNy > 0], na.rm = TRUE)
# number of No maybe
median(people.df$numcNm[!is.na(people.df$Experienced) & people.df$numcNm > 0], na.rm = TRUE)
median(people.df$numcNm[people.df$Experienced == "Experienced" & people.df$numcNm > 0], na.rm = TRUE)
median(people.df$numcNm[people.df$Experienced == "Student" & people.df$numcNm > 0], na.rm = TRUE)
# number of No no
median(people.df$numcNn[!is.na(people.df$Experienced) & people.df$numcNn > 0], na.rm = TRUE)
median(people.df$numcNn[people.df$Experienced == "Experienced" & people.df$numcNn > 0], na.rm = TRUE)
median(people.df$numcNn[people.df$Experienced == "Student" & people.df$numcNn > 0], na.rm = TRUE)


## percentage accuracy by size (based on mean diameter)
# calculating percentage correct
people.df$fracCorr400 <- NA
people.df$fracCorr200 <- NA
people.df$fracCorr125 <- NA

sum(completeIDs[, people.col[2+5]] == completeIDs$DefinitiveID & sp.size$MeanDia < 200, na.rm = TRUE) / sum(sp.size$MeanDia < 200 & completeIDs[, people.col[2+5]] != "lost", na.rm = TRUE)

for (i in 1:length(p.conf.col)){
  people.df$fracCorr125[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sp.size$MeanDia < 200, na.rm = TRUE) / sum(sp.size$MeanDia < 200 & completeIDs[, people.col[2+i]] != "lost", na.rm = TRUE)
  people.df$fracCorr200[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sp.size$MeanDia < 400 & sp.size$MeanDia > 200, na.rm = TRUE) / sum(sp.size$MeanDia < 400 & sp.size$MeanDia > 200 & completeIDs[, people.col[2+i]] != "lost", na.rm = TRUE)
  people.df$fracCorr400[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & sp.size$MeanDia > 400, na.rm = TRUE) / sum(sp.size$MeanDia > 400 & completeIDs[, people.col[2+i]] != "lost", na.rm = TRUE)
}

png("Figures/Percent accuracy size.png")
plot(NULL, bty = "l", ylim = c(0, 1), xlim = c(0, 4.8), xaxt = "n", cex.axis = 1.5, xlab = "", ylab = "", las = 2)
with(people.df, boxplot(fracCorr125 ~ Experienced, boxwex = 0.2, col = c("steelblue3", "plum2"), at = c(1.8, 3.3), add = TRUE, xaxt = "n", yaxt = "n", xlim = c(0, 4.8)))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorr125, add = TRUE, at = 0.3, boxwex = 0.4, col = "indianred2", xaxt = "n", yaxt = "n"))
axis(1, c(0.6,2.1,3.6), c("Full", levels(people.df$Experienced)), cex.axis = 1.5)

with(people.df, boxplot(fracCorr200 ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("royalblue3", "darkorchid2"), cex.axis = 1.5, at = c(2.1, 3.6), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorr200, add = TRUE, at = 0.6, boxwex = 0.4, col = "firebrick2", xaxt = "n", yaxt = "n"))

with(people.df, boxplot(fracCorr400 ~ Experienced, xaxt = "n", boxwex = 0.2, col = c("blue3", "purple2"), cex.axis = 1.5, at = c(2.4, 3.9), add = TRUE, yaxt = "n"))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorr400, add = TRUE, at = 0.9, boxwex = 0.4, col = "red2", xaxt = "n", yaxt = "n"))

legend(3.4, 0.16, c(expression(paste(">400 ", mu, "m")), expression(paste("200-400 ", mu, "m")), expression(paste("125-200 ", mu, "m"))), fill = c("black", "grey50", "grey80"), bty = "n", cex = 1.2)
dev.off()

# percent correct if 125-200
median(people.df$fracCorr125[!is.na(people.df$Experienced)])
median(people.df$fracCorr125[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorr125[people.df$Experienced == "Student"], na.rm = TRUE)

# percent correct if 200-400
median(people.df$fracCorr200[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorr200[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorr200[people.df$Experienced == "Student"], na.rm = TRUE)

# percent correct if >400
median(people.df$fracCorr400[!is.na(people.df$Experienced)], na.rm = TRUE)
median(people.df$fracCorr400[people.df$Experienced == "Experienced"], na.rm = TRUE)
median(people.df$fracCorr400[people.df$Experienced == "Student"], na.rm = TRUE)

# 3c. Plotting confidence -------------------------------------------------
plot(c(1,3), c(0, 1), type = "n", xaxt = "n", xlab = "", ylab = "", las = 2, xlim = c(0.5, 3.5), main = "Specimen confidence")
axis(1, 1:3, labels = c("no", "maybe", "yes"))
for (i in 3:nrow(people.df)) {
  tmp.y <- people.df[i, c("fracCorrNwoUID", "fracCorrM", "fracCorrY")]
  points((1:3)[!is.na(tmp.y)], tmp.y[!is.na(tmp.y)], type = "b", col = rep(1:5, 10)[i])
}

plot(c(1,3), c(0, 1), type = "n", xaxt = "n", xlab = "", ylab = "", las = 2, xlim = c(0.5, 3.5), main = "Species confidence")
axis(1, 1:3, labels = c("no", "maybe", "yes"))
for (i in 1:nrow(people.df)) {
  tmp.y <- people.df[i, c("fracCorrSN", "fracCorrSM", "fracCorrSY")]
  points((1:3)[!is.na(tmp.y)], tmp.y[!is.na(tmp.y)], type = "b", col = rep(1:5, 10)[i])
}


# 4. Percentage accuracy per specimen ------------------------------------
sum(completeIDs[2, people.col[2:length(people.col)]] == as.character(completeIDs$DefinitiveID[2]) & completeIDs[2, people.col[2:length(people.col)]] != "lost") / sum(completeIDs[2, people.col[2:length(people.col)]] != "lost")

completeIDs$fracCorr <- apply(completeIDs[, people.col[2:length(people.col)]], 1, function(x) sum(x == x[1] & x != "lost") / sum(x != "lost"))

hist(completeIDs$fracCorr)


par(mfrow = c(2,1))
plot(factor(completeIDs$DefinitiveID), completeIDs$fracCorr, las = 2)
par(mfrow = c(1,1))

# 5. Dendrogram ---------------------------------------------------------
# n.b. use daisy rather than hist, as the data are factors
png("Figures/dendrogram.png")
plot(hclust(daisy(data.frame(t(completeIDs[, people.col])))))
dev.off()


# 6. Reshaping the dataset for modelling -----------------------------------

# 6a. Creating a long form dataframe --------------------------------------
head(completeIDs[,!(colnames(completeIDs) %in% p.conf.col) & colnames(completeIDs) != "fracCorr"])
IDs.long <- reshape(completeIDs[,!(colnames(completeIDs) %in% p.conf.col) & colnames(completeIDs) != "fracCorr"], varying = list(people.col[-1]), direction = "long", idvar = c("SpecNumber", "DefinitiveID"), times = people.col[-1], timevar = "Person")
IDs.long$Conf <- factor(c(rep(NA, 100), as.character(unlist(completeIDs[, p.conf.col]))))
names(IDs.long)[names(IDs.long) == "OriginalID"] <- "ID"
rownames(IDs.long) <- 1:nrow(IDs.long)
IDs.long$Person <- gsub("ID", "", IDs.long$Person)
head(IDs.long)
tail(IDs.long)

# correct ID?
IDs.long$Corr <- as.numeric(IDs.long$DefinitiveID == IDs.long$ID)

# 6b. Adding in the specimen level size data ----------------------------------
head(sp.size)
sp.size$logMeanDia <- log(sp.size$MeanDia)
size.col <- c("SpecNum", "MeanDia", "logMeanDia")
IDs.long <- merge(IDs.long, sp.size[, size.col], by.x = "SpecNumber", by.y = "SpecNum")
rm(size.col)

head(IDs.long)
tail(IDs.long)

# linear model of size against accuracy
tmp <- glm(Corr~ logMeanDia, data = IDs.long[!is.na(IDs.long$Corr),], family = binomial)
tmp2 <- predict(tmp, newdata = data.frame(logMeanDia = seq(5, 7, length.out = 100)), type = "response")
summary(tmp2)
plot(IDs.long$Corr~ IDs.long$logMeanDia)
lines(log(seq(180, 1050, length.out = 100)), tmp2)
rm(tmp, tmp2)

# 6c. Adding in the species level data -------------------------------
head(SpeciesData)
IDs.long <- merge(IDs.long, SpeciesData[, 2:ncol(SpeciesData)], by.x = "DefinitiveID", by.y = "Species")
head(IDs.long)

# 6d. Adding in the person level data -------------------------------------
# convert how long into an ordered factor for modelling based on quantiles
summary(personData$"How long have you worked with forams (years)")
quantile(personData$"How long have you worked with forams (years)", c(0.25, 0.5, 0.75))

personData$HowLong <- 0
personData$HowLong[personData$"How long have you worked with forams (years)" >= quantile(personData$"How long have you worked with forams (years)", c(0.25)) & personData$"How long have you worked with forams (years)" < quantile(personData$"How long have you worked with forams (years)", c(0.5))] <- 1
personData$HowLong[personData$"How long have you worked with forams (years)" >= quantile(personData$"How long have you worked with forams (years)", c(0.5)) & personData$"How long have you worked with forams (years)" < quantile(personData$"How long have you worked with forams (years)", c(0.75))] <- 2
personData$HowLong[personData$"How long have you worked with forams (years)" >= quantile(personData$"How long have you worked with forams (years)", c(0.75))] <- 3
cbind(personData$"How long have you worked with forams (years)", personData$HowLong)

# adding in person level data
head(personData)
IDs.long <- merge(IDs.long, personData)
head(IDs.long)

# adding in the species level confidence by person
IDs.long$SpConf <- apply(IDs.long, 1, function (x) as.character(x[x["DefinitiveID"]]))
IDs.long$SpConf <- factor(IDs.long$SpConf)


# 7. Modelling ------------------------------------------------------------

# 7a. Create a dataframe for modelling ------------------------------------
# remove the genus level identifications from the analysis dataframe
IDs.long.mod <- IDs.long[-grep("g", IDs.long$Conf),c("Corr", "HowLong", "logMeanDia", "Taught", "ED", "Conf", "SpConf", "Experience", "Gender", "DefinitiveID", "Person", "SpecNumber")]

# remove the juvenile
IDs.long.mod <- IDs.long.mod[-which(IDs.long.mod$SpecNumber == 53), ]
# remove the NAs
IDs.long.mod <- na.omit(IDs.long.mod)
IDs.long.mod$HowLong <- ordered(IDs.long.mod$HowLong)
IDs.long.mod$Taught <- factor(IDs.long.mod$Taught, levels = c("n", "m", "y"), ordered = TRUE)
IDs.long.mod$Conf <- factor(IDs.long.mod$Conf, levels = c("n", "m", "y"), ordered = TRUE)
IDs.long.mod$SpConf <- factor(IDs.long.mod$SpConf, levels = c("n", "m", "y"), ordered = TRUE)
IDs.long.mod$Experience <- ordered(IDs.long.mod$Experience)
IDs.long.mod$Person <- factor(IDs.long.mod$Person)
IDs.long.mod$SpecNumber <- factor(IDs.long.mod$SpecNumber)
IDs.long.mod$csLogMeanDia <- (IDs.long.mod$logMeanDia - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia)
IDs.long.mod$csED <- (IDs.long.mod$ED - mean(IDs.long.mod$ED)) / sd(IDs.long.mod$ED)
summary(IDs.long.mod)

# checking for correlations
summary(IDs.long.mod[, c("Corr", "HowLong", "logMeanDia", "Taught", "ED", "Conf", "SpConf", "Experience", "Gender")])
str(IDs.long.mod[, c("Corr", "HowLong", "logMeanDia", "Taught", "ED", "Conf", "SpConf", "Experience", "Gender")])
pairs(IDs.long.mod[, c("HowLong", "logMeanDia", "Taught", "ED", "Conf", "SpConf", "Experience", "Gender")])


# 7b. Fixed model ---------------------------------------------------------
# run a fixed model first with everything, to get an approximate idea
m1.fix <- glm(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender + DefinitiveID) + Person + SpecNumber, family = binomial, data = IDs.long.mod) # you get convergence warnings on this (I think because they aren't fitted as randoms.)
par(mfrow = c(2,2))
plot(m1.fix)
par(mfrow = c(1,1))
summary(m1.fix)

# use stepped to see what drops out immediately
stepped <- step(m1.fix)
summary(stepped)

sjp.glm(stepped, type = "pred", vars = c("csLogMeanDia"), geom.colors = matlab.like(30))


# 7c. Optimal model -------------------------------------------------------
# optimal random structure
glmcon <- glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=3e5)) # stop convergence issues
m1.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|DefinitiveID) + (1|Person) + (1|SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

sjp.glm(m1.me, type = "pred", vars = c("csLogMeanDia", "DefinitiveID"), geom.colors = matlab.like(30))

m2.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|DefinitiveID) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

m3.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|DefinitiveID) + (1|SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

m4.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|DefinitiveID), family = binomial, data = IDs.long.mod, control = glmcon)

m5.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|DefinitiveID) + (1|Person) + (1|SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

m6.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|DefinitiveID) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

m7.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|DefinitiveID) + (1|SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

m8.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|DefinitiveID), family = binomial, data = IDs.long.mod, control = glmcon)

m9.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|Person) + (1|SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

m10.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

m11.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

anova(m1.me, m2.me, m3.me, m4.me, m5.me, m6.me, m7.me, m8.me, m9.me, m10.me, m11.me)

# best model is m5.me (i.e. no random slopes only random intercepts)
summary(m5.me)

re.m5.me <- residuals(m5.me, type = "response")
fi.m5.me <- fitted(m5.me)
plot(fi.m5.me, re.m5.me)

# run model simplification
m5a.me <- update(m5.me, .~. -csLogMeanDia:Gender)
anova(m5.me, m5a.me)
summary(m5a.me)
m5b.me <- update(m5a.me, .~. -csLogMeanDia:SpConf)
anova(m5b.me, m5a.me)
summary(m5b.me)
m5c.me <- update(m5b.me, .~. -csLogMeanDia:csED)
anova(m5b.me, m5c.me)
summary(m5c.me)
m5d.me <- update(m5c.me, .~. -csED)
anova(m5d.me, m5c.me)
summary(m5d.me)
m5e.me <- update(m5d.me, .~. -csLogMeanDia:Conf)
anova(m5d.me, m5e.me)
summary(m5e.me)
m5f.me <- update(m5e.me, .~. -Gender)
anova(m5f.me, m5e.me)
summary(m5f.me)

anova(m5f.me)
Anova(m5f.me)
r.squaredGLMM(m5f.me)

write.csv(Anova(m5f.me), file = "Outputs/anova_m5f_me.csv")


# 7d. marginal effects of the different EVs -------------------------------
summary(m5f.me)
r.squaredGLMM(m5f.me) # r2m = 0.4366040

# drop each EV in turn to get the marginal effect
# log size
m5f.me.cLMD <- update(m5f.me, ~. - csLogMeanDia - csLogMeanDia:HowLong - csLogMeanDia:Taught - csLogMeanDia:Experience)
summary(m5f.me.cLMD)
r.squaredGLMM(m5f.me.cLMD) # r2m = 0.4022265
0.4366040 - 0.4022265
# how long
m5f.me.HL <- update(m5f.me, ~. - HowLong - csLogMeanDia:HowLong)
summary(m5f.me.HL)
r.squaredGLMM(m5f.me.HL) # r2m = 0.3912333
0.4366040 - 0.3912333
# taught
m5f.me.T <- update(m5f.me, ~. - Taught - csLogMeanDia:Taught)
summary(m5f.me.T)
r.squaredGLMM(m5f.me.T) # r2m = 0.2620388 
0.4366040 - 0.2620388
# Experience
m5f.me.E <- update(m5f.me, ~. - Experience - csLogMeanDia:Experience)
summary(m5f.me.E)
r.squaredGLMM(m5f.me.E) # r2m = 0.4193970
0.4366040 - 0.4193970
# Specimen confidence
m5f.me.C <- update(m5f.me, ~. - Conf)
summary(m5f.me.C)
r.squaredGLMM(m5f.me.C) # r2m = 0.3822157
0.4366040 - 0.3822157
# Species confidence
m5f.me.SC <- update(m5f.me, ~. - SpConf)
summary(m5f.me.SC)
r.squaredGLMM(m5f.me.SC) # r2m = 0.4299028
0.4366040 - 0.4299028

# 7e. Model plotting ------------------------------------------------------
# can't simplify further, so move to plotting
pdf("Figures/GLM_effects.pdf")
sjp.glmer(m5f.me)
sjp.glmer(m5f.me, type = "fe")
sjp.glmer(m5f.me, type = "eff")
sjp.int(m5f.me, type = "eff")

p <- sjp.glmer(m5f.me, type = "eff", facet.grid = FALSE, 
             show.ci = TRUE, prnt.plot = FALSE)$plot.list
# plot all marginal effects, as grid, proper x-axes
# also, set margins for this example
plot_grid(p, margin = c(0.3, 0.3, 0.3, 0.3))
sjp.glmer(m5f.me, type = "pred", vars = "csLogMeanDia")
sjp.glmer(m5f.me, type = "pred", vars = "HowLong")
sjp.glmer(m5f.me, type = "pred", vars = "Taught")
sjp.glmer(m5f.me, type = "pred", vars = "Conf")
sjp.glmer(m5f.me, type = "pred", vars = "SpConf")
sjp.glmer(m5f.me, type = "pred", vars = "Experience")
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "HowLong"), title = "How Long")
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Taught"), title = "Taught")
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Conf"), title = "Specimen Confidence")
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "SpConf"), title = "Species confidence")
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Experience"), title = "Experience")
dev.off()

png("Figures/GLM_size_sp.png", 1500, 1200, res = 120)
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "DefinitiveID"), geom.colors = rep("blue4", 26), show.ci = TRUE, cex = 10, cex.axis = 10, cex.labels = 10)
dev.off()

for (i in 1:26){
  png(paste("Figures/GLM_size_", i, ".png", sep = ""), 900, 500)
  sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "DefinitiveID"), facet.grid = FALSE, pch = ".", geom.colors = c(rep("grey90",i-1), "black", rep("grey90", 27-i)))
  dev.off()
}

# same thing but saved into one pdf
pdf("Figures/GLM_size.pdf", 10,7)
for (i in 1:26){
  sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "DefinitiveID"), facet.grid = FALSE, pch = ".", geom.colors = c(rep("grey90",i-1), "black", rep("grey90", 27-i)))
}
dev.off()

# how do the sjp.glmer plots work?
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Conf"), title = "Conf", show.ci = TRUE, facet.grid = FALSE)

aug <- augment(m5f.me, IDs.long.mod)
aug$fit.rs <- 1/(1 + exp(-aug$.fitted)) # back transformation(inv.logit)
aug$fix.rs <- 1/(1 + exp(-aug$.fixed))

inv.logit <- function(x) 1/(1 + exp(x))

# generates a linear model of the fitted data, so not a good match to the sjp plots
p4 <- ggplot(aug, aes(csLogMeanDia, fit.rs, color=Conf, fill=Conf)) + geom_smooth(method = lm) + ylab("Corr") 
plot(p4)

# how the sjp.glmer function does it:
tmp <- merTools::predictInterval(m5f.me, newdata = IDs.long.mod, which = "full", type = "probability", level = 0.95) # n.b. this is an approximation, so answers will be subtly different each time

aug$pi.fit <- tmp$fit
aug$pi.upr <- tmp$upr
aug$pi.lwr <- tmp$lwr

mp <- ggplot(aug, aes(x = csLogMeanDia, y = pi.fit, colour = Conf, fill = Conf)) + stat_smooth(method = "glm", method.args = list(family = m5f.me@resp$family), se = TRUE, fullrange = TRUE, level = 0.95, alpha = 0.3) + coord_cartesian(ylim = c(0,1)) 
plot(mp)

# so considering the plots for the paper
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Conf"), title = "Conf", show.ci = TRUE, facet.grid = FALSE)

plotLmerMeans(m5f.me, var.of.interest = "Conf", var.interaction = "csLogMeanDia")
tmp <- plotLmerMeans(m5f.me, var.of.interest = "Conf", var.interaction = "csLogMeanDia", transform_fun = inv.logit, ylims = c(0,1), meancol = 1:3, return_quietly = FALSE)

png("Figures/Size_Conf_sjp_glmer.png")
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Conf"), title = "Conf", show.ci = TRUE, facet.grid = FALSE)
dev.off()

png("Figures/Size_Conf_plotLmerMeans.png")
plotLmerMeans(m5f.me, var.of.interest = "Conf", var.interaction = "csLogMeanDia", transform_fun = inv.logit, ylims = c(0,1), meancol = 1:3)
legend("topright", legend = levels(IDs.long.mod$Conf), col = 1:3, lty = 1, lwd = 3)
dev.off()

png("Figures/Size_Conf_plotLmerMeansRepeat.png")
plotLmerMeansRepeat(m5f.me, var.of.interest = "Conf", var.interaction = "csLogMeanDia", transform_fun = inv.logit, ylims = c(0,1), meancol = 1:3, mainlabel = "Confidence")
legend("topright", legend = levels(IDs.long.mod$Conf), col = 1:3, lty = 1, lwd = 3)
dev.off()

png("Figures/Size_Taught_plotLmerMeansRepeat.png")
plotLmerMeansRepeat(m5f.me, var.of.interest = "Taught", var.interaction = "csLogMeanDia", transform_fun = inv.logit, ylims = c(0,1), meancol = 1:3, mainlabel = "Taught")
legend("topright", legend = levels(IDs.long.mod$Conf), col = 1:3, lty = 1, lwd = 3)
dev.off()

# using Adriana's code (but edited to highest level as baseline)
pdf("Figures/Terms_plotLmerMeansRepeat.pdf")
plotLmerMeansRepeat(m5f.me, var.of.interest = "csLogMeanDia", transform_fun = inv.logit, ylims = c(0,1), mainlabel = "Log Mean Diameter")

plotLmerMeansRepeat(m5f.me, var.of.interest = "Conf", transform_fun = inv.logit, ylims = c(0,1), mainlabel = "Specimen Confidence")

plotLmerMeansRepeat(m5f.me, var.of.interest = "SpConf", transform_fun = inv.logit, ylims = c(0,1), mainlabel = "Species Confidence")

plotLmerMeansRepeat(m5f.me, var.of.interest = "Experience", var.interaction = "csLogMeanDia", transform_fun = inv.logit, ylims = c(0,1), meancol = c("red3", "blue2", "green4"), mainlabel= "Experience", xlabel = "Log Mean Diameter")
legend("bottomright", legend = levels(IDs.long.mod$Experience), col = c("red3", "blue2", "green4"), lty = 1, lwd = 3)

plotLmerMeansRepeat(m5f.me, var.of.interest = "HowLong", var.interaction = "csLogMeanDia", transform_fun = inv.logit, ylims = c(0,1), meancol = c("red3", "blue2", "green4", "purple"), mainlabel= "How Long", xlabel = "Log Mean Diameter")
legend("bottomright", legend = levels(IDs.long.mod$HowLong), col = c("red3", "blue2", "green4", "purple"), lty = 1, lwd = 3)

plotLmerMeansRepeat(m5f.me, var.of.interest = "Taught", var.interaction = "csLogMeanDia", transform_fun = inv.logit, ylims = c(0,1), meancol = c("red3", "blue2", "green4"), mainlabel= "Taught", xlabel = "Log Mean Diameter")
legend("topright", legend = levels(IDs.long.mod$Conf), col = c("red3", "blue2", "green4"), lty = 1, lwd = 3)
dev.off()


# sjp.glmer code
pdf("Figures/Terms_SJP_Glmer.pdf")
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia"), title = "Log Mean Diameter", show.ci = TRUE, facet.grid = FALSE)

sjp.glmer(m5f.me, type = "pred", vars = c("Conf"), title = "Specimen Confidence", show.ci = TRUE, facet.grid = FALSE)

sjp.glmer(m5f.me, type = "pred", vars = c("SpConf"), title = "Species Confidence", show.ci = TRUE, facet.grid = FALSE)

sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Experience"), title = "Experience", show.ci = TRUE, facet.grid = FALSE)

sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "HowLong"), title = "How Long", show.ci = TRUE, facet.grid = FALSE)

sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Taught"), title = "Taught", show.ci = TRUE, facet.grid = FALSE)
dev.off()

# paper plots, we decided to just use sjp.glmer
# the scale for the x-axis is currently not very informative (centred / scaled log mean diameter). See if I can work out what this is
# (IDs.long.mod$logMeanDia - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia)
-1 * sd(IDs.long.mod$logMeanDia) + mean(IDs.long.mod$logMeanDia) # 5.297101
0 * sd(IDs.long.mod$logMeanDia) + mean(IDs.long.mod$logMeanDia) # 5.793243
1 * sd(IDs.long.mod$logMeanDia) + mean(IDs.long.mod$logMeanDia) # 6.289385
2 * sd(IDs.long.mod$logMeanDia) + mean(IDs.long.mod$logMeanDia) # 6.785527
# or unlogged:
exp(-1 * sd(IDs.long.mod$logMeanDia) + mean(IDs.long.mod$logMeanDia)) # 199.7568
exp(0 * sd(IDs.long.mod$logMeanDia) + mean(IDs.long.mod$logMeanDia)) # 328.0751
exp(1 * sd(IDs.long.mod$logMeanDia) + mean(IDs.long.mod$logMeanDia)) # 538.8216
exp(2 * sd(IDs.long.mod$logMeanDia) + mean(IDs.long.mod$logMeanDia)) # 884.9459
# so if I want them at 200, 400, 600, 800
(log(200) - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia) # -0.9975478
(log(400) - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia) # 0.3995266
(log(600) - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia) # 1.216763
(log(800) - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia) # 1.796601

png("Figures/Size_Taught_sjp_glmer.png")
sjp.setTheme(theme_classic(), axis.title.size = 1.4, axis.textsize = 1.2)
p <- sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Taught"), title = "", show.ci = TRUE, facet.grid = FALSE, axis.title = c(expression(paste("Mean Diameter / ", mu, "m")), "Percentage Correct"), geom.size = 1.2, point.alpha = 0.5, prnt.plot = FALSE)
p$plot + scale_x_continuous(breaks = c(-0.9975478, 0.3995266, 1.216763, 1.796601), labels = c(200, 400, 600, 800))
dev.off()

png("Figures/Size_Exp_sjp_glmer.png")
sjp.setTheme(theme_classic(), axis.title.size = 1.4, axis.textsize = 1.2, panel.gridcol = "white", legend.size = 1.2, legend.item.backcol = "white")
p1 <- sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Experience"), title = "", show.ci = TRUE, facet.grid = FALSE, axis.title = c(expression(paste("Mean Diameter / ", mu, "m")), "Percentage Correct"), geom.size = 1.2, point.alpha = 0.5, prnt.plot = FALSE)
p1$plot + scale_x_continuous(breaks = c(-0.9975478, 0.3995266, 1.216763, 1.796601), labels = c(200, 400, 600, 800))
dev.off()

png("Figures/Size_HowLong_sjp_glmer.png")
sjp.setTheme(theme_classic(), axis.title.size = 1.4, axis.textsize = 1.2, panel.gridcol = "white", legend.size = 1.2, legend.item.backcol = "white")
p2 <- sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "HowLong"), title = "", show.ci = TRUE, facet.grid = FALSE, axis.title = c(expression(paste("Mean Diameter / ", mu, "m")), "Percentage Correct"), geom.size = 1.2, point.alpha = 0.5, prnt.plot = FALSE)
p2$plot + scale_x_continuous(breaks = c(-0.9975478, 0.3995266, 1.216763, 1.796601), labels = c(200, 400, 600, 800))
dev.off()

# 8. Try running a model with unordered factors ----------------------------------------------
# unordered factors (sorted to make the base level the best value)
IDs.long.mod2 <- IDs.long.mod
IDs.long.mod2$HowLong <- factor(as.character(IDs.long.mod2$HowLong), levels = c("3", "2", "1", "0"))
IDs.long.mod2$Taught <- factor(as.character(IDs.long.mod2$Taught), levels = c("y", "m", "n"))
IDs.long.mod2$Conf <- factor(as.character(IDs.long.mod2$Conf), levels = c("y", "m", "n"))
IDs.long.mod2$SpConf <- factor(as.character(IDs.long.mod2$SpConf), levels = c("y", "m", "n"))
IDs.long.mod2$Experience <- factor(as.character(IDs.long.mod2$Experience), levels = c("2", "1", "0"))
summary(IDs.long.mod2)

m5f.me2 <- update(m5f.me, data = IDs.long.mod2)
summary(m5f.me2)

inv.logit(summary(m5f.me2)$coefficients[1, 1])
inv.logit(summary(m5f.me2)$coefficients[1, 1] + summary(m5f.me2)$coefficients[3:13, 1])

# simple linear model


# 9. Create a confusion matrix --------------------------------------------
# full confusion matrix
sp.idd <- read.csv("Data/Specieslist.csv")$DefinitiveID
conf.sp <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$ID != "lost"], levels = sp.idd), factor(IDs.long$ID[IDs.long$ID != "lost"], levels = sp.idd))
conf.splist <- sp.idd[rowSums(conf.sp$table) > 0]
par(mfrow = c(1,1))
png("Figures/confusion_full.png", 800, 600)
par(mar = c(15, 15, 2, 2))
image(t(conf.sp$table / rowSums(conf.sp$table)), col = matlab.like(1000), ylim = c((28*((1+1/43))/44) - (1+1/43)/44/2, -(1+1/43)/44/2 ), axes = FALSE)
axis(1, seq(0,1, length.out = length(sp.idd)), sp.idd, las = 2)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(sp.idd[rowSums(conf.sp$table) > 0])), sp.idd[rowSums(conf.sp$table) > 0], las = 1)
dev.off()

png("Figures/confusion_full_grey.png", 800, 600)
par(mar = c(15, 15, 2, 2))
image(t(conf.sp$table / rowSums(conf.sp$table)), col = c("grey70", matlab.like(1000)), ylim = c((28*((1+1/43))/44) - (1+1/43)/44/2, -(1+1/43)/44/2 ), axes = FALSE)
axis(1, seq(0,1, length.out = length(sp.idd)), sp.idd, las = 2)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(sp.idd[rowSums(conf.sp$table) > 0])), sp.idd[rowSums(conf.sp$table) > 0], las = 1)
par(par.def)
dev.off()

# full confusion matrix with key
png("Figures/confusion_full_grey_key.png", 1000, 700)
par(fig = c(0, 0.9, 0, 1))
par(mar = c(15, 15, 2, 2))
# confusion matrix
image(t(conf.sp$table / rowSums(conf.sp$table)), col = c("grey70", matlab.like(1000)), ylim = c((28*((1+1/43))/44) - (1+1/43)/44/2, -(1+1/43)/44/2 ), axes = FALSE)

title(xlab = expression(bold("Individual ID")), ylab = expression(bold("Definitive ID")), line = 13, cex.axis = 1.2)

# x axis names
xaxis.names <- sp.idd
axis(1, seq(0,1, length.out = length(xaxis.names))[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], xaxis.names[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], las = 2, font = 3, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[-grep("^[A-Z]", xaxis.names)], xaxis.names[-grep("^[A-Z]", xaxis.names)], las = 2, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[str_count(xaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", xaxis.names[str_count(xaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", xaxis.names[str_count(xaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Individual ID
axis(3, seq(0,1, length.out = length(xaxis.names)) - 0.002, table(IDs.long$ID)[match(xaxis.names, names(table(IDs.long$ID)))], cex.axis = 0.8, las = 2, tick = FALSE, mgp = c(3, 0.5, 0))

# y axis names
yaxis.names <- sp.idd[rowSums(conf.sp$table) > 0]
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)), labels = FALSE)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 1], yaxis.names[str_count(yaxis.names, " ") == 1], las = 1, font = 3, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[-grep("^[A-Z]", yaxis.names)], yaxis.names[-grep("^[A-Z]", yaxis.names)], las = 1, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", yaxis.names[str_count(yaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", yaxis.names[str_count(yaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Definitive ID
axis(4, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)) - 0.002, (table(IDs.long$DefinitiveID)/23)[match(yaxis.names, names(table(IDs.long$DefinitiveID)))], cex.axis = 0.8, las = 1, tick = FALSE, mgp = c(3, 0.5, 0))

# add key
par(fig = c(0.9, 1, 0.35, 1), new = TRUE)
par(mai = c(1, 0.4, 0.8, 0.5))
names.key <- rep("", 101)
names.key[seq(1, 101, by = 20)] <- seq(0.0, 1.0, by = 0.2)
key <- rep(1, 101)
barplot(key, main = "\nFraction of\nSpecimens", horiz = TRUE, space = 0, border = NA, col = c("grey70", matlab.like(100)), fg = "white", cex.main = 1.1, font.main = 1, xaxt = "n", yaxt = "n")
axis(4, at = 1:101, names.key, las = 1, mgp = c(0, 1, 0), xaxt = "n", cex.axis = 1, tick = FALSE)
par(par.def)
dev.off() 

conf.sp$overall


# confusion matrix for experienced workers
conf.exp <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Experienced == "Experienced"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Experienced == "Experienced"], levels = sp.idd))
conf.exp.frac <- t(conf.exp$table / rowSums(conf.exp$table))
conf.exp.frac[conf.exp.frac == "NaN"] <- 0 # otherwise white patches appear across the plots

png("Figures/confusion_exp_grey_key.png", 1000, 700)
par(fig = c(0, 0.9, 0, 1))
par(mar = c(15, 15, 2, 2))
# confusion matrix
image(conf.exp.frac, col = c("grey70", matlab.like(1000)), ylim = c((28*((1+1/43))/44) - (1+1/43)/44/2, -(1+1/43)/44/2 ), axes = FALSE, bg = "grey70")

title(xlab = expression(bold("Individual ID")), ylab = expression(bold("Definitive ID")), line = 13, cex.axis = 1.2)

# x axis names
xaxis.names <- sp.idd
axis(1, seq(0,1, length.out = length(xaxis.names))[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], xaxis.names[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], las = 2, font = 3, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[-grep("^[A-Z]", xaxis.names)], xaxis.names[-grep("^[A-Z]", xaxis.names)], las = 2, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[str_count(xaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", xaxis.names[str_count(xaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", xaxis.names[str_count(xaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Individual ID
axis(3, seq(0,1, length.out = length(xaxis.names)) - 0.002, table(IDs.long$ID[IDs.long$Experienced == "Experienced"])[match(xaxis.names, names(table(IDs.long$ID)))], cex.axis = 0.8, las = 2, tick = FALSE, mgp = c(3, 0.5, 0))

# y axis names
yaxis.names <- sp.idd[rowSums(conf.sp$table) > 0]
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)), labels = FALSE)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 1], yaxis.names[str_count(yaxis.names, " ") == 1], las = 1, font = 3, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[-grep("^[A-Z]", yaxis.names)], yaxis.names[-grep("^[A-Z]", yaxis.names)], las = 1, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", yaxis.names[str_count(yaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", yaxis.names[str_count(yaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Definitive ID
axis(4, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)) - 0.002, (table(IDs.long$DefinitiveID[IDs.long$Experienced == "Experienced"])/4)[match(yaxis.names, names(table(IDs.long$DefinitiveID)))], cex.axis = 0.8, las = 1, tick = FALSE, mgp = c(3, 0.5, 0))

# add key
par(fig = c(0.9, 1, 0.35, 1), new = TRUE)
par(mai = c(1, 0.4, 0.8, 0.5))
names.key <- rep("", 101)
names.key[seq(1, 101, by = 20)] <- seq(0.0, 1.0, by = 0.2)
key <- rep(1, 101)
barplot(key, main = "\nFraction of\nSpecimens", horiz = TRUE, space = 0, border = NA, col = c("grey70", matlab.like(100)), fg = "white", cex.main = 1.1, font.main = 1, xaxt = "n", yaxt = "n")
axis(4, at = 1:101, names.key, las = 1, mgp = c(0, 1, 0), xaxt = "n", cex.axis = 1, tick = FALSE)
par(par.def)
dev.off() 

conf.exp$overall


# confusion matrix for students
conf.stu <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Experienced == "Student"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Experienced == "Student"], levels = sp.idd))
conf.stu.frac <- t(conf.stu$table / rowSums(conf.stu$table))
conf.stu.frac[conf.stu.frac == "NaN"] <- 0 # otherwise white patches appear across the plots

png("Figures/confusion_stu_grey_key.png", 1000, 700)
par(fig = c(0, 0.9, 0, 1))
par(mar = c(15, 15, 2, 2))
# confusion matrix
image(conf.stu.frac, col = c("grey70", matlab.like(1000)), ylim = c((28*((1+1/43))/44) - (1+1/43)/44/2, -(1+1/43)/44/2 ), axes = FALSE, bg = "grey70")

title(xlab = expression(bold("Individual ID")), ylab = expression(bold("Definitive ID")), line = 13, cex.axis = 1.2)

# x axis names
xaxis.names <- sp.idd
axis(1, seq(0,1, length.out = length(xaxis.names))[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], xaxis.names[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], las = 2, font = 3, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[-grep("^[A-Z]", xaxis.names)], xaxis.names[-grep("^[A-Z]", xaxis.names)], las = 2, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[str_count(xaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", xaxis.names[str_count(xaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", xaxis.names[str_count(xaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Individual ID
axis(3, seq(0,1, length.out = length(xaxis.names)) - 0.002, table(IDs.long$ID[IDs.long$Experienced == "Student"])[match(xaxis.names, names(table(IDs.long$ID)))], cex.axis = 0.8, las = 2, tick = FALSE, mgp = c(3, 0.5, 0))

# y axis names
yaxis.names <- sp.idd[rowSums(conf.sp$table) > 0]
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)), labels = FALSE)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 1], yaxis.names[str_count(yaxis.names, " ") == 1], las = 1, font = 3, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[-grep("^[A-Z]", yaxis.names)], yaxis.names[-grep("^[A-Z]", yaxis.names)], las = 1, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", yaxis.names[str_count(yaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", yaxis.names[str_count(yaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Definitive ID
axis(4, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)) - 0.002, (table(IDs.long$DefinitiveID[IDs.long$Experienced == "Student"])/19)[match(yaxis.names, names(table(IDs.long$DefinitiveID)))], cex.axis = 0.8, las = 1, tick = FALSE, mgp = c(3, 0.5, 0))

# add key
par(fig = c(0.9, 1, 0.35, 1), new = TRUE)
par(mai = c(1, 0.4, 0.8, 0.5))
names.key <- rep("", 101)
names.key[seq(1, 101, by = 20)] <- seq(0.0, 1.0, by = 0.2)
key <- rep(1, 101)
barplot(key, main = "\nFraction of\nSpecimens", horiz = TRUE, space = 0, border = NA, col = c("grey70", matlab.like(100)), fg = "white", cex.main = 1.1, font.main = 1, xaxt = "n", yaxt = "n")
axis(4, at = 1:101, names.key, las = 1, mgp = c(0, 1, 0), xaxt = "n", cex.axis = 1, tick = FALSE)
par(par.def)
dev.off()
  
  
# confusion matrix for confident IDs
conf.conf <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Conf == "y"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Conf == "y"], levels = sp.idd))
conf.conf.frac <- t(conf.conf$table / rowSums(conf.conf$table))
conf.conf.frac[conf.conf.frac == "NaN"] <- 0 # otherwise white patches appear across the plots

png("Figures/confusion_conf_grey_key.png", 1000, 700)
par(fig = c(0, 0.9, 0, 1))
par(mar = c(15, 15, 2, 2))
# confusion matrix
image(conf.conf.frac, col = c("grey70", matlab.like(1000)), ylim = c((28*((1+1/43))/44) - (1+1/43)/44/2, -(1+1/43)/44/2 ), axes = FALSE, bg = "grey70")

title(xlab = expression(bold("Individual ID")), ylab = expression(bold("Definitive ID")), line = 13, cex.axis = 1.2)

# x axis names
xaxis.names <- sp.idd
axis(1, seq(0,1, length.out = length(xaxis.names))[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], xaxis.names[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], las = 2, font = 3, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[-grep("^[A-Z]", xaxis.names)], xaxis.names[-grep("^[A-Z]", xaxis.names)], las = 2, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[str_count(xaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", xaxis.names[str_count(xaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", xaxis.names[str_count(xaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Individual ID
axis(3, seq(0,1, length.out = length(xaxis.names)) - 0.002, table(IDs.long$ID[IDs.long$Conf == "y"])[match(xaxis.names, names(table(IDs.long$ID)))], cex.axis = 0.8, las = 2, tick = FALSE, mgp = c(3, 0.5, 0))

# y axis names
yaxis.names <- sp.idd[rowSums(conf.sp$table) > 0]
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)), labels = FALSE)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 1], yaxis.names[str_count(yaxis.names, " ") == 1], las = 1, font = 3, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[-grep("^[A-Z]", yaxis.names)], yaxis.names[-grep("^[A-Z]", yaxis.names)], las = 1, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", yaxis.names[str_count(yaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", yaxis.names[str_count(yaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Definitive ID
axis(4, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)) - 0.002, (table(IDs.long$DefinitiveID)/23)[match(yaxis.names, names(table(IDs.long$DefinitiveID)))], cex.axis = 0.8, las = 1, tick = FALSE, mgp = c(3, 0.5, 0))

# add key
par(fig = c(0.9, 1, 0.35, 1), new = TRUE)
par(mai = c(1, 0.4, 0.8, 0.5))
names.key <- rep("", 101)
names.key[seq(1, 101, by = 20)] <- seq(0.0, 1.0, by = 0.2)
key <- rep(1, 101)
barplot(key, main = "\nFraction of\nSpecimens", horiz = TRUE, space = 0, border = NA, col = c("grey70", matlab.like(100)), fg = "white", cex.main = 1.1, font.main = 1, xaxt = "n", yaxt = "n")
axis(4, at = 1:101, names.key, las = 1, mgp = c(0, 1, 0), xaxt = "n", cex.axis = 1, tick = FALSE)
par(par.def)
dev.off() 

conf.conf$overall


conf.unconf <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Conf == "n"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Conf == "n"], levels = sp.idd))
image(conf.unconf$table / rowSums(conf.unconf$table), axes = FALSE, col = matlab.like(100))
axis(1, seq(0,1, length.out = length(sp.idd)), sp.idd, las = 2)
axis(2, seq(0,1, length.out = length(sp.idd)), sp.idd, las = 1)
conf.unconf$overall

# confusion matrix for unconfident IDs
conf.unconf <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Conf == "n"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Conf == "n"], levels = sp.idd))
conf.unconf.frac <- t(conf.unconf$table / rowSums(conf.unconf$table))
conf.unconf.frac[conf.unconf.frac == "NaN"] <- 0 # otherwise white patches appear across the plots

png("Figures/confusion_unconf_grey_key.png", 1000, 700)
par(fig = c(0, 0.9, 0, 1))
par(mar = c(15, 15, 2, 2))
# confusion matrix
image(conf.unconf.frac, col = c("grey70", matlab.like(1000)), ylim = c((28*((1+1/43))/44) - (1+1/43)/44/2, -(1+1/43)/44/2 ), axes = FALSE, bg = "grey70")

title(xlab = expression(bold("Individual ID")), ylab = expression(bold("Definitive ID")), line = 13, cex.axis = 1.2)

# x axis names
xaxis.names <- sp.idd
axis(1, seq(0,1, length.out = length(xaxis.names))[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], xaxis.names[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], las = 2, font = 3, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[-grep("^[A-Z]", xaxis.names)], xaxis.names[-grep("^[A-Z]", xaxis.names)], las = 2, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[str_count(xaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", xaxis.names[str_count(xaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", xaxis.names[str_count(xaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Individual ID
axis(3, seq(0,1, length.out = length(xaxis.names)) - 0.002, table(IDs.long$ID[IDs.long$Conf == "n"])[match(xaxis.names, names(table(IDs.long$ID)))], cex.axis = 0.8, las = 2, tick = FALSE, mgp = c(3, 0.5, 0))

# y axis names
yaxis.names <- sp.idd[rowSums(conf.sp$table) > 0]
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)), labels = FALSE)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 1], yaxis.names[str_count(yaxis.names, " ") == 1], las = 1, font = 3, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[-grep("^[A-Z]", yaxis.names)], yaxis.names[-grep("^[A-Z]", yaxis.names)], las = 1, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", yaxis.names[str_count(yaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", yaxis.names[str_count(yaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Definitive ID
axis(4, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)) - 0.002, (table(IDs.long$DefinitiveID)/23)[match(yaxis.names, names(table(IDs.long$DefinitiveID)))], cex.axis = 0.8, las = 1, tick = FALSE, mgp = c(3, 0.5, 0))

# add key
par(fig = c(0.9, 1, 0.35, 1), new = TRUE)
par(mai = c(1, 0.4, 0.8, 0.5))
names.key <- rep("", 101)
names.key[seq(1, 101, by = 20)] <- seq(0.0, 1.0, by = 0.2)
key <- rep(1, 101)
barplot(key, main = "\nFraction of\nSpecimens", horiz = TRUE, space = 0, border = NA, col = c("grey70", matlab.like(100)), fg = "white", cex.main = 1.1, font.main = 1, xaxt = "n", yaxt = "n")
axis(4, at = 1:101, names.key, las = 1, mgp = c(0, 1, 0), xaxt = "n", cex.axis = 1, tick = FALSE)
par(par.def)
dev.off() 

# confusion matrix values
conf.sp$overall
conf.exp$overall
conf.stu$overall
conf.conf$overall
conf.unconf$overall

# Accuracy is the fraction of correct identifications
# Accuracy null is (I think) the expected accuracy - which is the accuracy that would be expected based on random.
# the kappa statistic is (observed - expected) / (1 - expected)



# confusion matrix for large specimens 
conf.large <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$MeanDia >= 200], levels = sp.idd), factor(IDs.long$ID[IDs.long$MeanDia >= 200], levels = sp.idd))
conf.large.frac <- t(conf.large$table / rowSums(conf.large$table))
conf.large.frac[conf.large.frac == "NaN"] <- 0 # otherwise white patches appear across the plots

png("Figures/confusion_large_grey_key.png", 1000, 700)
par(fig = c(0, 0.9, 0, 1))
par(mar = c(15, 15, 2, 2))
# confusion matrix
image(conf.large.frac, col = c("grey70", matlab.like(1000)), ylim = c((28*((1+1/43))/44) - (1+1/43)/44/2, -(1+1/43)/44/2 ), axes = FALSE, bg = "grey70")

title(xlab = expression(bold("Individual ID")), ylab = expression(bold("Definitive ID")), line = 13, cex.axis = 1.2)

# x axis names
xaxis.names <- sp.idd
axis(1, seq(0,1, length.out = length(xaxis.names))[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], xaxis.names[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], las = 2, font = 3, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[-grep("^[A-Z]", xaxis.names)], xaxis.names[-grep("^[A-Z]", xaxis.names)], las = 2, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[str_count(xaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", xaxis.names[str_count(xaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", xaxis.names[str_count(xaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Individual ID
axis(3, seq(0,1, length.out = length(xaxis.names)) - 0.002, table(IDs.long$ID[IDs.long$MeanDia >= 200])[match(xaxis.names, names(table(IDs.long$ID)))], cex.axis = 0.8, las = 2, tick = FALSE, mgp = c(3, 0.5, 0))

# y axis names
yaxis.names <- sp.idd[rowSums(conf.sp$table) > 0]
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)), labels = FALSE)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 1], yaxis.names[str_count(yaxis.names, " ") == 1], las = 1, font = 3, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[-grep("^[A-Z]", yaxis.names)], yaxis.names[-grep("^[A-Z]", yaxis.names)], las = 1, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", yaxis.names[str_count(yaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", yaxis.names[str_count(yaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Definitive ID
axis(4, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)) - 0.002, ifelse(is.na(tapply(sp.size$SpecNum[sp.size$MeanDia >= 200], completeIDs$DefinitiveID[sp.size$MeanDia >= 200], length)[match(yaxis.names, names(table(IDs.long$DefinitiveID)))]), 0, tapply(sp.size$SpecNum[sp.size$MeanDia >= 200], completeIDs$DefinitiveID[sp.size$MeanDia >= 200], length)[match(yaxis.names, names(table(IDs.long$DefinitiveID)))]), cex.axis = 0.8, las = 1, tick = FALSE, mgp = c(3, 0.5, 0))

# add key
par(fig = c(0.9, 1, 0.35, 1), new = TRUE)
par(mai = c(1, 0.4, 0.8, 0.5))
names.key <- rep("", 101)
names.key[seq(1, 101, by = 20)] <- seq(0.0, 1.0, by = 0.2)
key <- rep(1, 101)
barplot(key, main = "\nFraction of\nSpecimens", horiz = TRUE, space = 0, border = NA, col = c("grey70", matlab.like(100)), fg = "white", cex.main = 1.1, font.main = 1, xaxt = "n", yaxt = "n")
axis(4, at = 1:101, names.key, las = 1, mgp = c(0, 1, 0), xaxt = "n", cex.axis = 1, tick = FALSE)
par(par.def)
dev.off() 

conf.large$overall

# confusion matrix for small specimens 
conf.small <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$MeanDia < 200], levels = sp.idd), factor(IDs.long$ID[IDs.long$MeanDia < 200], levels = sp.idd))
conf.small.frac <- t(conf.small$table / rowSums(conf.small$table))
conf.small.frac[conf.small.frac == "NaN"] <- 0 # otherwise white patches appear across the plots

png("Figures/confusion_small_grey_key.png", 1000, 700)
par(fig = c(0, 0.9, 0, 1))
par(mar = c(15, 15, 2, 2))
# confusion matrix
image(conf.small.frac, col = c("grey70", matlab.like(1000)), ylim = c((28*((1+1/43))/44) - (1+1/43)/44/2, -(1+1/43)/44/2 ), axes = FALSE, bg = "grey70")

title(xlab = expression(bold("Individual ID")), ylab = expression(bold("Definitive ID")), line = 13, cex.axis = 1.2)

# x axis names
xaxis.names <- sp.idd
axis(1, seq(0,1, length.out = length(xaxis.names))[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], xaxis.names[intersect(which(str_count(xaxis.names, " ") != 2),grep("^[A-Z]", xaxis.names))], las = 2, font = 3, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[-grep("^[A-Z]", xaxis.names)], xaxis.names[-grep("^[A-Z]", xaxis.names)], las = 2, cex.axis = 1.05)
axis(1, seq(0,1, length.out = length(xaxis.names))[str_count(xaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", xaxis.names[str_count(xaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", xaxis.names[str_count(xaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Individual ID
axis(3, seq(0,1, length.out = length(xaxis.names)) - 0.002, table(IDs.long$ID[IDs.long$MeanDia < 200])[match(xaxis.names, names(table(IDs.long$ID)))], cex.axis = 0.8, las = 2, tick = FALSE, mgp = c(3, 0.5, 0))

# y axis names
yaxis.names <- sp.idd[rowSums(conf.sp$table) > 0]
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)), labels = FALSE)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 1], yaxis.names[str_count(yaxis.names, " ") == 1], las = 1, font = 3, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[-grep("^[A-Z]", yaxis.names)], yaxis.names[-grep("^[A-Z]", yaxis.names)], las = 1, cex.axis = 1.05)
axis(2, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names))[str_count(yaxis.names, " ") == 2], labels = parse(text = paste(paste("italic('", gsub("^([^ ]* [^ ]*) (.*$)", "\\1", yaxis.names[str_count(yaxis.names, " ") == 2]), "')", sep = ""), gsub("^([^ ]* [^ ]*) (.*$)", "\\2", yaxis.names[str_count(yaxis.names, " ") == 2]), sep = "~")), las = 2, cex.axis = 1.05)

# number of specimens in Definitive ID
axis(4, seq(0,(27*((1+1/43))/44), length.out = length(yaxis.names)) - 0.002, ifelse(is.na(tapply(sp.size$SpecNum[sp.size$MeanDia < 200], completeIDs$DefinitiveID[sp.size$MeanDia < 200], length)[match(yaxis.names, names(table(IDs.long$DefinitiveID)))]), 0, tapply(sp.size$SpecNum[sp.size$MeanDia < 200], completeIDs$DefinitiveID[sp.size$MeanDia < 200], length)[match(yaxis.names, names(table(IDs.long$DefinitiveID)))]), cex.axis = 0.8, las = 1, tick = FALSE, mgp = c(3, 0.5, 0))

# add key
par(fig = c(0.9, 1, 0.35, 1), new = TRUE)
par(mai = c(1, 0.4, 0.8, 0.5))
names.key <- rep("", 101)
names.key[seq(1, 101, by = 20)] <- seq(0.0, 1.0, by = 0.2)
key <- rep(1, 101)
barplot(key, main = "\nFraction of\nSpecimens", horiz = TRUE, space = 0, border = NA, col = c("grey70", matlab.like(100)), fg = "white", cex.main = 1.1, font.main = 1, xaxt = "n", yaxt = "n")
axis(4, at = 1:101, names.key, las = 1, mgp = c(0, 1, 0), xaxt = "n", cex.axis = 1, tick = FALSE)
par(par.def)
dev.off() 


# 10. Producing summaries for people ---------------------------------------
write.csv(completeIDs[, c("SpecNumber", "DefinitiveID", "ChealesID", "ChealesC", "fracCorr")], file = "Outputs/Repeatability_Cheales.csv", row.names = FALSE)

# create personal summaries
for (i in people[3:length(people)]) {
  # extract the data
  tmp <- completeIDs[, c("SpecNumber", "DefinitiveID", grep(i, names(completeIDs), value = TRUE), "fracCorr")]
  tmp$Correct <- ifelse(tmp$DefinitiveID == tmp[, grep(paste(i, "ID", sep = ""), names (tmp))], "Yes", "No")
  tmp$Correct <- ifelse(tmp[, grep(paste(i, "ID", sep = ""), names (tmp))] == "lost", NA, tmp$Correct)
  write.csv(tmp, file = paste("Outputs/Repeatability_", i, ".csv", sep = ""), row.names = FALSE)
}

# 11. Modelling mortality -------------------------------------------------
mm <- glm((ID == "lost") ~ logMeanDia*DefinitiveID, data = IDs.long, family = "binomial")

