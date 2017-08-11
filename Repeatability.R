# Repeatability analysis
# Project: Repeatability
# Author: Isabel Fenton
# Date created: 9/3/2017
# Date last edited: 2/8/2017
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

rm(list = ls())

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

# 6b. adding in the specimen level size data ----------------------------------
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

# 6c. adding in the species level data -------------------------------
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

m1.fix <- glm(Corr ~ logMeanDia*(HowLong + Taught + ED + Conf + SpConf + Experience + Gender + DefinitiveID) + Person + SpecNumber, family = binomial, data = IDs.long.mod) 
par(mfrow = c(2,2))
plot(m1.fix)
par(mfrow = c(1,1))
summary(m1.fix)

stepped <- step(m1.fix)
summary(stepped)

sjp.glm(stepped, type = "pred", vars = c("logMeanDia", "DefinitiveID"), geom.colors = matlab.like(30))

# optimal random structure
glmcon <- glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=3e5))
m1.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|DefinitiveID) + (1|Person) + (1|SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

sjp.glm(m1.me, type = "pred", vars = c("logMeanDia", "DefinitiveID"), geom.colors = matlab.like(30))

m2.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (csLogMeanDia|DefinitiveID) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

m3.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (csLogMeanDia|DefinitiveID) + (1|SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

m4.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (csLogMeanDia|DefinitiveID), family = binomial, data = IDs.long.mod, control = glmcon)

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
summary(m5a.me)
anova(m5.me, m5a.me)
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

# can't simplify further, so move to plotting
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
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "HowLong"))
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Taught"))
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Conf"))
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "SpConf"))
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "Experience"))
png("Figures/GLM_size_sp.png", 1500, 1200, res = 120)
sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "DefinitiveID"), geom.colors = rep("blue4", 26), show.ci = TRUE, cex = 10, cex.axis = 10, cex.labels = 10)
dev.off()

for (i in 1:26){
  png(paste("Figures/GLM_size_", i, ".png", sep = ""), 900, 500)
  sjp.glmer(m5f.me, type = "pred", vars = c("csLogMeanDia", "DefinitiveID"), facet.grid = FALSE, pch = ".", geom.colors = c(rep("grey90",i-1), "black", rep("grey90", 27-i)))
  dev.off()
}

sjp.int(m5f.me, type = "eff")

# 8. Create a confusion matrix --------------------------------------------
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
dev.off()



conf.sp$overall

conf.exp <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Conf == "y"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Conf == "y"], levels = sp.idd))
image(conf.exp$table / rowSums(conf.exp$table), axes = FALSE, col = matlab.like(100))
axis(1, seq(0,1, length.out = length(sp.idd)), sp.idd, las = 2)
axis(2, seq(0,1, length.out = length(sp.idd)), sp.idd, las = 1)
conf.exp$overall


conf.conf <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Conf == "y"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Conf == "y"], levels = sp.idd))
image(conf.conf$table / rowSums(conf.conf$table), axes = FALSE, col = matlab.like(100))
axis(1, seq(0,1, length.out = length(sp.idd)), sp.idd, las = 2)
axis(2, seq(0,1, length.out = length(sp.idd)), sp.idd, las = 1)
conf.conf$overall
 
conf.unconf <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Conf == "n"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Conf == "n"], levels = sp.idd))
image(conf.unconf$table / rowSums(conf.unconf$table), axes = FALSE, col = matlab.like(100))
axis(1, seq(0,1, length.out = length(sp.idd)), sp.idd, las = 2)
axis(2, seq(0,1, length.out = length(sp.idd)), sp.idd, las = 1)
conf.unconf$overall


# 9. Producing summaries for people ---------------------------------------
write.csv(completeIDs[, c("SpecNumber", "DefinitiveID", "ChealesID", "ChealesC", "fracCorr")], file = "Outputs/Repeatability_Cheales.csv", row.names = FALSE)


