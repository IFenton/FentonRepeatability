# Repeatability analysis
# Project: Repeatability
# Author: Isabel Fenton
# Date created: 8/3/2018
# Date last edited: 23/4/2018
# 
# Code for analysis of the repeatability paper
# Fenton et al on foram repeatability
# 
# Previous file: Repeatability.R

rm(list = ls())

# Inputs ------------------------------------------------------------------
# assumed to be in a folder /Data
# Fenton_CompleteIDs.csv - the raw data
# Fenton_SpeciesData.csv - 
# Fenton_PersonData.csv - information about the participants
# Fenton_MeasurementsFull.csv - the size measurements
# Fenton_Specieslist.csv - list of species

# Outputs -----------------------------------------------------------------
# F_Figures are produced in a folder /F_F_Figures

# Source files / libraries ------------------------------------------------
library(cluster) # for the dendrogram
library(caret) # for the confusion matrix
library(colorRamps) # colours
library(ade4) # for distinctiveness
library(sjPlot) # plotting glms
library(lme4) # for modelling
library(stringr) # for confusion matrix axes
library(car) # p-values in the anova
library(MuMIn) # r2 values for the model
source("Code/Confusion_matrix.R")

# 1. Load in the data -----------------------------------------------------
SpeciesData <- read.csv("Data/Fenton_SpeciesData.csv")
str(SpeciesData)

completeIDs <- read.csv("Data/Fenton_CompleteIDs.csv")
str(completeIDs)

personData <- read.csv("Data/Fenton_PersonData.csv", check.names = FALSE) # to keep the spaces etc. in the column names
str(personData)

sp.size <- read.csv("Data/Fenton_MeasurementsFull.csv")
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
png("F_Figures/FD_dend.png", 800, 800)
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

rm(sp.fd, sp.dend, sp.dist)

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
sum(completeIDs$Exp3ID != "lost") # exclude lost specimens
sum(completeIDs$Exp3ID != "lost"& completeIDs$Exp3ID %in% SpeciesData$Species) # excluded lost and those not at species level
people.df$numLost <- apply(completeIDs[, people.col], 2, function(x) sum(x == "lost"))
people.df$TotalIndSp <- apply(completeIDs[, people.col], 2, function(x) sum(x != "lost" & x %in% SpeciesData$Species))
people.df$TotalIndGe <- apply(completeIDs[, people.col], 2, function(x) sum(x != "lost"))

# Fraction of species identified
people.df$fracSpecIDd <- people.df$NumSpID / people.df$NumSpID[people.df$Name == "Definitive"]

# Fraction of genera identified
people.df$fracGenIDd <- people.df$NumGenID / people.df$NumGenID[people.df$Name == "Definitive"]


## Percentage accuracy per person
sum(completeIDs[,grep("Exp1", names(completeIDs))[1]] == completeIDs$DefinitiveID) / sum(completeIDs[,grep("Exp1", names(completeIDs))[1]] != "lost") # fraction that match the definitive ID and weren't lost

# species
people.df$fracCorrSp <- apply(completeIDs[, people.col], 2, function(x) sum(x == completeIDs$DefinitiveID) / sum(x != "lost"))
mean(people.df$fracCorrSp) # 0.63
median(people.df$fracCorrSp[!is.na(people.df$Experienced)]) # 0.59
median(people.df$fracCorrSp[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.79
median(people.df$fracCorrSp[people.df$Experienced == "Student"], na.rm = TRUE) # 0.57

# genus
gsub(" .*", "", completeIDs$DefinitiveID) # list the genus

sum(gsub(" .*", "", completeIDs$Exp1ID) == gsub(" .*", "", completeIDs$DefinitiveID)) / sum(completeIDs[,grep("Exp1", names(completeIDs))[1]] != "lost") # fraction of specimens that match at genus level and weren't lost

people.df$fracCorrGen <- apply(completeIDs[, people.col], 2, function(x) sum(gsub(" .*", "", x) == gsub(" .*", "", completeIDs$DefinitiveID)) / sum(x != "lost"))
mean(people.df$fracCorrGen) # 0.76

# number of specimens unidentified
people.df$numUnID <- NA
people.df$numUnID <- apply(completeIDs[, people.col], 2, function(x) sum(x == "unIDd", na.rm = TRUE))

# box plot of percentage accuracy split by students and Experienced
png("F_Figures/Percent accuracy.png")
with(people.df, boxplot(fracCorrSp ~ Experienced, ylim = c(0, 1), xlim = c(-0.5, 2.5), xaxt = "n", boxwex = 0.4, col = c("blue3", "purple2"), cex.axis = 1.5, las = 1))
with(people.df[!is.na(people.df$Experienced), ], boxplot(fracCorrSp, add = TRUE, at = 0, boxwex = 0.8, col = "red2", xaxt = "n", yaxt = "n"))
axis(1, c(0,1,2), c("Full", levels(people.df$Experienced)), cex.axis = 1.5)
text(0:2, 0.95, paste("(", c(sum(!is.na(people.df$Experienced)), sum(people.df$Experienced == "Experienced", na.rm = TRUE), sum(people.df$Experienced == "Student", na.rm = TRUE)), ")", sep = ""))
dev.off()

## Split by confidence 
# number of identifications by confidence (working at specimen level)
sum(completeIDs$Exp1C == "y", na.rm = TRUE)

people.df$numCy <- NA
people.df$numCy[3:nrow(people.df)] <- apply(completeIDs[, p.conf.col], 2, function(x) sum(x == "y", na.rm = TRUE))

people.df$numCm <- NA
people.df$numCm[3:nrow(people.df)] <- apply(completeIDs[, p.conf.col], 2, function(x) sum(x == "m", na.rm = TRUE))

people.df$numCn <- NA
people.df$numCn[3:nrow(people.df)] <- apply(completeIDs[, p.conf.col], 2, function(x) sum(x == "n", na.rm = TRUE))

head(people.df)

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
rm(i)
# where there are no specimens classified as 'n', this is currently giving NaN
people.df[is.na(people.df)] <- NA

hist(people.df$fracCorrY)
hist(people.df$fracCorrM)
hist(people.df$fracCorrN)
hist(people.df$fracCorrNwoUID)

# percentage accuracy by confidence
png("F_Figures/Percent accuracy conf.png")
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
median(people.df$fracCorrY[!is.na(people.df$Experienced)]) # 0.77
median(people.df$fracCorrY[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.93
median(people.df$fracCorrY[people.df$Experienced == "Student"], na.rm = TRUE) # 0.75

# percent correct if maybe
median(people.df$fracCorrM[!is.na(people.df$Experienced)]) # 0.44
median(people.df$fracCorrM[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.67
median(people.df$fracCorrM[people.df$Experienced == "Student"], na.rm = TRUE) # 0.41

# percent correct if no
median(people.df$fracCorrNwoUID[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.25
median(people.df$fracCorrNwoUID[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.25
median(people.df$fracCorrNwoUID[people.df$Experienced == "Student"], na.rm = TRUE) # 0.27


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
rm(i)

png("F_Figures/Percent accuracy species level.png")
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
median(people.df$fracCorrSY[!is.na(people.df$Experienced)]) # 0.77
median(people.df$fracCorrSY[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.85
median(people.df$fracCorrSY[people.df$Experienced == "Student"], na.rm = TRUE) # 0.75

# percent correct if maybe
median(people.df$fracCorrSM[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.78
median(people.df$fracCorrSM[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.69
median(people.df$fracCorrSM[people.df$Experienced == "Student"], na.rm = TRUE) # 0.78

# percent correct if no
median(people.df$fracCorrSN[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.32
median(people.df$fracCorrSN[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.33
median(people.df$fracCorrSN[people.df$Experienced == "Student"], na.rm = TRUE) # 0.31


## accuracy based on species and specimen level confidence
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
rm(i)

png("F_Figures/Percent accuracy species level.png")
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
median(people.df$fracCorrYy[!is.na(people.df$Experienced)]) # 0.84
median(people.df$fracCorrYy[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.97
median(people.df$fracCorrYy[people.df$Experienced == "Student"], na.rm = TRUE) # 0.84
# percent correct if Yes maybe
median(people.df$fracCorrYm[!is.na(people.df$Experienced)]) # 0.6
median(people.df$fracCorrYm[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.68
median(people.df$fracCorrYm[people.df$Experienced == "Student"], na.rm = TRUE) # 0.6
# percent correct if Yes no
median(people.df$fracCorrYn[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.05
median(people.df$fracCorrYn[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0
median(people.df$fracCorrYn[people.df$Experienced == "Student"], na.rm = TRUE) # 0.15

# percent correct if Maybe yes
median(people.df$fracCorrMy[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.47
median(people.df$fracCorrMy[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.47
median(people.df$fracCorrMy[people.df$Experienced == "Student"], na.rm = TRUE) # 0.44
# percent correct if Maybe maybe
median(people.df$fracCorrMm[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.86
median(people.df$fracCorrMm[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.80
median(people.df$fracCorrMm[people.df$Experienced == "Student"], na.rm = TRUE) # 0.86
# percent correct if Maybe no
median(people.df$fracCorrMn[!is.na(people.df$Experienced)], na.rm = TRUE) # 1
median(people.df$fracCorrMn[people.df$Experienced == "Experienced"], na.rm = TRUE) # NA
median(people.df$fracCorrMn[people.df$Experienced == "Student"], na.rm = TRUE) # 1

# percent correct if No yes
median(people.df$fracCorrNy[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.5
median(people.df$fracCorrNy[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.94
median(people.df$fracCorrNy[people.df$Experienced == "Student"], na.rm = TRUE) # 0.48
# percent correct if No maybe
median(people.df$fracCorrNm[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.33
median(people.df$fracCorrNm[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.5
median(people.df$fracCorrNm[people.df$Experienced == "Student"], na.rm = TRUE) # 0.32
# percent correct if No no
median(people.df$fracCorrNn[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.12
median(people.df$fracCorrNn[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.25
median(people.df$fracCorrNn[people.df$Experienced == "Student"], na.rm = TRUE) # 0.11

# how many people does this apply to
sum(!is.na(people.df$fracCorrYy)) # 23
sum(!is.na(people.df$fracCorrYm)) # 23
sum(!is.na(people.df$fracCorrYn)) # 18
sum(!is.na(people.df$fracCorrMy)) # 4
sum(!is.na(people.df$fracCorrMm)) # 5
sum(!is.na(people.df$fracCorrMn)) # 1
sum(!is.na(people.df$fracCorrNy)) # 21
sum(!is.na(people.df$fracCorrNm)) # 21
sum(!is.na(people.df$fracCorrNn)) # 20
sum(!is.na(people.df$fracCorrYy[people.df$Experienced == "Experienced"])) # 4
sum(!is.na(people.df$fracCorrYm[people.df$Experienced == "Experienced"])) # 4
sum(!is.na(people.df$fracCorrYn[people.df$Experienced == "Experienced"])) # 4
sum(!is.na(people.df$fracCorrMy[people.df$Experienced == "Experienced"])) # 2
sum(!is.na(people.df$fracCorrMm[people.df$Experienced == "Experienced"])) # 2
sum(!is.na(people.df$fracCorrMn[people.df$Experienced == "Experienced"])) # 0
sum(!is.na(people.df$fracCorrNy[people.df$Experienced == "Experienced"])) # 2
sum(!is.na(people.df$fracCorrNm[people.df$Experienced == "Experienced"])) # 3
sum(!is.na(people.df$fracCorrNn[people.df$Experienced == "Experienced"])) # 2
sum(!is.na(people.df$fracCorrYy[people.df$Experienced == "Student"])) # 19
sum(!is.na(people.df$fracCorrYm[people.df$Experienced == "Student"])) # 19
sum(!is.na(people.df$fracCorrYn[people.df$Experienced == "Student"])) # 14
sum(!is.na(people.df$fracCorrMy[people.df$Experienced == "Student"])) # 2
sum(!is.na(people.df$fracCorrMm[people.df$Experienced == "Student"])) # 3
sum(!is.na(people.df$fracCorrMn[people.df$Experienced == "Student"])) # 1
sum(!is.na(people.df$fracCorrNy[people.df$Experienced == "Student"])) # 19
sum(!is.na(people.df$fracCorrNm[people.df$Experienced == "Student"])) # 18
sum(!is.na(people.df$fracCorrNn[people.df$Experienced == "Student"])) # 18

# what is the median of the number of specimens
# number of Yes yes
median(people.df$numcYy[!is.na(people.df$Experienced) & people.df$numcYy > 0]) # 38
median(people.df$numcYy[people.df$Experienced == "Experienced" & people.df$numcYy > 0], na.rm = TRUE) # 46.5
median(people.df$numcYy[people.df$Experienced == "Student" & people.df$numcYy > 0], na.rm = TRUE) # 33
# number of Yes maybe
median(people.df$numcYm[!is.na(people.df$Experienced) & people.df$numcYm > 0]) # 13
median(people.df$numcYm[people.df$Experienced == "Experienced" & people.df$numcYm > 0], na.rm = TRUE) # 16.5
median(people.df$numcYm[people.df$Experienced == "Student" & people.df$numcYm > 0], na.rm = TRUE) # 11
# number of Yes no
median(people.df$numcYn[!is.na(people.df$Experienced) & people.df$numcYn > 0], na.rm = TRUE) # 3
median(people.df$numcYn[people.df$Experienced == "Experienced" & people.df$numcYn > 0], na.rm = TRUE) # 1.5
median(people.df$numcYn[people.df$Experienced == "Student" & people.df$numcYn > 0], na.rm = TRUE) # 4.5

# number of Maybe yes
median(people.df$numcMy[!is.na(people.df$Experienced) & people.df$numcMy > 0], na.rm = TRUE) # 4
median(people.df$numcMy[people.df$Experienced == "Experienced" & people.df$numcMy > 0], na.rm = TRUE) # 4
median(people.df$numcMy[people.df$Experienced == "Student" & people.df$numcMy > 0], na.rm = TRUE) # 5
# number of Maybe maybe
median(people.df$numcMm[!is.na(people.df$Experienced) & people.df$numcMm > 0], na.rm = TRUE) # 7
median(people.df$numcMm[people.df$Experienced == "Experienced" & people.df$numcMm > 0], na.rm = TRUE) # 8
median(people.df$numcMm[people.df$Experienced == "Student" & people.df$numcMm > 0], na.rm = TRUE) # 7
# number of Maybe no
median(people.df$numcMn[!is.na(people.df$Experienced) & people.df$numcMn > 0], na.rm = TRUE) # 1
median(people.df$numcMn[people.df$Experienced == "Experienced" & people.df$numcMn > 0], na.rm = TRUE) # NA
median(people.df$numcMn[people.df$Experienced == "Student" & people.df$numcMn > 0], na.rm = TRUE) # 1

# number of No yes
median(people.df$numcNy[!is.na(people.df$Experienced) & people.df$numcNy > 0], na.rm = TRUE) # 8
median(people.df$numcNy[people.df$Experienced == "Experienced" & people.df$numcNy > 0], na.rm = TRUE) # 4.5
median(people.df$numcNy[people.df$Experienced == "Student" & people.df$numcNy > 0], na.rm = TRUE) # 8
# number of No maybe
median(people.df$numcNm[!is.na(people.df$Experienced) & people.df$numcNm > 0], na.rm = TRUE) # 9
median(people.df$numcNm[people.df$Experienced == "Experienced" & people.df$numcNm > 0], na.rm = TRUE) # 2
median(people.df$numcNm[people.df$Experienced == "Student" & people.df$numcNm > 0], na.rm = TRUE) # 9.5
# number of No no
median(people.df$numcNn[!is.na(people.df$Experienced) & people.df$numcNn > 0], na.rm = TRUE) # 8
median(people.df$numcNn[people.df$Experienced == "Experienced" & people.df$numcNn > 0], na.rm = TRUE) # 4
median(people.df$numcNn[people.df$Experienced == "Student" & people.df$numcNn > 0], na.rm = TRUE) # 8.5


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
rm(i)

png("F_Figures/Percent accuracy size.png")
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
median(people.df$fracCorr125[!is.na(people.df$Experienced)]) # 0.43
median(people.df$fracCorr125[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.63
median(people.df$fracCorr125[people.df$Experienced == "Student"], na.rm = TRUE) # 0.37

# percent correct if 200-400
median(people.df$fracCorr200[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.55
median(people.df$fracCorr200[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.75
median(people.df$fracCorr200[people.df$Experienced == "Student"], na.rm = TRUE) # 0.53

# percent correct if >400
median(people.df$fracCorr400[!is.na(people.df$Experienced)], na.rm = TRUE) # 0.76
median(people.df$fracCorr400[people.df$Experienced == "Experienced"], na.rm = TRUE) # 0.96
median(people.df$fracCorr400[people.df$Experienced == "Student"], na.rm = TRUE) # 0.76

# 3c. Plotting confidence -------------------------------------------------
plot(c(1,3), c(0, 1), type = "n", xaxt = "n", xlab = "", ylab = "", las = 2, xlim = c(0.5, 3.5), main = "Specimen confidence")
axis(1, 1:3, labels = c("no", "maybe", "yes"))
for (i in 3:nrow(people.df)) {
  tmp.y <- people.df[i, c("fracCorrNwoUID", "fracCorrM", "fracCorrY")]
  points((1:3)[!is.na(tmp.y)], tmp.y[!is.na(tmp.y)], type = "b", col = rep(1:5, 10)[i])
}
rm(i, tmp.y)

plot(c(1,3), c(0, 1), type = "n", xaxt = "n", xlab = "", ylab = "", las = 2, xlim = c(0.5, 3.5), main = "Species confidence")
axis(1, 1:3, labels = c("no", "maybe", "yes"))
for (i in 1:nrow(people.df)) {
  tmp.y <- people.df[i, c("fracCorrSN", "fracCorrSM", "fracCorrSY")]
  points((1:3)[!is.na(tmp.y)], tmp.y[!is.na(tmp.y)], type = "b", col = rep(1:5, 10)[i])
}
rm(i, tmp.y)


# 4. Percentage accuracy per specimen ------------------------------------
sum(completeIDs[2, people.col[2:length(people.col)]] == as.character(completeIDs$DefinitiveID[2]) & completeIDs[2, people.col[2:length(people.col)]] != "lost") / sum(completeIDs[2, people.col[2:length(people.col)]] != "lost")

completeIDs$fracCorr <- apply(completeIDs[, people.col[2:length(people.col)]], 1, function(x) sum(x == x[1] & x != "lost") / sum(x != "lost"))

hist(completeIDs$fracCorr)

# how does accuracy depend on the species?
par(mfrow = c(2,1))
plot(factor(completeIDs$DefinitiveID), completeIDs$fracCorr, las = 2)
par(mfrow = c(1,1))

# 5. Reshaping the dataset for modelling -----------------------------------

# 5a. Creating a long form dataframe --------------------------------------
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

# 5b. Adding in the specimen level size data ----------------------------------
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

# 5c. Adding in the species level data -------------------------------
head(SpeciesData)
IDs.long <- merge(IDs.long, SpeciesData[, 2:ncol(SpeciesData)], by.x = "DefinitiveID", by.y = "Species")
head(IDs.long)

# 5d. Adding in the person level data -------------------------------------
# convert how long into an ordered factor for modelling based on quantiles
summary(personData$"How long have you worked with forams (years)")
quantile(personData$"How long have you worked with forams (years)", c(0.25, 0.5, 0.75))

personData$HowLong <- 0
personData$HowLong[personData$"How long have you worked with forams (years)" >= quantile(personData$"How long have you worked with forams (years)", c(0.25)) & personData$"How long have you worked with forams (years)" < quantile(personData$"How long have you worked with forams (years)", c(0.5))] <- 1
personData$HowLong[personData$"How long have you worked with forams (years)" >= quantile(personData$"How long have you worked with forams (years)", c(0.5)) & personData$"How long have you worked with forams (years)" < quantile(personData$"How long have you worked with forams (years)", c(0.75))] <- 2
personData$HowLong[personData$"How long have you worked with forams (years)" >= quantile(personData$"How long have you worked with forams (years)", c(0.75))] <- 3
# check this worked
cbind(personData$"How long have you worked with forams (years)", personData$HowLong)

# adding in person level data
head(personData)
IDs.long <- merge(IDs.long, personData)
head(IDs.long)

# adding in the species level confidence by person
IDs.long$SpConf <- apply(IDs.long, 1, function (x) as.character(x[x["DefinitiveID"]]))
IDs.long$SpConf <- factor(IDs.long$SpConf)


# 6. Modelling ------------------------------------------------------------

# 6a. Create a dataframe for modelling ------------------------------------
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
pairs(IDs.long.mod[, c("HowLong", "logMeanDia", "Taught", "ED", "Conf", "SpConf", "Experience", "Gender")]) # nothing is obviously correlated


# 6b. Fixed model ---------------------------------------------------------
# run a fixed model first with everything, to get an approximate idea
m1.fix <- glm(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender + DefinitiveID) + Person + SpecNumber, family = binomial, data = IDs.long.mod) # you get convergence warnings on this (I think because they aren't fitted as randoms.)
par(mfrow = c(2,2))
plot(m1.fix)
par(mfrow = c(1,1))
summary(m1.fix)
rm(m1.fix)

# 6c. Optimal model -------------------------------------------------------
# optimal random structure
opt.RS <- list()
glmcon <- glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=3e5)) # stop convergence issues

# test the different combination of random effects
opt.RS$m1.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|DefinitiveID/SpecNumber) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

sjp.glm(opt.RS$m1.me, type = "pred", vars = c("csLogMeanDia", "DefinitiveID"), geom.colors = matlab.like(30))

opt.RS$m2.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|DefinitiveID/SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m3.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m4.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|DefinitiveID/SpecNumber) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m5.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|DefinitiveID/SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m6.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|DefinitiveID) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m7.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|DefinitiveID), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m8.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|DefinitiveID) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m9.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|DefinitiveID), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m10.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|SpecNumber) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m11.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1 + csLogMeanDia|SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m12.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|SpecNumber) + (1|Person), family = binomial, data = IDs.long.mod, control = glmcon)

opt.RS$m13.me <- glmer(Corr ~ csLogMeanDia*(HowLong + Taught + csED + Conf + SpConf + Experience + Gender) + (1|SpecNumber), family = binomial, data = IDs.long.mod, control = glmcon)


anova(opt.RS$m1.me, opt.RS$m2.me, opt.RS$m3.me, opt.RS$m4.me, opt.RS$m5.me, opt.RS$m6.me, opt.RS$m7.me, opt.RS$m8.me, opt.RS$m9.me, opt.RS$m10.me, opt.RS$m11.me, opt.RS$m12.me, opt.RS$m13.me)

# best model is opt.RS$m1.me
summary(opt.RS$m1.me)

# output this result
write.csv(anova(opt.RS$m1.me, opt.RS$m2.me, opt.RS$m3.me, opt.RS$m4.me, opt.RS$m5.me, opt.RS$m6.me, opt.RS$m7.me, opt.RS$m8.me, opt.RS$m9.me, opt.RS$m10.me, opt.RS$m11.me, opt.RS$m12.me, opt.RS$m13.me)
, "F_Outputs/randomErrorComp.csv")

opt.RS$re.m1.me <- residuals(opt.RS$m1.me, type = "response")
opt.RS$fi.m1.me <- fitted(opt.RS$m1.me)
plot(opt.RS$fi.m1.me, opt.RS$re.m1.me)

# run model simplification
opt.MS <- list()
opt.MS$m1a.me <- update(opt.RS$m1.me, .~. -csLogMeanDia:Gender)
anova(opt.RS$m1.me, opt.MS$m1a.me)
summary(opt.MS$m1a.me)
opt.MS$m1b.me <- update(opt.MS$m1a.me, .~. -csLogMeanDia:SpConf)
anova(opt.MS$m1b.me, opt.MS$m1a.me)
summary(opt.MS$m1b.me)
opt.MS$m1c.me <- update(opt.MS$m1b.me, .~. -csLogMeanDia:csED)
anova(opt.MS$m1b.me, opt.MS$m1c.me)
summary(opt.MS$m1c.me)
opt.MS$m1d.me <- update(opt.MS$m1c.me, .~. -csED)
anova(opt.MS$m1d.me, opt.MS$m1c.me)
summary(opt.MS$m1d.me)
opt.MS$m1e.me <- update(opt.MS$m1d.me, .~. -csLogMeanDia:Conf)
anova(opt.MS$m1d.me, opt.MS$m1e.me)
summary(opt.MS$m1e.me)
opt.MS$m1f.me <- update(opt.MS$m1e.me, .~. -Gender)
anova(opt.MS$m1f.me, opt.MS$m1e.me)
summary(opt.MS$m1f.me)

anova(opt.MS$m1f.me)
Anova(opt.MS$m1f.me)
r.squaredGLMM(opt.MS$m1f.me) # R2m: 0.426, R2c: 0.571

write.csv(Anova(opt.MS$m1f.me), file = "F_Outputs/anova_m1f_me.csv")


# 6d. Marginal effects of the different EVs -------------------------------
summary(opt.MS$m1f.me)
r.squaredGLMM(opt.MS$m1f.me) # r2m = 0.4263517

# drop each EV in turn to get the marginal effect
mar.eff <- list()

# create a table of this data
mar.eff.dat <- data.frame(EV = c("FullModel", "csLogMeanDia", "HowLong", "Taught", "Experience", "Conf", "SpConf"))
mar.eff.dat$r2m <- NA
mar.eff.dat$r2m[mar.eff.dat$EV == "FullModel"] <- r.squaredGLMM(opt.MS$m1f.me)["R2m"]

# log size
mar.eff$m1f.me.cLMD <- update(opt.MS$m1f.me, ~. - csLogMeanDia - csLogMeanDia:HowLong - csLogMeanDia:Taught - csLogMeanDia:Experience)
summary(mar.eff$m1f.me.cLMD)
(mar.eff.dat$r2m[mar.eff.dat$EV == "csLogMeanDia"] <- r.squaredGLMM(mar.eff$m1f.me.cLMD)["R2m"]) # r2m = 0.3840455
0.4263517 - 0.3840455 # 0.0423062
# how long
mar.eff$m1f.me.HL <- update(opt.MS$m1f.me, ~. - HowLong - csLogMeanDia:HowLong)
summary(mar.eff$m1f.me.HL)
(mar.eff.dat$r2m[mar.eff.dat$EV == "HowLong"] <- r.squaredGLMM(mar.eff$m1f.me.HL)["R2m"]) # r2m = 0.3807146
0.4263517 - 0.3807146 # 0.0456371
# taught
mar.eff$m1f.me.T <- update(opt.MS$m1f.me, ~. - Taught - csLogMeanDia:Taught)
summary(mar.eff$m1f.me.T)
(mar.eff.dat$r2m[mar.eff.dat$EV == "Taught"] <- r.squaredGLMM(mar.eff$m1f.me.T)["R2m"]) # r2m = 0.2512889
0.4263517 - 0.2512889 # 0.1750628
# Experience
mar.eff$m1f.me.E <- update(opt.MS$m1f.me, ~. - Experience - csLogMeanDia:Experience)
summary(mar.eff$m1f.me.E)
(mar.eff.dat$r2m[mar.eff.dat$EV == "Experience"] <- r.squaredGLMM(mar.eff$m1f.me.E)["R2m"]) # r2m = 0.4083618
0.4263517 - 0.4083618 # 0.0179899
# Specimen confidence
mar.eff$m1f.me.C <- update(opt.MS$m1f.me, ~. - Conf)
summary(mar.eff$m1f.me.C)
(mar.eff.dat$r2m[mar.eff.dat$EV == "Conf"] <- r.squaredGLMM(mar.eff$m1f.me.C)["R2m"]) # r2m = 0.3733382
0.4263517 - 0.3733382 # 0.0544
# Species confidence
mar.eff$m1f.me.SC <- update(opt.MS$m1f.me, ~. - SpConf)
summary(mar.eff$m1f.me.SC)
(mar.eff.dat$r2m[mar.eff.dat$EV == "SpConf"] <- r.squaredGLMM(mar.eff$m1f.me.SC)["R2m"]) # r2m = 0.4198742
0.4263517 - 0.4198742 # 0.0064775

mar.eff.dat$d.r2m <- mar.eff.dat$r2m[mar.eff.dat$EV == "FullModel"] - mar.eff.dat$r2m
write.csv(mar.eff.dat, "F_Outputs/marginalEffects.csv", row.names = FALSE)

# 6e. Model plotting ------------------------------------------------------
# can't simplify further, so move to plotting
pdf("F_Figures/GLM_effects.pdf")
sjp.glmer(opt.MS$m1f.me)
sjp.glmer(opt.MS$m1f.me, type = "fe")
sjp.glmer(opt.MS$m1f.me, type = "eff")
sjp.int(opt.MS$m1f.me, type = "eff")

p <- sjp.glmer(opt.MS$m1f.me, type = "eff", facet.grid = FALSE, 
             show.ci = TRUE, prnt.plot = FALSE)$plot.list
# plot all marginal effects, as grid, proper x-axes
# also, set margins for this example
plot_grid(p, margin = c(0.3, 0.3, 0.3, 0.3))
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = "csLogMeanDia")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = "HowLong")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = "Taught")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = "Conf")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = "SpConf")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = "Experience")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "HowLong"), title = "How Long")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "Taught"), title = "Taught")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "Conf"), title = "Specimen Confidence")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "SpConf"), title = "Species confidence")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "Experience"), title = "Experience")
dev.off()
rm(p)

# influence of size by species 
png("F_Figures/GLM_size_sp.png", 1500, 1200, res = 120)
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "DefinitiveID"), geom.colors = rep("blue4", 26), show.ci = TRUE, cex = 10, cex.axis = 10, cex.labels = 10)
dev.off()

pdf("F_Figures/GLM_size.pdf", 10,7)
for (i in 1:26){
  sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "DefinitiveID"), facet.grid = FALSE, pch = ".", geom.colors = c(rep("grey90",i-1), "black", rep("grey90", 27-i)))
}
dev.off()
rm(i)

# size by confidence
png("F_Figures/Size_Conf_sjp_glmer.png")
set_theme(theme_classic(), axis.title.size = 1.4, axis.textsize = 1.2, panel.gridcol = "white")
sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "Conf"), title = "Conf", show.ci = TRUE, facet.grid = FALSE, show.legend = FALSE)
dev.off()

# paper plots
# the scale for the x-axis is currently not very informative (centred / scaled log mean diameter)
# so if I want them at 200, 400, 600, 800
(log(200) - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia) # -0.9975478
(log(400) - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia) # 0.3995266
(log(600) - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia) # 1.216763
(log(800) - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia) # 1.796601

png("F_Figures/Size_Taught_sjp_glmer.png")
set_theme(theme_classic(), axis.title.size = 1.4, axis.textsize = 1.2, panel.gridcol = "white")
p <- sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "Taught"), title = "", show.ci = TRUE, facet.grid = FALSE, axis.title = c(expression(paste("Mean Diameter / ", mu, "m")), "Percentage Correct"), geom.size = 1.2, point.alpha = 0.5, prnt.plot = FALSE)
p$plot + scale_x_continuous(breaks = c(-0.9975478, 0.3995266, 1.216763, 1.796601), labels = c(200, 400, 600, 800))
dev.off()
rm(p)

png("F_Figures/Size_Exp_sjp_glmer.png")
set_theme(theme_classic(), axis.title.size = 1.4, axis.textsize = 1.2, panel.gridcol = "white")
p1 <- sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "Experience"), title = "", show.ci = TRUE, facet.grid = FALSE, axis.title = c(expression(paste("Mean Diameter / ", mu, "m")), "Percentage Correct"), geom.size = 1.2, point.alpha = 0.5, prnt.plot = FALSE)
p1$plot + scale_x_continuous(breaks = c(-0.9975478, 0.3995266, 1.216763, 1.796601), labels = c(200, 400, 600, 800))
dev.off()
rm(p1)

png("F_Figures/Size_HowLong_sjp_glmer.png")
set_theme(theme_classic(), axis.title.size = 1.4, axis.textsize = 1.2, panel.gridcol = "white")
p2 <- sjp.glmer(opt.MS$m1f.me, type = "pred", vars = c("csLogMeanDia", "HowLong"), title = "", show.ci = TRUE, facet.grid = FALSE, axis.title = c(expression(paste("Mean Diameter / ", mu, "m")), "Percentage Correct"), geom.size = 1.2, point.alpha = 0.5, prnt.plot = FALSE)
p2$plot + scale_x_continuous(breaks = c(-0.9975478, 0.3995266, 1.216763, 1.796601), labels = c(200, 400, 600, 800))
dev.off()
rm(p2)


# 6f. Accuracy under a higher powered microscope ---------------------------
new.IDs.long.mod <- IDs.long.mod

# calculate sizes as if we were working on 50x rather than 35x magnification
tmp <- new.IDs.long.mod$logMeanDia/35*50
# rescale using the old mean / SD to make sure the data is on the right scale

new.IDs.long.mod$csLogMeanDia <- (tmp - mean(IDs.long.mod$logMeanDia)) / sd(IDs.long.mod$logMeanDia)
summary(IDs.long.mod$csLogMeanDia)
rm(tmp)

# predict the correctness using this new data
pred.val.rs <- predict(opt.MS$m1f.me, newdata = new.IDs.long.mod, type = "response")
summary(pred.val.rs)

# compare this new data to the old data
pred.val.orig <- predict(opt.MS$m1f.me, type = "response")
summary(pred.val.orig)

## to understand this data, compare the percentage correct from the smaller dataset used in the modelling to the full dataset
# full data
tmp <- people.df$fracCorrSp[3:nrow(people.df)]
names(tmp) <- people.df$Name[3:nrow(people.df)]
tmp <- tmp[order(names(tmp))]
tmp
# modelled dataset
tapply(IDs.long.mod$Corr, IDs.long.mod$Person, function(x) sum(x) / length(x))
plot(tmp, tapply(IDs.long.mod$Corr, IDs.long.mod$Person, function(x) sum(x) / length(x)), xlab = "Full dataset", ylab = "Modelling dataset", pch = 16, col = grepl("Exp", levels(IDs.long.mod$Person)) + 1)
abline(0, 1)
# some differences, which is to be expected as they are slightly different datasets, but generally similar. 

# how do these compare with the model predictions
tapply(pred.val.orig, IDs.long.mod$Person, mean)
plot(tapply(IDs.long.mod$Corr, IDs.long.mod$Person, function(x) sum(x) / length(x)), tapply(pred.val.orig, IDs.long.mod$Person, mean), xlab = "Modelling dataset", ylab = "Predicted accuracy", pch = 16, col = grepl("Exp", levels(IDs.long.mod$Person)) + 1)
abline(0, 1)
# again some differences, but generally similar (again to be expected, as these are modelled results)

# so how much better would the results be if a larger magnification was used. 
mean(pred.val.rs)
mean(pred.val.orig)

# or split by person
tapply(pred.val.rs, IDs.long.mod$Person, mean)
tapply(pred.val.orig, IDs.long.mod$Person, mean)
plot(tapply(pred.val.orig, IDs.long.mod$Person, mean), tapply(pred.val.rs, IDs.long.mod$Person, mean), xlab = "Original predicted accuracy", ylab = "Rescaled predicted accuracy", pch = 16, col = grepl("Exp", levels(IDs.long.mod$Person)) + 1)
abline(0, 1)

# or split by groups
tapply(pred.val.rs, grepl("Exp", IDs.long.mod$Person), mean)
tapply(pred.val.orig, grepl("Exp", IDs.long.mod$Person), mean)
rm(tmp)

# 7. Create confusion matrices --------------------------------------------
# full confusion matrix
sp.idd <- read.csv("Data/Fenton_Specieslist.csv")$DefinitiveID
conf.mat <- list()

png("F_Figures/confusion_full_key.png", 1000, 700)
conf_mat(IDs.long, "ID", "DefinitiveID", axes.same = FALSE, abb.end = c("juvenile", "nonmacro", "unIDd"), sp.exc = "lost", xlab = "Individual ID")
dev.off()
conf.mat$sp <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$ID != "lost"], levels = sp.idd), factor(IDs.long$ID[IDs.long$ID != "lost"], levels = sp.idd))

# confusion matrix for experienced workers
png("F_Figures/confusion_exp_key.png", 900, 700)
conf_mat(IDs.long, "ID", "DefinitiveID", axes.same = FALSE, abb.end = c("juvenile", "nonmacro", "unIDd"), sp.exc = "lost", subset.col = "Experienced", subset.lev = "Experienced", xlab = "Individual ID")
dev.off()
conf.mat$exp <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Experienced == "Experienced"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Experienced == "Experienced"], levels = sp.idd))

# confusion matrix for students
png("F_Figures/confusion_stu_key.png", 1000, 700)
conf_mat(IDs.long, "ID", "DefinitiveID", axes.same = FALSE, abb.end = c("juvenile", "nonmacro", "unIDd"), sp.exc = "lost", xlab = "Individual ID")
dev.off()
conf.mat$stu <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Experienced == "Student"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Experienced == "Student"], levels = sp.idd))

# confusion matrix for confident IDs
png("F_Figures/confusion_conf_key.png", 900, 720)
conf_mat(IDs.long, "ID", "DefinitiveID", axes.same = FALSE, abb.end = c("juvenile", "nonmacro", "unIDd"), sp.exc = "lost", subset.col = "Conf", subset.lev = "y", xlab = "Individual ID")
dev.off()
conf.mat$conf <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Conf == "y"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Conf == "y"], levels = sp.idd))

# confusion matrix for unconfident IDs
png("F_Figures/confusion_unconf_key.png", 950, 700)
conf_mat(IDs.long, "ID", "DefinitiveID", axes.same = FALSE, abb.end = c("juvenile", "nonmacro", "unIDd"), sp.exc = "lost", subset.col = "Conf", subset.lev = "n", xlab = "Individual ID")
dev.off()
conf.mat$unconf <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Conf == "n"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Conf == "n"], levels = sp.idd))

# confusion matrix for large specimens 
# add a column with the cutoffs
IDs.long$Size <- NA
IDs.long$Size[IDs.long$MeanDia >= 200] <- "Large"
IDs.long$Size[IDs.long$MeanDia < 200] <- "Small"

png("F_Figures/confusion_large_key.png", 1000, 650)
conf_mat(IDs.long, "ID", "DefinitiveID", axes.same = FALSE, abb.end = c("juvenile", "nonmacro", "unIDd"), sp.exc = "lost", subset.col = "Size", subset.lev = "Large", xlab = "Individual ID")
dev.off()
conf.mat$large <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Size == "Large"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Size == "Large"], levels = sp.idd))

# confusion matrix for small specimens 
png("F_Figures/confusion_small_key.png", 900, 500)
conf_mat(IDs.long, "ID", "DefinitiveID", axes.same = FALSE, abb.end = c("juvenile", "nonmacro", "unIDd"), sp.exc = "lost", subset.col = "Size", subset.lev = "Small", xlab = "Individual ID")
dev.off()
conf.mat$small <- confusionMatrix(factor(IDs.long$DefinitiveID[IDs.long$Size == "Small"], levels = sp.idd), factor(IDs.long$ID[IDs.long$Size == "Small"], levels = sp.idd))

# confusion matrix summaries (esp. kappa)
# Accuracy is the fraction of correct identifications
# Accuracy null is (I think) the expected accuracy - which is the accuracy that would be expected based on random.
# the kappa statistic is (observed - expected) / (1 - expected)

conf.mat$sp$overall # kappa : 0.58
conf.mat$exp$overall # kappa : 0.78
conf.mat$stu$overall # kappa : 0.54
conf.mat$conf$overall # kappa : 0.76
conf.mat$unconf$overall # kappa : 0.21
conf.mat$large$overall # kappa : 0.64
conf.mat$small$overall # kappa : 0.38

