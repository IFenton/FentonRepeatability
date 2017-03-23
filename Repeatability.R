# Repeatability analysis
# Project: Repeatability
# Author: Isabel Fenton
# Date created: 9/3/2017
# Date last edited: 23/3/2017
# 
# Code for analysis of the repeatability results. For more details see Lab book Repeatability.docx
# 
# Previous file: size histogram.R
# Next file:

# Inputs ------------------------------------------------------------------
# CompleteIDs.csv
# Checklist.csv

# Outputs -----------------------------------------------------------------

# Source files / libraries ------------------------------------------------
library(cluster) # for the dendrogram
library(caret) # for the confusion matrix

rm(list = ls())

# 1. Load in the data -----------------------------------------------------
checklist <- read.csv("Data/Checklist.csv")
str(checklist)

completeIDs <- read.csv("Data/CompleteIDs.csv")
str(completeIDs)

# set the levels of each ID column to the checklist
for (i in grep("ID", names(completeIDs), value = TRUE)) {
  completeIDs[,i] <- factor(completeIDs[,i], levels = unique(c(levels(checklist$Species), levels(checklist$ChecklistGenus))))
}
rm(i)
str(completeIDs)

# 2. Create dataframes for storing the results -------------------------------

# 2a. Person level dataframe ----------------------------------------------
# list of people in the analysis
people <- gsub("ID", "", grep("ID", names(completeIDs), value = TRUE))
people.col <- grep("ID", names(completeIDs), value = TRUE)
p.conf.col <- grep("C$", names(completeIDs), value = TRUE)

people.df <- data.frame(Name = people)
str(people.df)

# 2b. Species level dataframe ---------------------------------------------
species <- grep(" ", checklist$Species, value = TRUE)
genera <- unique(grep("^[[:upper:]]", checklist$ChecklistGenus, value = TRUE))
  
species.df <- data.frame(Species = checklist$Species)


# 3. Person level statistics -----------------------------

# 3a. All individuals -----------------------------------------------------
# Number of species identified
people.df$NumSpID <- apply(completeIDs[, people.col], 2, function(x) sum(unique(x) %in% species))

# Number of genera identified
unique(grep("^[[:upper:]]", unlist(strsplit(as.character(completeIDs$DefinitiveID), " ")), value = TRUE))
people.df$NumGenID <- apply(completeIDs[, people.col], 2, function(x) length(unique(grep("^[[:upper:]]", unlist(strsplit(as.character(x), " ")), value = TRUE))))

# Number of specimens lost
sum(completeIDs$PooleID != "lost")
sum(completeIDs$PooleID != "lost"& completeIDs$PooleID %in% checklist$Species)
people.df$numLost <- apply(completeIDs[, people.col], 2, function(x) sum(x == "lost"))
people.df$TotalIndSp <- apply(completeIDs[, people.col], 2, function(x) sum(x != "lost" & x %in% checklist$Species))
people.df$TotalIndGe <- apply(completeIDs[, people.col], 2, function(x) sum(x != "lost"))

# Fraction of species identified
people.df$numUnIDd <- apply(completeIDs[, people.col], 2, function(x) sum(x == "unIDd"))

people.df$fracSpecIDd <- people.df$NumSpID / people.df$NumSpID[people.df$Name == "Definitive"]

# Fraction of genera identified
people.df$fracGenIDd <- people.df$NumGenID / people.df$NumGenID[people.df$Name == "Definitive"]


# Percentage accuracy per person
sum(completeIDs[,grep("Fenton", names(completeIDs))[1]] == completeIDs$DefinitiveID) / sum(completeIDs[,grep("Fenton", names(completeIDs))[1]] != "lost")
sum(completeIDs$DefinitiveID == completeIDs$DefinitiveID) / people.df$TotalInd[people.df$Name == "Definitive"]

# species
people.df$fracCorrSp <- apply(completeIDs[, people.col], 2, function(x) sum(x == completeIDs$DefinitiveID) / sum(x != "lost"))
mean(people.df$fracCorrSp)
median(people.df$fracCorrSp)
median(people.df$fracCorrSp[7:nrow(people.df)])
median(people.df$fracCorrSp[3:6])

# genus
gsub(" .*", "", completeIDs$DefinitiveID)

sum(gsub(" .*", "", completeIDs$FentonID) == gsub(" .*", "", completeIDs$DefinitiveID)) / sum(completeIDs[,grep("Fenton", names(completeIDs))[1]] != "lost")

people.df$fracCorrGen <- apply(completeIDs[, people.col], 2, function(x) sum(gsub(" .*", "", x) == gsub(" .*", "", completeIDs$DefinitiveID)) / sum(x != "lost"))
mean(people.df$fracCorrGen)

# 3b. Split by confidence -------------------------------------------------
# number of identifications by confidence (working at species level)
sum(completeIDs$FentonC == "y", na.rm = TRUE)

people.df$numCy <- NA
people.df$numCy[3:nrow(people.df)] <- apply(completeIDs[, p.conf.col], 2, function(x) sum(x == "y", na.rm = TRUE))

people.df$numCm <- NA
people.df$numCm[3:nrow(people.df)] <- apply(completeIDs[, p.conf.col], 2, function(x) sum(x == "m", na.rm = TRUE))

people.df$numCn <- NA
people.df$numCn[3:nrow(people.df)] <- apply(completeIDs[, p.conf.col], 2, function(x) sum(x == "n", na.rm = TRUE))

people.df

# fraction of confident identifications
sum(completeIDs[, people.col[2+1]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[1]] == "y", na.rm = TRUE) / sum(completeIDs[, p.conf.col[1]] == "y", na.rm = TRUE)

sum(completeIDs[, people.col[2+1]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[1]] == "m", na.rm = TRUE) / sum(completeIDs[, p.conf.col[1]] == "m", na.rm = TRUE)

sum(completeIDs[, people.col[2+1]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[1]] == "n", na.rm = TRUE) / sum(completeIDs[, p.conf.col[1]] == "n", na.rm = TRUE)

people.df$fracCorrY <- NA
people.df$fracCorrM <- NA
people.df$fracCorrN <- NA
for (i in 1:length(p.conf.col)){
  people.df$fracCorrY[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[i]] == "y", na.rm = TRUE) / sum(completeIDs[, p.conf.col[i]] == "y", na.rm = TRUE)
  people.df$fracCorrM[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[i]] == "m", na.rm = TRUE) / sum(completeIDs[, p.conf.col[i]] == "m", na.rm = TRUE)
  people.df$fracCorrN[2 + i] <- sum(completeIDs[, people.col[2+i]] == completeIDs$DefinitiveID & completeIDs[, p.conf.col[i]] == "n", na.rm = TRUE) / sum(completeIDs[, p.conf.col[i]] == "n", na.rm = TRUE)
}

hist(people.df$fracCorrY)
hist(people.df$fracCorrM)
hist(people.df$fracCorrN)

# 4. Percentage accuracy per specimen ------------------------------------
sum(completeIDs[2, people.col[2:length(people.col)]] == as.character(completeIDs$DefinitiveID[2]) & completeIDs[2, people.col[2:length(people.col)]] != "lost") / sum(completeIDs[2, people.col[2:length(people.col)]] != "lost")

completeIDs$fracCorr <- apply(completeIDs[, people.col[2:length(people.col)]], 1, function(x) sum(x == x[1] & x != "lost") / sum(x != "lost"))

hist(completeIDs$fracCorr)

# 5. Create a confusion matrix --------------------------------------------
image(confusionMatrix(completeIDs$DefinitiveID, completeIDs$PurvisID)$table, axes=FALSE)

image(confusionMatrix(rep(completeIDs$DefinitiveID,2), factor(c(as.character(completeIDs$FentonID), as.character(completeIDs$PurvisID)), levels = levels(completeIDs$DefinitiveID)))$table, axes=FALSE)


# 6. Dendrogram ---------------------------------------------------------
# n.b. use daisy rather than hist, as the data are factors
png("Figures/dendrogram.png")
plot(hclust(daisy(data.frame(t(completeIDs[, people.col])))))
dev.off()



par(mfrow = c(2,1))
plot(factor(completeIDs$DefinitiveID), completeIDs$fracCorr, las = 2)