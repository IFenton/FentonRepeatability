# Repeatability analysis
# Project: Repeatability
# Author: Isabel Fenton
# Date created: 9/3/2017
# Date last edited: 16/3/2017
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

people.df <- data.frame(Name = people)
str(people.df)

# 2b. Species level dataframe ---------------------------------------------
species <- grep(" ", checklist$Species, value = TRUE)
genera <- unique(grep("^[[:upper:]]", checklist$ChecklistGenus, value = TRUE))
  
species.df <- data.frame(Species = checklist$Species)


# 3. Person level statistics -----------------------------

# Number of species identified
people.df$NumSpecID <- apply(completeIDs[, people.col], 2, function(x) sum(unique(x) %in% species))

# Number of genera identified
unique(grep("^[[:upper:]]", unlist(strsplit(as.character(completeIDs$DefinitiveID), " ")), value = TRUE))
people.df$NumGenID <- apply(completeIDs[, people.col], 2, function(x) length(unique(grep("^[[:upper:]]", unlist(strsplit(as.character(x), " ")), value = TRUE))))

# Number of specimens lost
sum(completeIDs$DefinitiveID == "lost")
people.df$numLost <- apply(completeIDs[, people.col], 2, function(x) sum(x == "lost"))
people.df$TotalInd <- apply(completeIDs[, people.col], 2, function(x) sum(x != "lost"))

# Fraction of species identified
people.df$numUnIDd <- apply(completeIDs[, people.col], 2, function(x) sum(x == "unIDd"))

people.df$fracSpecIDd <- people.df$NumSpecID / people.df$NumSpecID[people.df$Name == "Definitive"]

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

# 4. Percentage accuracy per specimen ------------------------------------

# 5. Create a confusion matrix --------------------------------------------


# 6. Dendrogram ---------------------------------------------------------
# n.b. use daisy rather than hist, as the data are factors
png("Figures/dendrogram.png")
plot(hclust(daisy(data.frame(t(completeIDs[, people.col])))))
dev.off()

