#Author: Eleonor Cayford

#Set working directory 
setwd("C:/Users/ecayf_000/PycharmProjects/IndependentProject")
#Opens the necessary libraries
library(stargazer)
library(sciplot)

# Reading in files---------------------------------------------------------------
#Reads in statistics on the files and removing the extra last column
rawStats <- read.csv('raw_data_stats.csv')
rawStats <- rawStats[,1:4]
processedStats <- read.csv('processed_data_stats.csv')
processedStats <- processedStats[,1:4]
#Reads in the quality scores 
rawAverageQualityScores <- read.csv("raw_data_quality_scores.csv", check.names = FALSE)
processedAverageQualityScores <- read.csv("processed_data_quality_scores.csv", check.names = FALSE)
#Reads in the base content by read while converting it to a transposed matrix
#This allows for easier plotting later on
rawBaseContent <- t(as.matrix(read.csv("raw_data_base_content.csv")))
processedBaseContent <- t(as.matrix(read.csv("processed_data_base_content.csv")))

#Generating table-------------------------------------------------------------------
row.names(rawStats) <- "Raw_Data"
row.names(processedStats) <- "Processed_Data"
totalStatistics <- rbind(rawStats, processedStats)
#Generates a table that gives the statistics 
stargazer(totalStatistics, title="Statistics", type="html", digits=3, out="StatisticsTable.html", summary = FALSE)

#Generating quality plots----------------------------------------------------------
#Creates a Fastq quality plot: a display of quality score across all bases
par(mfrow=c(2,1))
boxplot.default(rawAverageQualityScores, main = " Raw data fastq quality plot", xlab = "Location in read", ylab = "Quality score", las = 1, col=("yellow"), cex=.3)
boxplot.default(processedAverageQualityScores, main = " Processed data fastq quality plot", xlab = "Location in read", ylab = "Quality score", las = 1, col=("yellow"), cex=.3)

#Generating base content plots------------------------------------------------------
par(mfrow=c(2,1))
#Creates a barplot showing the percentage of different bases at each read location
barplot(rawBaseContent, main="Percentage of each base by location on raw reads", xlab="Location on read", ylab = "Percentage of total bases", col=c("blue","red", "yellow", "green", "black"), legend = rownames(rawBaseContent), las=1)
axis(1, at = c(0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120), las=1, labels=c('0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'))
barplot(processedBaseContent, main="Percentage of each base by location on processed reads", xlab="Location on read", ylab = "Percentage of total bases", col=c("blue","red", "yellow", "green", "black"), legend = rownames(processedBaseContent), las=1, cex=.4)
axis(1, at = c(0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120), las=1, labels=c('0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'))
