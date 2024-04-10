#Advanced bioinformatics assessment, m2305834

#task 3.1
sum(5:55)
#task 3.2
sumfun <- function(n) {sum(5:n)}
sumfun(10)
sumfun(20)
sumfun(100)
#task 3.3
a <- 1
b <- 1
fibonacci <- function(n) { cat(a, b, sep=", ")

for(i in 3:n) {next_fib <- a + b
cat(", ", next_fib, sep = "")
a <- b
b <- next_fib}
cat("\n")}
fibonacci(12)

#task 3.4
library(ggplot2)
ggplot(mtcars, aes(x = factor(gear), y = mpg, fill = factor(gear))) +
  geom_boxplot() +
  labs(x = "Number of Gears", y = "Miles per Gallon") +
  ggtitle("Boxplot of Miles per Gallon by Number of Gears") +
  scale_fill_discrete(name = "Number of Gears")

#task 3.5
data(cars)
lm(dist ~ speed, data = cars)
attr(cars$speed, "unit")
attr(cars$dist, "unit")


data(cars)
model <- lm(dist ~ speed, data = cars)
summary(model)
units_speed <- attr(cars$speed, "unit")
units_dist <- attr(cars$dist, "unit")
cat("Fitted Slope (Speed): ", coef(model)[2], "\n")
cat("Intercept: ", coef(model)[1], "\n")
cat("Standard Error of Slope: ", summary(model)$coefficients[2, "Std. Error"], "\n")
cat("Standard Error of Intercept: ", summary(model)$coefficients[1, "Std. Error"], "\n")
cat("Units for Speed: ", units_speed, "\n")
cat("Units for Distance: ", units_dist, "\n")

#task 3.6
library(ggplot2)
model <- lm(dist ~ speed, data = cars)
plot3.6 <- ggplot(cars, aes(x = speed, y = dist)) +
  geom_point() +  # Add the data points
  geom_smooth(method = "lm", se = FALSE, color = "purple") + 
  labs(x = "Speed (mph)", y = "Braking Distance (feet)") +
  ggtitle("Linear Relationship Between Speed and Braking Distance") + 
  theme_minimal()
print(plot3.6)

#task 3.7 

#linear regression to estimate breaking reaction time
reaction_time_estimate <- coef(model_reaction_time)[2]
print(paste("Estimated Average Reaction Time (seconds):", reaction_time_estimate))
reaction_time <- lm(dist ~ 0 + speed_squared + I(speed^2), data = cars)
summary(reaction_time)
#defining conversion for mph to ft/s
conversion_factor <- 1.46667
#reaction time estimation result
reaction_time_estimate <- coef(reaction_time)[1] / conversion_factor
print(paste("Estimated Average Reaction Time (seconds):", reaction_time_estimate))

#plotting regression for Breaking reaction time

library(ggplot2)

plot3.7 <- ggplot(cars, aes(x = speed, y = dist)) +
  geom_point() + 
  geom_smooth(method = "lm", formula = y ~ I(x^2), se = FALSE, color = "purple") +  
  labs(x = "Speed (ft/s)", y = "Braking Distance (feet)") +  
  ggtitle("Breaking Reaction Time") + 
  theme_minimal()  

print(plot3.7)

#task 3.8

setwd("/Users/HannahGross/Desktop/A. Bioinf/Day5/LMS_RNAseq_short-master-2023-final/course")

count_data <- read.csv("exercises/data/exercise1_counts.csv", header = T, row.names = 1)
head(count_data)

sample_description <- read.table("exercises/data/exercise1_sample_description.info", sep = "\t", header = TRUE )
head(sample_description)

#task 3.9

col_data <- data.frame(sample = sample_description$sample, condition = sample_description$condition, batch = sample_description$batch)
head(col_data)

all(colnames(count_data) == col_data$name)

#task 3.10

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design =~ condition)
dds

#task 3.11

rlog3.11 <- rlog(dds)
vsd3.11 <- varianceStabilizingTransformation(dds)
rlog3.11
vsd3.11

#task 3.12
library("pheatmap")
dds <- DESeq(dds)
dds_counts <- counts(dds, normalized = TRUE)
head(dds_counts)

select <- order(rowMeans(dds_counts), decreasing = TRUE)[1:40]
head(select)

pheatmap(assay(rlog3.11)[select, ])
pheatmap(assay(vsd3.11)[select, ])

#task 3.13
sample_dist <- dist(t(assay(rlog3.11)))
class(sample_dist)
sdm <- as.matrix(sample_dist)
class(sdm)

library("RColorBrewer")
colors <- colorRampPalette(rev(brewer.pal(9, "Purples")))(255)
rownames(sdm) <- rlog3.11$Group
colnames(sdm) <- NULL

pheatmap(sdm,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         col = colors)


#task 3.14

plotPCA(rlog3.11, intgroup = "condition")


#task 3.15

plotPCA(vsd3.11, intgroup = "condition")


#3.16 
#task 3.16
setwd("/Users/HannahGross/Desktop/A. Bioinf/Day5/LMS_ChIPseq_short-master-2023-final/course")

library(GenomicRanges)

melPeak_Rep1 <- read.delim("data/MacsPeaks/mycmelrep1_peaks.xls",sep="\t",comment.char = "#")
melPeak_Rep2 <- read.delim("data/MacsPeaks/mycmelrep2_peaks.xls",sep="\t",comment.char = "#")
melRep1_GR <- GRanges(
  seqnames=melPeak_Rep1[,"chr"],
  IRanges(melPeak_Rep1[,"start"],
          melPeak_Rep1[,"end"]
  )
)

mcols(melRep1_GR) <- melPeak_Rep1[,c("abs_summit", "fold_enrichment")]

melRep2_GR <- GRanges(
  seqnames=melPeak_Rep2[,"chr"],
  IRanges(melPeak_Rep2[,"start"],
          melPeak_Rep2[,"end"]
  )
)

mcols(melRep2_GR) <- melPeak_Rep2[,c("abs_summit", "fold_enrichment")]

table(melRep1_GR %over% melRep2_GR)
length(melRep1_GR[melRep1_GR %over% melRep2_GR])
table(!melRep1_GR %over% melRep2_GR)
length(melRep1_GR[!melRep1_GR %over% melRep2_GR])
commonMelPeaks <- melRep1_GR[melRep1_GR %over% melRep2_GR]
head(commonMelPeaks)
melRep1_GRSummits <- melRep1_GR
start(melRep1_GRSummits) <- end(melRep1_GRSummits) <- melRep1_GR$abs_summit
melRep1_GRSummits

#3.17
sorted_commonMelPeaks <- melRep1_GRSummits[order(mcols(melRep1_GRSummits)$fold_enrichment, decreasing = TRUE)]
top_500_peaks <- sorted_commonMelPeaks[1:500]
center <- as.integer(width(top_500_peaks) / 2)
flank <- 100 
resized_peaks <- resize(top_500_peaks, width = 2*flank)
head(resized_peaks)

#3.18 
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm9)

genome <- BSgenome.Mmusculus.UCSC.mm9

resized_peaks_gr <- GRanges(seqnames=resized_peaks$seqnames,
                            ranges=IRanges(start=resized_peaks$start - flank,
                                           end=resized_peaks$start + flank))

sequences <- getSeq(genome, resized_peaks_gr)

writeXStringSet(sequences, file="top_500_peaks.fa")





