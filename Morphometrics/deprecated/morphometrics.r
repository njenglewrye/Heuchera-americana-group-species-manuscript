library(dplyr)

# Load dataset
combined <- read.csv("morphometric_data_final.csv")

combined$Species <- as.factor(combined$Species)
head(combined)

# Code numbers for variegation column
combined$Variegation <- gsub("None", 0, combined$Variegation)
combined$Variegation <- gsub("Low", 1, combined$Variegation)
combined$Variegation <- gsub("Medium", 2, combined$Variegation)
combined$Variegation <- gsub("High", 3, combined$Variegation)

combined$Variegation <- as.numeric(combined$Variegation)


combined$Notes <- NULL
combined$Thyrse_length <- NULL
combined$Pedicel_lengths <- as.numeric(combined$Pedicel_lengths) # Some sort of hidden nonnumeric character causes the column to be dropped in aggregation; this coerces to NA
combined$Variegation <- NULL
combined$Style_length <- NULL
#combined$Upper_lamina_hair <- NULL


# Save old dataset
combined.nonnormalized <- combined

# Normalize data matrix
combined <- rapply(combined, scale, c("numeric","integer"), how="replace")


# Aggregate -- combine per-sample replicate measurements by averaging so there is one data point per specimen
# Statistics below can be run on individual measurements by using the combined rather than combined.formatted data frame 
# But run this line first
# combined <- na.omit(combined)
# Initial comparisons demonstrated that jackknife success is greatly improved by aggregation
library(dplyr)
combined %>% group_by(Barcode, Species) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> combined.formatted
combined.formatted <- as.data.frame(combined.formatted)
combined.formatted <- na.omit(combined.formatted)

combined.formatted <- combined.formatted[combined.formatted$Species %in% c("AMERICANA", "CALYCOSA", "BREVIPETALA", "HISPIDA"), ]


# Build the discriminant
library(MASS)
#discriminant <- lda(Species ~ Variegation + Long_petiole_hair + Short_petiole_hair + Upper_lamina_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Style_length + Pedicel_lengths + Cymule_branches, data = combined.formatted, na.action="na.omit")
discriminant <- lda(Species ~ Long_petiole_hair + Short_petiole_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Pedicel_lengths + Cymule_branches, data = combined.formatted, na.action="na.omit")


# Classification success
#discriminant.jackknife <- lda(Species ~ Variegation + Long_petiole_hair + Short_petiole_hair + Upper_lamina_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Style_length + Pedicel_lengths + Cymule_branches, data = combined.formatted, na.action="na.omit", CV = TRUE)
discriminant.jackknife <- lda(Species ~  Long_petiole_hair + Short_petiole_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Pedicel_lengths + Cymule_branches, data = combined.formatted, na.action="na.omit", CV = TRUE)
ct <- table(combined.formatted$Species, discriminant.jackknife$class)
sum(diag(prop.table(ct)))

# Predict species by the discriminant function
discriminant.prediction <- predict(discriminant)

# Create dataframe for plotting
plotdata <- data.frame(type = combined.formatted$Species, lda = discriminant.prediction$x)

library(ggplot2)
ggplot(plotdata) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 2.5)

# Multivariate MANOVA
#res.man <- manova(cbind(Variegation, Long_petiole_hair, Short_petiole_hair, Upper_lamina_hair, Tooth_length, Tooth_Width, Flower_length, Hypanthium_long, Hypanthium_short, Flower_width, Stamen_length, Style_length, Pedicel_lengths, Cymule_branches) ~ Species, data = combined)
res.man <- manova(cbind(Long_petiole_hair, Short_petiole_hair, Tooth_length, Tooth_Width, Flower_length, Hypanthium_long, Hypanthium_short, Flower_width, Stamen_length, Pedicel_lengths, Cymule_branches) ~ Species, data = combined)
summary(res.man)

# Break down variable importance
summary.aov(res.man)

# Assess species pairwise significance
# You must drop perfectly correlated values or you will get a rank deficiency error
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#pairwise.adonis(combined.formatted[,c("Variegation", "Long_petiole_hair", "Short_petiole_hair", "Upper_lamina_hair", "Tooth_length", "Tooth_Width", "Flower_length", "Hypanthium_long", "Hypanthium_short", "Flower_width", "Stamen_length", "Style_length", "Pedicel_lengths", "Cymule_branches")], combined.formatted$Species, sim.method = "euclidean", p.adjust.m = "hochberg", perm = 10000)
pairwise.adonis(combined.formatted[,c("Long_petiole_hair", "Short_petiole_hair", "Tooth_length", "Tooth_Width", "Flower_length", "Hypanthium_long", "Hypanthium_short", "Flower_width", "Stamen_length", "Pedicel_lengths", "Cymule_branches")], combined.formatted$Species, sim.method = "euclidean", p.adjust.m = "hochberg", perm = 10000)

## Confusion matrix
#library(class)
## K-nearest neighbor classification
## This is a bit ad hoc as the training and test datasets are the same
#nn_classification <- knn(log(combined[,c("zygomorphylong","zygomorphyshort", "sinus", "ovary", "sinusovaryratio", "flowerlength", "branchlength", "thyrselength")]), log(combined[,c("zygomorphylong","zygomorphyshort", "sinus", "ovary", "sinusovaryratio", "flowerlength", "branchlength", "thyrselength")]), combined$species.x, k = 3)
#confusion_matrix = table(nn_classification, combined$species.x)
#confusion_matrix
#accuracy = sum(nn_classification == combined$species.x)/length(combined$species.x)
#accuracy

# Univariate boxplots

p1 <- ggplot(combined.nonnormalized, aes(y=Long_petiole_hair, x=Species)) + geom_boxplot() + xlab("") + ylab("Longest petiole hair length (mm)") + theme(axis.text.x = element_text(angle = 45))
p2 <- ggplot(combined.nonnormalized, aes(y=Tooth_length, x=Species)) + geom_boxplot() + xlab("") + ylab("Tooth length (mm)") + theme(axis.text.x = element_text(angle = 45))
p3 <- ggplot(combined.nonnormalized, aes(y=Flower_length, x=Species)) + geom_boxplot() + xlab("") + ylab("Flower length (mm)") + theme(axis.text.x = element_text(angle = 45))
p4 <- ggplot(combined.nonnormalized, aes(y=Stamen_length, x=Species)) + geom_boxplot() + xlab("") + ylab("Stamen length (mm)") + theme(axis.text.x = element_text(angle = 45))
p5 <- ggplot(combined.nonnormalized, aes(y=as.numeric(Pedicel_lengths), x=Species)) + geom_boxplot() + xlab("") + ylab("Pedicel length (mm)") + theme(axis.text.x = element_text(angle = 45))
p6 <- ggplot(combined.nonnormalized, aes(y=Hypanthium_long/Hypanthium_short, x=Species)) + geom_boxplot() + xlab("") + ylab("Zygomorphy ratio") + coord_cartesian(ylim = c(0,3)) + theme(axis.text.x = element_text(angle = 45)) #+ geom_text(aes(label=Barcode), hjust=0, vjust=0)


library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)


## PCA analysis
pca <- prcomp(combined.formatted[,c("Variegation", "Long_petiole_hair", "Short_petiole_hair", "Upper_lamina_hair", "Tooth_length", "Tooth_Width", "Flower_length", "Hypanthium_long", "Hypanthium_short", "Flower_width", "Stamen_length", "Style_length", "Pedicel_lengths", "Cymule_branches")], center = TRUE, scale = TRUE)
#
## Plot PCA
library(ggfortify)
## PCA with 90% confidence ellipses
autoplot(pca, data = combined.formatted, colour = 'Species', frame = TRUE, frame.type = 'norm', level = 0.9) # PCA plot
## PCA with plots of loadings
#autoplot(pca, data = combined, colour = 'species.x', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) # PCA plot with loadings



## Univariate ANOVA
#summary(aov(branchlength ~ species.x, data = combined))
#summary(aov(sinusovaryratio ~ species.x, data = combined))
#summary(aov(zygomorphyratio ~ species.x, data = combined))
#summary(aov(pedicellength ~ species.x, data = combined))
#summary(aov(flowerlength ~ species.x, data = combined))
#summary(aov(thyrselength ~ species.x, data = combined))