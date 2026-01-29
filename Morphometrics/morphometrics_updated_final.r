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

########################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################
############################################################## Trying to add all taxa below #######################################################################################################################################################
########################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################

# Aggregate -- combine per-sample replicate measurements by averaging so there is one data point per specimen
# Statistics below can be run on individual measurements by using the combined rather than combined.formatted data frame 
# But run this line first
# combined <- na.omit(combined)
# Initial comparisons demonstrated that jackknife success is greatly improved by aggregation
library(dplyr)
combined %>% group_by(Barcode, Species) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> combined.formatted
combined.formatted <- as.data.frame(combined.formatted)
combined.formatted <- na.omit(combined.formatted)

combined.formatted <- combined.formatted[combined.formatted$Species %in% c("AMERICANA", "CALYCOSA", "BREVIPETALA", "HISPIDA", "GRAYANA", "RICHARDSONII", "HIRSUTICAULIS", "FUMOSIMONTANA"), ]


# Build the discriminant
library(MASS)
#discriminant <- lda(Species ~ Variegation + Long_petiole_hair + Short_petiole_hair + Upper_lamina_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Style_length + Pedicel_lengths + Cymule_branches, data = combined.formatted, na.action="na.omit")
discriminant <- lda(Species ~ Long_petiole_hair + Short_petiole_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Pedicel_lengths + Cymule_branches, data = combined.formatted, na.action="na.omit")


# Classification success
#discriminant.jackknife <- lda(Species ~ Variegation + Long_petiole_hair + Short_petiole_hair + Upper_lamina_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Style_length + Pedicel_lengths + Cymule_branches, data = combined.formatted, na.action="na.omit", CV = TRUE)
discriminant.jackknife <- lda(Species ~  Long_petiole_hair + Short_petiole_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Pedicel_lengths + Cymule_branches, data = combined.formatted, na.action="na.omit", CV = TRUE)
ct <- table(combined.formatted$Species, discriminant.jackknife$class)
sum(diag(prop.table(ct)))
# 0.7032967

# Predict species by the discriminant function
discriminant.prediction <- predict(discriminant)

# Create dataframe for plotting
plotdata <- data.frame(type = combined.formatted$Species, barcode = combined.formatted$Barcode, lda = discriminant.prediction$x)

#Pre fumosimontana
#library(ggplot2)
#ggplot(plotdata) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 2.5)

library(ggplot2)

ggplot(plotdata) +
  geom_point(
    aes(lda.LD1, lda.LD2, colour = type),
    size = 2.5
  ) +
  #geom_text(
  #  aes(lda.LD1, lda.LD2, label = barcode),
  #  vjust = -0.5
  #) +
  scale_colour_manual(
    values = c(
      "AMERICANA"      = "#E69F00",
      "BREVIPETALA"    = "#56B4E9",
      "CALYCOSA"       = "#009E73",
      "FUMOSIMONTANA"  = "#0072B2",
      "GRAYANA"        = "#F0E442",
      "HIRSUTICAULIS"  = "#D55E00",
      "HISPIDA"        = "#CC79A7",
      "RICHARDSONII"   = "#999999"
    )
  )



### PCA analysis (All taxa [removed "Variegation", "Style_length"] )
#pca <- prcomp(combined.formatted[,c("Long_petiole_hair", "Short_petiole_hair", "Upper_lamina_hair", "Tooth_length", "Tooth_Width", "Flower_length", "Hypanthium_long", "Hypanthium_short", "Flower_width", "Stamen_length", "Pedicel_lengths", "Cymule_branches")], center = TRUE, scale = TRUE)
##
### Plot PCA
#library(ggfortify)
### PCA with 90% confidence ellipses
#autoplot(pca, data = combined.formatted, colour = 'Species', frame = TRUE, frame.type = 'norm', level = 0.9) # PCA plot
### PCA with plots of loadings
##autoplot(pca, data = combined, colour = 'species.x', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) # PCA plot with loadings




########################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################
############################################################## Trying to add all taxa above #######################################################################################################################################################
########################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################


############################################################## Downstream results are subset #########################################################################################################################################################

# Aggregate -- combine per-sample replicate measurements by averaging so there is one data point per specimen
# Statistics below can be run on individual measurements by using the combined rather than combined.formatted data frame 
# But run this line first
# combined <- na.omit(combined)
# Initial comparisons demonstrated that jackknife success is greatly improved by aggregation
library(dplyr)
combined %>% group_by(Barcode, Species) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> combined.formatted
combined.formatted <- as.data.frame(combined.formatted)
combined.formatted <- na.omit(combined.formatted)

combined.formatted.subset <- combined.formatted[combined.formatted$Species %in% c("AMERICANA", "CALYCOSA", "BREVIPETALA", "HISPIDA"), ]


# Build the discriminant
library(MASS)
#discriminant <- lda(Species ~ Variegation + Long_petiole_hair + Short_petiole_hair + Upper_lamina_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Style_length + Pedicel_lengths + Cymule_branches, data = combined.formatted, na.action="na.omit")
discriminant <- lda(Species ~ Long_petiole_hair + Short_petiole_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Pedicel_lengths + Cymule_branches, data = combined.formatted.subset, na.action="na.omit")


# Classification success
#discriminant.jackknife <- lda(Species ~ Variegation + Long_petiole_hair + Short_petiole_hair + Upper_lamina_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Style_length + Pedicel_lengths + Cymule_branches, data = combined.formatted, na.action="na.omit", CV = TRUE)
discriminant.jackknife <- lda(Species ~  Long_petiole_hair + Short_petiole_hair + Tooth_length + Tooth_Width + Flower_length + Hypanthium_long + Hypanthium_short + Flower_width + Stamen_length + Pedicel_lengths + Cymule_branches, data = combined.formatted.subset, na.action="na.omit", CV = TRUE)
ct <- table(combined.formatted.subset$Species, discriminant.jackknife$class)
sum(diag(prop.table(ct)))
# 0.6666667

# Predict species by the discriminant function
discriminant.prediction <- predict(discriminant)

# Create dataframe for plotting
plotdata <- data.frame(type = combined.formatted.subset$Species, lda = discriminant.prediction$x)

library(ggplot2)

ggplot(plotdata) +
  geom_point(
    aes(x = lda.LD1, y = lda.LD2, colour = type),
    size = 2.5
  ) +
  scale_colour_manual(
    values = c(
      "AMERICANA"   = "#E69F00",
      "BREVIPETALA" = "#56B4E9",
      "CALYCOSA"    = "#009E73",
      "HISPIDA"     = "#CC79A7"
    )
  )
# theme_classic() # uncomment if you want all white background

#library(ggplot2)
#ggplot(plotdata) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 2.5)



# Multivariate MANOVA (all taxa)
#res.man <- manova(cbind(Variegation, Long_petiole_hair, Short_petiole_hair, Upper_lamina_hair, Tooth_length, Tooth_Width, Flower_length, Hypanthium_long, Hypanthium_short, Flower_width, Stamen_length, Style_length, Pedicel_lengths, Cymule_branches) ~ Species, data = combined)
res.man <- manova(cbind(Long_petiole_hair, Short_petiole_hair, Tooth_length, Tooth_Width, Flower_length, Hypanthium_long, Hypanthium_short, Flower_width, Stamen_length, Pedicel_lengths, Cymule_branches) ~ Species, data = combined)
summary(res.man)

# Df Pillai approx F num Df den Df    Pr(>F)    
# Species     7 2.6341   14.534     77   1855 < 2.2e-16 ***
#   Residuals 269                                            
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Break down variable importance
summary.aov(res.man)

# Response 1 :
#   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Species       7 247.977  35.425   140.6 < 2.2e-16 ***
#   Residuals   269  67.775   0.252                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response 2 :
#   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Species       7 243.230  34.747  117.08 < 2.2e-16 ***
#   Residuals   269  79.838   0.297                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response 3 :
#   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Species       7  44.747  6.3925  7.4728 3.201e-08 ***
#   Residuals   269 230.111  0.8554                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response 4 :
#   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Species       7  90.992 12.9989  18.086 < 2.2e-16 ***
#   Residuals   269 193.338  0.7187                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response 5 :
#   Df Sum Sq Mean Sq F value    Pr(>F)    
# Species       7 162.08 23.1547  49.001 < 2.2e-16 ***
#   Residuals   269 127.11  0.4725                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response 6 :
#   Df Sum Sq Mean Sq F value    Pr(>F)    
# Species       7 173.35 24.7644   59.26 < 2.2e-16 ***
#   Residuals   269 112.41  0.4179                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response 7 :
#   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Species       7  98.955 14.1364  20.993 < 2.2e-16 ***
#   Residuals   269 181.143  0.6734                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response 8 :
#   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Species       7  68.847  9.8353   12.21 1.509e-13 ***
#   Residuals   269 216.685  0.8055                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response 9 :
#   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Species       7  89.265 12.7522  22.461 < 2.2e-16 ***
#   Residuals   269 152.721  0.5677                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response 10 :
#   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Species       7  42.417  6.0596  6.8627 1.632e-07 ***
#   Residuals   269 237.520  0.8830                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Response 11 :
#   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Species       7  27.119  3.8742   3.962 0.0003871 ***
#   Residuals   269 263.039  0.9778                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 174 observations deleted due to missingness


# Assess species pairwise significance (all taxa)
# You must drop perfectly correlated values or you will get a rank deficiency error
# At this point, three variables were dropped from the final publication analysis due to species monomorphy or rank deficiency. The commented lines document the full variable set.
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#pairwise.adonis(combined.formatted[,c("Variegation", "Long_petiole_hair", "Short_petiole_hair", "Upper_lamina_hair", "Tooth_length", "Tooth_Width", "Flower_length", "Hypanthium_long", "Hypanthium_short", "Flower_width", "Stamen_length", "Style_length", "Pedicel_lengths", "Cymule_branches")], combined.formatted$Species, sim.method = "euclidean", p.adjust.m = "hochberg", perm = 10000)
pairwise.adonis(combined.formatted[,c("Long_petiole_hair", "Short_petiole_hair", "Tooth_length", "Tooth_Width", "Flower_length", "Hypanthium_long", "Hypanthium_short", "Flower_width", "Stamen_length", "Pedicel_lengths", "Cymule_branches")], combined.formatted$Species, sim.method = "euclidean", p.adjust.m = "hochberg", perm = 10000)



# Assess species pairwise significance (subset americana group)
# You must drop perfectly correlated values or you will get a rank deficiency error
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#pairwise.adonis(combined.formatted.subset[,c("Variegation", "Long_petiole_hair", "Short_petiole_hair", "Upper_lamina_hair", "Tooth_length", "Tooth_Width", "Flower_length", "Hypanthium_long", "Hypanthium_short", "Flower_width", "Stamen_length", "Style_length", "Pedicel_lengths", "Cymule_branches")], combined.formatted.subset$Species, sim.method = "euclidean", p.adjust.m = "hochberg", perm = 10000)
pairwise.adonis(combined.formatted.subset[,c("Long_petiole_hair", "Short_petiole_hair", "Tooth_length", "Tooth_Width", "Flower_length", "Hypanthium_long", "Hypanthium_short", "Flower_width", "Stamen_length", "Pedicel_lengths", "Cymule_branches")], combined.formatted.subset$Species, sim.method = "euclidean", p.adjust.m = "hochberg", perm = 10000)



## Confusion matrix
#library(class)
## K-nearest neighbor classification
## This is a bit ad hoc as the training and test datasets are the same
#nn_classification <- knn(log(combined[,c("zygomorphylong","zygomorphyshort", "sinus", "ovary", "sinusovaryratio", "flowerlength", "branchlength", "thyrselength")]), log(combined[,c("zygomorphylong","zygomorphyshort", "sinus", "ovary", "sinusovaryratio", "flowerlength", "branchlength", "thyrselength")]), combined$species.x, k = 3)
#confusion_matrix = table(nn_classification, combined$species.x)
#confusion_matrix
#accuracy = sum(nn_classification == combined$species.x)/length(combined$species.x)
#accuracy

# Univariate boxplots (all taxa)

p1 <- ggplot(combined.nonnormalized, aes(y=Long_petiole_hair, x=Species)) + geom_boxplot() + xlab("") + ylab("Longest petiole hair length (mm)") + theme(axis.text.x = element_text(angle = 45))
p2 <- ggplot(combined.nonnormalized, aes(y=Tooth_length, x=Species)) + geom_boxplot() + xlab("") + ylab("Tooth length (mm)") + theme(axis.text.x = element_text(angle = 45))
p3 <- ggplot(combined.nonnormalized, aes(y=Upper_lamina_hair, x=Species)) + geom_boxplot() + xlab("") + ylab("Upper laminar hair length (mm)") + theme(axis.text.x = element_text(angle = 45)) # This one was hard to measure
p4 <- ggplot(combined.nonnormalized, aes(y=Flower_length, x=Species)) + geom_boxplot() + xlab("") + ylab("Flower length (mm)") + theme(axis.text.x = element_text(angle = 45))
p5 <- ggplot(combined.nonnormalized, aes(y=Stamen_length, x=Species)) + geom_boxplot() + xlab("") + ylab("Stamen length (mm)") + theme(axis.text.x = element_text(angle = 45))
p6 <- ggplot(combined.nonnormalized, aes(y=as.numeric(Pedicel_lengths), x=Species)) + geom_boxplot() + xlab("") + ylab("Pedicel length (mm)") + theme(axis.text.x = element_text(angle = 45))
p7 <- ggplot(combined.nonnormalized, aes(y=Hypanthium_long/Hypanthium_short, x=Species)) + geom_boxplot() + xlab("") + ylab("Zygomorphy ratio") + coord_cartesian(ylim = c(0,3)) + theme(axis.text.x = element_text(angle = 45)) #+ geom_text(aes(label=Barcode), hjust=0, vjust=0)
p8 <- ggplot(combined.nonnormalized, aes(y=Flower_length/Flower_width, x=Species)) + geom_boxplot() + xlab("") + ylab("Flower length/width ratio") + coord_cartesian(ylim = c(0,3)) + theme(axis.text.x = element_text(angle = 45)) #+ geom_text(aes(label=Barcode), hjust=0, vjust=0)


library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2)

### PCA analysis (subset americana group [removed "Variegation", "Style_length"] )
#pca <- prcomp(combined.formatted.subset[,c("Long_petiole_hair", "Short_petiole_hair", "Upper_lamina_hair", "Tooth_length", "Tooth_Width", "Flower_length", "Hypanthium_long", "Hypanthium_short", "Flower_width", "Stamen_length", "Pedicel_lengths", "Cymule_branches")], center = TRUE, scale = TRUE)
##
### Plot PCA
#library(ggfortify)
### PCA of subset americana group with 90% confidence ellipses
#autoplot(pca, data = combined.formatted.subset, colour = 'Species', frame = TRUE, frame.type = 'norm', level = 0.9) # PCA plot
### PCA with plots of loadings
##autoplot(pca, data = combined.formatted.subset, colour = 'species.x', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) # PCA plot with loadings



## Univariate ANOVA
#summary(aov(branchlength ~ species.x, data = combined))
#summary(aov(sinusovaryratio ~ species.x, data = combined))
#summary(aov(zygomorphyratio ~ species.x, data = combined))
#summary(aov(pedicellength ~ species.x, data = combined))
#summary(aov(flowerlength ~ species.x, data = combined))
#summary(aov(thyrselength ~ species.x, data = combined))