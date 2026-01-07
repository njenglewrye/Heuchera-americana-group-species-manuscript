#Calculating population genetics statistics and indices in R

###Converting my file to hierfstat format

library(adegenet)
library(vcfR)
library(hierfstat)
library(dplyr)

americanaV <- read.vcfR("americana.SNPall.pruned.vcf.recode.vcf")

#I need a genind object first
americanaGenind <- vcfR2genind(americanaV)

#Assign populations to the strata slot of my genlight object

americanapopdata <- read.table("pca_popfile.txt")
americanapopdata <- as.data.frame(americanapopdata)

strata(americanaGenind) <- americanapopdata

#Setting populations to the pop slot
# Create a factor vector for population assignments
americanapopfactor <- factor(americanapopdata$V2)

# Assign the factor to the 'pop' slot of the genlight object
pop(americanaGenind) <- americanapopfactor

# Convert genind to hierfstat format
americanahs <- genind2hierfstat(americanaGenind)


####################################
## Across taxon Fst + heatmap
####################################


#############################################Calculate Fixation index (FST) according to Nei (1987)

NeiFST <- pairwise.neifst(americanahs, diploid = TRUE)

#############################################Calculate Fixation index (FST) according to Weir and Cockerham (1984)

americanaWCFST <- pairwise.WCfst(americanahs, diploid = TRUE)
write.csv(americanaWCFST, "americanapairwiseWCfst.csv")  



#plot heatmap of pairwise Fst
americanapairwiseWCfst <- read.csv("americanapairwiseWCfst.csv", header=TRUE, row.names = 1)
americanapairwiseWCfst <-as.matrix(americanapairwiseWCfst)

heatmap(americanapairwiseWCfst, margins=c(10,10), revC=T, col=rev(heat.colors(256)))

heatmap(americanapairwiseWCfst, margins=c(10,10), revC=T, col=rev(heat.colors(400)))



library(pheatmap)

pheatmap(americanapairwiseWCfst,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("white", "navy", "firebrick3"))(100))

pheatmap(americanapairwiseWCfst,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("white", "yellow", "red"))(256))

pheatmap(americanapairwiseWCfst,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("white", "turquoise", "orange"))(100))
         
         
##########################################
## Across taxon assorted other stats
##########################################


#############################################Calculate Allelic richness per population

americanaAr <- allelic.richness(americanahs)$Ar
americanameanAr <-colMeans(americanaAr, na.rm = TRUE) 
americanameanAr

#AlRi <- allelic.richness(gen_data)
#print(AlRi$min.all)
#print(AlRi$Ar)

#############################################Calculate Observed heterozygosity per population (Ho)
americanahet_obs <- basic.stats(americanahs)$Ho 
americanameanHo <-colMeans(americanahet_obs, na.rm = TRUE) 

#############################################Calculate within population gene diversity (Hs)
americanagendiv <- basic.stats(americanahs)$Hs
americanameanHs <-colMeans(americanagendiv, na.rm = TRUE) 

#############################################Calculate Wright's inbreeding coefficient (FIS)
americanaInbrcoef <- basic.stats(americanahs)$Fis
americanameanFis <-colMeans(americanaInbrcoef, na.rm = TRUE) 
