#####Exploratory PCA in R
#You will need the (1) vcf file and (2) txt file with a sample id column and an assigned population name column next to it

library(adegenet)
library(vcfR)

######Example for single species H. americana
#convert vcf to genlight object
americanaGenlight <- read.vcfR("americana.SNPall.pruned.vcf.recode.vcf")
americanaGenlight <- vcfR2genlight(americanaGenlight) 

americanapca <- glPca(americanaGenlight)

#Assign populations to the strata slot of my genlight object

americanapopdata <- read.table("pca_popfile.txt")
americanapopdata <- as.data.frame(americanapopdata)

strata(americanaGenlight) <- americanapopdata

#Setting populations to the pop slot
# Create a factor vector for population assignments
americanapopfactor <- factor(americanapopdata$V2)

# Assign the factor to the 'pop' slot of the genlight object
pop(americanaGenlight) <- americanapopfactor


#individual names
scatter(americanapca)

library(ggplot2)

americanapca_coords <- as.data.frame(americanapca$scores)
americanapca_coords$group <- pop(americanaGenlight)

cb_palette <- c(
  "ACEROIDES"      = "black",
  "ALBA"           = "black",
  "AMERICANA"      = "#E69F00",
  "BREVIPETALA"    = "#56B4E9",
  "CALYCOSA"       = "#009E73",
  "CAROLINIANA"    = "black",
  "GRAYANA"        = "#F0E442",
  "HETERADENIA"    = "#0072B2",
  "HIRSUTICAULIS"  = "#D55E00",
  "HISPIDA"        = "#CC79A7",
  "LONGIFLORA"     = "black",
  "PUBESCENS"      = "black",
  "RICHARDSONII"   = "#999999"
)

#color
ggplot(americanapca_coords, aes(x = PC1, y = PC2, fill = group)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = cb_palette) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))

  
#####
# Stopped here, rest of lines not used for figures
#####
  

#access percent explaining variance
# Access eigenvalues
americanaeigenvalues <- americanapca$eig

# Calculate proportion of variance for each PC
americanaproportion_variance <- americanaeigenvalues / sum(americanaeigenvalues) * 100
americanaproportion_variance[1:2] #only for the first 2 PCs

####Discriminant Analysis of Principal Components (optional)

#DAPC
#testing optimum number of PCs (try n.pca up to 50, n.da=k-1)
americanadapc <- dapc(americanaGenlight, n.pca=50, n.da=13)
americanatest <- optim.a.score(americanadapc)

#then run with the optimum number
americanadapc <- dapc(americanaGenlight, n.pca=1, n.da=1)


#Plot distribution
par(cex=0.7)
scatter(americanadapc, scree.da=FALSE, bg="white", posi.leg="topleft", legend=TRUE,)

#Composition plot
devtools::install_github("zkamvar/ggcompoplot")

#adegenet compoplot
compoplot(americanadapc, lab="")

#ggcompoplots
library(ggcompoplot)
library(ggplot2)

rainbov <- setNames(rainbow(nPop(VillosaeGenlight)), popNames(VillosaeGenlight))

#3columns
ggcompoplot(americanadapc, americanaGenlight, cols = 3, pal = rainbov) + theme(axis.text.x = element_blank())

#2 columns
ggcompoplot(americanadapc, americanaGenlight, cols = 2, pal = rainbov) + theme(axis.text.x = element_blank(),axis.text.y = element_text(size=5) )

#2columns funky palette
ggcompoplot(americanadapc, americanaGenlight, cols = 2, pal = funky) + theme(axis.text.x = element_blank(),axis.text.y = element_text(size=6) )
