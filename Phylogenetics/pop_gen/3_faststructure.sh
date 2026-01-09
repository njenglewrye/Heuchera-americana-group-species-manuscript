####fastSTRUCTURE#####
# create conda environment
# only do this once:
#conda create -n fastStructure

conda activate fastStructure

# only do this once:
#conda install bioconda::faststructure

# ask ryan how to get VCF

########################################
# Generate popfile from renaming script
########################################

rm americana_popfile.txt missing.txt
cat americana.SNPall.pruned.str | sed 's/\t.*//g' | uniq > tmp

while read f; do
if grep -q " ${f} " new_names_putative_species.txt; then
grep " ${f} " new_names_putative_species.txt >> americana_popfile.txt
elif grep -P -q " ${f}$" new_names_putative_species.txt; then
grep -P " ${f}$" new_names_putative_species.txt >> americana_popfile.txt
else
echo "$f" >> missing.txt
fi
done < tmp

cat missing.txt # Fix new_names_putative_species.txt with any deficiencies 

sed -i 's/.*alba.*/ALBA/g' americana_popfile.txt
sed -i 's/.*pubescens.*/PUBESCENS/g' americana_popfile.txt
sed -i 's/.*hispida.*/HISPIDA/g' americana_popfile.txt
sed -i 's/.*caroliniana.*/CAROLINIANA/g' americana_popfile.txt
sed -i 's/.*aceroides.*/ACEROIDES/g' americana_popfile.txt
sed -i 's/.*longiflora.*/LONGIFLORA/g' americana_popfile.txt
sed -i 's/.*HIRSUTICAULIS.*/HIRSUTICAULIS/g' americana_popfile.txt
sed -i 's/.*AMERICANA.*/AMERICANA/g' americana_popfile.txt
sed -i 's/.*CALYCOSA.*/CALYCOSA/g' americana_popfile.txt
sed -i 's/.*BREVIPETALA.*/BREVIPETALA/g' americana_popfile.txt
sed -i 's/.*GRAYANA.*/GRAYANA/g' americana_popfile.txt
sed -i 's/.*RICHARDSONII.*/RICHARDSONII/g' americana_popfile.txt
sed -i 's/.*HETERADENIA.*/HETERADENIA/g' americana_popfile.txt


#run faststructure, make sure you have a .str file of your input file (i converted .vcf to .str file using PGDSpider)
structure.py -K 2 --input=americana.SNPall.pruned --output=americanasimple --format=str

#then run multiple Ks
structure.py -K 3 --input=americana.SNPall.pruned --output=americanasimple --format=str
structure.py -K 4 --input=americana.SNPall.pruned --output=americanasimple --format=str
structure.py -K 5 --input=americana.SNPall.pruned --output=americanasimple --format=str
structure.py -K 6 --input=americana.SNPall.pruned --output=americanasimple --format=str
structure.py -K 7 --input=americana.SNPall.pruned --output=americanasimple --format=str
structure.py -K 8 --input=americana.SNPall.pruned --output=americanasimple --format=str
structure.py -K 9 --input=americana.SNPall.pruned --output=americanasimple --format=str
structure.py -K 10 --input=americana.SNPall.pruned --output=americanasimple --format=str
structure.py -K 11 --input=americana.SNPall.pruned --output=americanasimple --format=str
structure.py -K 12 --input=americana.SNPall.pruned --output=americanasimple --format=str
structure.py -K 13 --input=americana.SNPall.pruned --output=americanasimple --format=str

#finding maximum likelihood K
chooseK.py --input=americanasimple
# Model complexity that maximizes marginal likelihood = 10
# Model components used to explain structure in data = 11
# K = 11

#distruct
export MPLBACKEND=Agg # Eliminate attempt to display in X windows
distruct.py -K 11 --input=americanasimple --output=americanasimple.svg
#popfile: just a single column file with population IDs in the same order as the samples in the input file.
# pop file must EXACTLY match this order:
cat americana.SNPall.pruned.str | sed 's/\t.*//g' | uniq > sample_list_reduced.txt
# if there is an index out of bounds error check for column label (not allowed) and trailing empty lines


distruct.py -K 11 --input=americanasimple --popfile=americana_popfile.txt --output=americanasimple_popfile.svg
distruct.py -K 11 --input=americanasimple --popfile=sample_list_reduced.txt --output=americanasimple_samplenames.svg


library(pophelper)
qlist <- readQ("./results/americanasimple.11.meanQ")
qlist <- lapply(qlist, as.data.frame)
ind_labels <- readLines("./sample_list_reduced.txt")
rownames(qlist[[1]]) <- ind_labels # If there is an invalid rownames length error make sure there is not a length mismatch
md <- read.delim("americana_popfile.txt", header=FALSE, stringsAsFactors=FALSE)

qlist_df <- lapply(qlist[1:2], function(x) as.data.frame(x))
plotQ(qlist = qlist, grplab = md, sortind = "all", ordergrp = TRUE, imgtype = "pdf", basesize = 1, exportpath = getwd())
plotQ(qlist = qlist, grplab = md, sortind = "all", ordergrp = TRUE, showindlab=TRUE, useindlab=TRUE, imgtype = "pdf", basesize = 1, exportpath = getwd())
plotQ(qlist, showindlab=TRUE, useindlab=TRUE, exportpath=getwd(),sortind="all", imgtype="pdf", basesize=2)

# Output pure individual and mixed individual, with max assignment probability 90% as the cutoff.
max_q <- apply(qlist_df[[1]], 1, max)
pure_individuals  <- names(max_q)[max_q >= 0.9]
mixed_individuals <- names(max_q)[max_q <  0.9]
# Percent pure individuals based on cutoff
length(pure_individuals)/(length(pure_individuals) + length(mixed_individuals))
# Percent mixed
1-length(pure_individuals)/(length(pure_individuals) + length(mixed_individuals))
write.table(pure_individuals, "pure_individuals.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(mixed_individuals, "mixed_individuals.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


