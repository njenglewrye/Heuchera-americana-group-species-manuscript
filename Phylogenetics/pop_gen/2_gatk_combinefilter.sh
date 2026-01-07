reference=H63-1.ref.fasta
prefix=americana

# Make gvcf list (must end in .list or the next step will throw an error)
ls *-g.vcf > vcf.list

# Combine and Jointly call GVCFs
gatk CombineGVCFs -R $reference --variant vcf.list --output "$prefix".cohort.g.vcf
gatk GenotypeGVCFs -R $reference -V "$prefix".cohort.g.vcf -O "$prefix".cohort.unfiltered.vcf

# Keep only SNPs passing a hard filter
time gatk SelectVariants -V "$prefix".cohort.unfiltered.vcf -R $reference -select-type-to-include SNP -O "$prefix".SNPall.vcf
#time gatk VariantFiltration -R $reference -V "$prefix".SNPall.vcf --filter-name "hardfilter" -O "$prefix".snp.filtered.vcf --filter-expression "QD < 5.0 & FS > 60.0 & MQ < 40.0 & MQRankSum < -12.5 & ReadPosRankSum < -8.0"
#gatk SelectVariants -V "$prefix".snp.filtered.vcf -O "$prefix".snp.filtered.nocall.vcf --set-filtered-gt-to-nocall

# LD estimation for next step
PopLDdecay -InVCF "$prefix".SNPall.vcf -OutStat ld_stats.txt
# I installed the plot scripts in a different way compared to documentation
PopLDdecay_Plot_OnePop.pl -inFile ld_stats.txt.stat.gz -output ld_fig

# Minor allele frequency -- 0.05 is a typical value but you can investigate with R package vcfR
# Thinning should be based on calculations of linkage disequilibrium e.g., in popLDdecay -- we eyeballed the graph for .
vcftools --vcf "$prefix".SNPall.vcf --maf 0.05 --max-missing 0.5 --thin 50 --recode --recode-INFO-all --out "$prefix".SNPall.pruned.vcf
# Kept 6733 out of 158683 SNPs

# This used PGDSpider. You will need to have that in a subfolder with the jar file noted below
java -Xmx1024m -Xms512m -jar ./PGDSpider_3.0.0.0/PGDSpider3-cli.jar -inFile "$prefix".SNPall.pruned.vcf.recode.vcf -inFormat VCF -outFile "$prefix".SNPall.pruned.str -outFormat STRUCTURE -spid template_VCF_STRUCTURE.spid
# The spid file is important. There is an example on github

# The output str file will have a blank first line and a second line for column names as genename_site -- remove these for structure
sed -i '1,2d' "$prefix".SNPall.pruned.str
