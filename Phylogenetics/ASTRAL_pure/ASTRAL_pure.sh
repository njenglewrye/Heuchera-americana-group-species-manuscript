# Filter gene trees by mixed individuals (<90% assignment to top gene pool in STRUCTURE) -- list is in FASTstructure script
# Then filter to at least 100 genes -- list generated from ASTRAL logs
cp /mnt/Tobtob/nick/nick_final_spruceup/*best* ./
for f in *best*; do
pxrmt -t $f -f mixed_individuals.txt > $f.pure
pxrmt -t $f.pure -f few_genes.txt > $f.pure.100genes # This one will return no matching tip labels a lot because few gene trees have the taxa.
done

cat *.pure.100genes > all_gene_trees.pure.tre
rm *best*
# Check that the gene tree count is correct
wc -l all_gene_trees.pure.tre
# Suppress branches with support less than 10 per ASTRAL documentation
nw_ed  all_gene_trees.pure.tre 'i & b<=10' o > all_gene_trees.pure.bs10.tre
# Download newest ASTRAL
git clone https://github.com/smirarab/ASTRAL
unzip ./ASTRAL/Astral*
rm -rf ASTRAL/
cp all_gene_trees.pure.bs10.tre ./Astral/
cd Astral
# Run ASTRAL
java -jar astral.5.7.8.jar -i all_gene_trees.pure.bs10.tre -o astral.nickfinal.tre 2> astral.nickfinal.tre.log
