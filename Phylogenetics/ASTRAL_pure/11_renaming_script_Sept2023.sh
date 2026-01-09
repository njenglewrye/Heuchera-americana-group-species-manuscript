for f in *.tre; do
g=`echo $f | sed 's/\.tre//g'`
pxrlt -t $f -c old_names.txt -n new_names_putative_species.txt > $g.renamed.tre
done