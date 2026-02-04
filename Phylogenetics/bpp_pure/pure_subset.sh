for f in ./alignments/*.reduced; do
pxrms -s $f -f mixed_individuals.txt | pxs2phy > $f.pure
done
