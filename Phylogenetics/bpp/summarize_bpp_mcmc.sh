# Deal with mcmc files
cp rep2_old/americana_all.mcmc.txt ./americana_all.mcmc.rep2partial.txt
cp rep4/americana_all.mcmc.txt ./americana_all.mcmc.rep4partial.txt
cat americana_all.mcmc.rep2partial.txt americana_all.mcmc.rep4partial.txt > americana_all.mcmc.combined.txt

# Perform subsequent steps on individual runs and combined
# Individual runs could conflict so combined would not be good
arr=(rep2partial rep4partial combined)

for f in "${arr[@]}"; do
# Extract species probabilities and report
echo $f
echo "Reporting things."
sed 's/.*; //g' americana_all.mcmc.${f}.txt > species_supported_vector.${f}.txt
echo "Species delimitation probabilities:" > ${f}.info.txt
awk '{ count[$1]++; total++ } END {for (k in count) printf "%s\t%.6f\n", k, count[k]/total}' species_supported_vector.${f}.txt | sort -n >> ${f}.info.txt

# Extract number of generations and report
echo "" >> ${f}.info.txt
echo "Number of generations:" >> ${f}.info.txt # Assumes the control file has sampfreq = 2
wc -l americana_all.mcmc.${f}.txt | awk '{print $1*2}' >> ${f}.info.txt

# Format to proper tree style
echo "Formatting things."
sed 's/;.*/;/g' americana_all.mcmc.${f}.txt > americana_all.mcmc.${f}.tre # Remove species number calls
sed -i 's/ //g' americana_all.mcmc.${f}.tre 
sed -i 's/#//g' americana_all.mcmc.${f}.tre # Remove branch label formatting

# Format will be θ:τ with θ not populated for all branches
# Treat τ as branch length

# Strip tip θ since newick style branch labels can't be used for tips
# Unfortunately this script assumes no numbers in species names
sed -E -i 's/([A-Za-z_]+)[0-9]+\.[0-9]+:/\1:/g' americana_all.mcmc.${f}.tre

# Select model with maximum species probability
# Get this from BPP logs
# This run we were lucky as the end is hanging on the best model
tail -n 1 americana_all.mcmc.${f}.tre > americana_all.mcmc.bestmodel.${f}.tre

# Collapse zero branches to represent the BPP algorithm (ape package in R)
# This is absolutely key because BPP handles a lumping operation as 0 coalescent branch length
# Careful with escapes.
echo "Collapsing zero branches (= BPP lumping)."
Rscript -e "
library(ape)

eps <- 1e-12   # bump to 1e-10 if needed

trs <- read.tree(\"americana_all.mcmc.${f}.tre\")

collapse0 <- function(tr, eps) {
	# If any edge lengths are NA, leave them alone; di2multi only uses numeric.
	# Collapse all (near-)zero-length internal edges into multifurcations:
	tr2 <- di2multi(tr, tol = eps)
	
	# Sometimes you need a second pass because collapsing can create new
	# near-zero internal edges in the updated structure.
	repeat {
    	prevN <- tr2\$Nnode
    	tr2 <- di2multi(tr2, tol = eps)
    	if (tr2\$Nnode == prevN) break
  		}
	tr2
	}

trs2 <- lapply(trs, collapse0, eps = eps)

write.tree(trs2, file='americana_all.mcmc.${f}.collapsed.tre', digits=16)"

# Python package dendropy tool for robust polytomy handling when calculating clade probabilities
echo "Computing clade probabilites."
sumtrees.py --target americana_all.mcmc.bestmodel.${f}.tre --decimals 0 --support-as-labels --output americana_all.mcmc.bestmodel.probabilities.${f}.tre americana_all.mcmc.${f}.collapsed.tre
done

# Generate the individual species probabilities (not directly in the tree generated)
for f in "${arr[@]}"; do
awk 'NR==1{ s=$0; while(match(s,/([A-Za-z_][A-Za-z0-9_.-]*):/,m)){ if(!(m[1] in seen)){ seen[m[1]]=1; ord[++n]=m[1] } s=substr(s,RSTART+RLENGTH) } }
{ s=$0; while(match(s,/([A-Za-z_][A-Za-z0-9_.-]*):([0-9.+-eE]+)/,m)){ t=m[1]; L=m[2]+0; tot[t]++; if(L!=0) nz[t]++; s=substr(s,RSTART+RLENGTH) } }
END{ for(i=1;i<=n;i++){ t=ord[i]; printf "%s\t%.6f\n", t, (tot[t]? nz[t]/tot[t] : 0) } }' americana_all.mcmc.rep2partial.collapsed.tre | sort
done

