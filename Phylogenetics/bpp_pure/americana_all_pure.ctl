          seed =  -1

       seqfile = alignments_reduced.phy
      Imapfile = lmap_reduced.txt
       jobname = americana_all

 speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
         speciestree = 1

*   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted



  species&tree = 8  CALYCOSA  BREVIPETALA  AMERICANA  HETERADENIA  HISPIDA  HIRSUTICAULIS  GRAYANA  RICHARDSONII
                    20  42  44  15  13  73  17  30
((CALYCOSA, (HISPIDA, (AMERICANA,HETERADENIA))), (BREVIPETALA, (HIRSUTICAULIS,(GRAYANA,RICHARDSONII))));

phase = 1 1 1 1 1 1 1 1

       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 277 * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.04 e  # Inv-gamma(a, b) for theta (integrated out by default; add E to also sample theta)
      tauprior = 3 0.2     # Inv-gamma(a, b) for root tau
    phiprior = 1 1  # Beta(a, b) for root tau & Dirichlet(a) for other tau's

      finetune =  1
      
      scaling = 1
	threads = 8
	checkpoint = 10000

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 10000
      sampfreq = 2
       nsample = 100000
