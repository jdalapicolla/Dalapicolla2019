* example file to simulate data: genes and traits        *
* run at the command line: ibpp-simul iBPP_example.ctl *

          seed = -30658

       seqfile = IBPP_LOCI_AN1.txt
      Imapfile = IBPP_MAP_AN1.txt
      traitfile = IBPP_MORPHO_AN1.txt
       outfile = ibpp_PRO_MOL_1_10_1_10.txt
      mcmcfile = mcmc_PRO_MOL_1_10_1_10.out

* speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 5    * speciesdelimitation algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1  * speciesdelimitation algorithm1 finetune (a m)
  speciesdelimitation = 1 1 2 1 0 1
* speciesdelimitation=1, algorithm=1, finetune (a m) = 2,1, diagnosis=0 (no), 
* startwithroot=1, i.e. the starting tree has to have speciation at the root

  uniformrootedtrees = 1         * 0 means uniform labeled histories
    species&tree =  3  HG	PRO	TEP       
                        2	74	1
                ((HG, TEP), PRO);
                
* compare to tree for simulation, where we need to specify 
* population sizes (like #0.002) and branch lengths (like :0.005):
* (((A #0.002, B) : 0.005 #.002, C #0.002) : 0.01 #.002, (D, E) :.015 #.002) : 0.02 #0.002;

*      usedata = 1    * 0: no data (prior); 1:seq & trait like
    useseqdata = 1    * 0: no seq data;     1:seq like
  usetraitdata = 0    * 0: no trait data;   1:trait like
         nloci = 250    * number of data sets in seqfile
       ntraits = 10    * number of trait variables
         nindT = 479    * total # individuals for which trait data is available
** required only if it differs from the number of individuals with genetic data

     cleandata = 0    * remove sites with ambiguity data? (1:yes, 0:no)

    thetaprior = 1 10    # gamma(a, b) for theta
      tauprior = 1 10   # gamma(a, b) for root tau & Dirichlet(a) for other tau's
           nu0 = 0         # parameters for prior on traits
        kappa0 = 0         # nu0=0 and kappa0=0 for non-informative prior

*      heredity = 1 4 4   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)
*    locusrate = 1 2.0    # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)
* sequenceerror = 0 0 0 0 0 : 0.05 1   # sequencing errors: gamma(a, b) prior

*      finetune = 1: 1 1 0.1 0.1 0.2 0.5 1.0
*      # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr, traitHsq
     finetune = 1: .01 .01 .01 .01 .01 .01 .01 .01 .01

          print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 20000
      sampfreq = 5
       nsample = 100000


*** Note: Make your window wider (144 columns) before running the program.
