seqfile = mcmc_prep_no_partition.dummy.phy  * A dummy alignment only allow to run MCMCtree
treefile = mcmc_prep_no_partition.rooted.nwk  * Rooted newick tree with annotated fossil/tip dates
mcmcfile = mcmc_prep_no_partitionrun1R.mcmctree.log  * MCMC log of parameters that can be examined in Tracer
outfile = mcmc_prep_no_partitionrun1R.mcmctree.out  * Output of the summerized results of MCMCtree

ndata = 1 * number of partitions
seqtype = 2    * 0 : nucleotides; 1: codons; 2: AAs (not required if the approximate likelihood method is used)
usedata = 2    * 0: sampling from priors with no data; 1: exact slow likelihood; 2: approximate likelihood
clock = 2     * 1: global clock with equal rates; 2: independent rates; 3: correlated rates

BDparas = 1 1 0.5 c    * birth rate, death rate, sampling priors for sampling times
finetune = 1: 0.1  0.1  0.1  0.01 .5  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErrprint = 1  * 1: normal output; 2: verbose output

*** These parameters are used for multi-loci partitioned data (ndata > 1), see dos Reis et al .(2013)

rgene_gamma = 2 2     * alpha and beta parameter of Dirichlet-gamma prior for mean rate across loci for clock=2 or 3
sigma2_gamma = 1 10   * alpha and beta parameter of Dirichlet-gamma prior for rate variance across loci for clock=2 or 3

*** These parameters control the MCMC run

burnin = 20000
sampfreq =  100
nsample =  20000

