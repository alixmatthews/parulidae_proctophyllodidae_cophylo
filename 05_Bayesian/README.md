## Bayesian analysis using PhyloBayes-MPI


### 0. Manually convert .fasta file to .phy file
      - use some savvy find/replace in BBEdit 

### 1. Install PhyloBayes-MPI on AHPCC
```
ssh pinnacle-l6
module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh
mamba create -n phylobayes-mpi phylobayes-mpi
conda activate phylobayes-mpi
bp_mpi --version # version 1.9
```

### 2. Run multiple chains for PhyloBayes (`concatenated` directory)
`04_PhyloBayes_chain*.slurm` 

### 3. Run multiple chains for PhyloBayes (COI) in `COI` directory
`COI_only/04_PhyloBayes_chain*.slurm` 
   - I ran these on comp01 and then restarted them as I was timed out. The `restart.slurm`s are the ones I ran to restart the chain from the previous stopping spot until I got to 30k MCMCs. comp72 and comp06 were too backed up or else would have just run it on there for a full 6 hours or so!

4. Check the .trace files to see how many cycles it has gone
      - aim for 30,000 MCMC cycles

```
tail chain*.trace
```

### 5. Check  `.trace` plots in Tracer for fun

### 6. Do tracecomp to compare traces in CLI
      - `-x` flag is related to how many burn-in cycles you want, based on for fun Tracer trace plots
      - can do multiple comparisons `chain1.trace chain2.trace chain3.trace`... etc... if you want to compare more independent MCMC runs... example below

```
module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh
conda activate phylobayes-mpi 

tracecomp -x 5000 chain1.trace chain2.trace
```

#### Check the results: 
   - rel diff < 0.1 and minimum effective size>300: good run
   - rel diff < 0.3 and minimum effective size > 50: acceptable run

### 7. Do bpcomp to compare trees and get .tre files
      - `-x` flag is again related to how many burn-in cycles you want, based on for fun Tracer trace plots
      - `-o` flag is what you want the output prefix
      - again can do multiple comparisons if you want, just add them in the command... example below

```
module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh
conda activate phylobayes-mpi 

bpcomp -x 5000 chain1 chain2 -o chain12
```

#### Check the results:
   - maxdiff < 0.1: good run
   - maxdiff < 0.3: acceptable: gives a good qualitative picture of the posterior consensus
   - 0.3 < maxdiff < 1: the sample is not yet sufficiently large, and the chains have not converged, but this is on the right track
maxdiff = 1 even after 10,000 points, this indicates that at least one of the runs is stuck in a local maximum
