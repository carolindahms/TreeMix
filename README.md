# TreeMix
Scripts to infer population aplits and mixture events from alelle frequency data using TreeMix by Pickrell & Pritchard (2012). This pipeline runs TreeMix with bootstrapping, helps choose number of migration events and creates a consensus tree. It plots the maximum likelihood tree with bootstrap values, drift and residuals and calculates statistics for every migration event, such as migration support, standard error and p-values.

*Based on scripts written by Vajana and Milanesi (2017) and R functions by Zecca, Labra and Grassi (2019).*

# Pipeline
## 1. Build consensus tree with multiple migration events

Assumes TreeMix input file, which can be created from plink file using scripts provided with the TreeMix software, or from a VCF by using Stacks with the populations module (*--treemix*).
Run `Step1_TreeMix.sh` in the command line, providing an input file, maximum number of cores, block size, outgroup (or 'noRoot' for unrooted trees), 
number of bootstrap replicates, path to PHYLIP consense program, output file name, range of migration events (m) and their number of replicates, for example:

`sh Step1_TreeMix.sh input.treemix.gz 10 100 Nipponicus 500 /appl/soft/phylip-3.697/exe/consense 3spine 1 10 10`
 
This builds a consensus tree from bootstraps and adds a specified range of m. 
Tree replicates will be stored in the *test_migrations* folder.
 
## 2. Test migration edges with OptM 

Set working directory to the *test_migrations* folder and run the R package OptM (step A) from the R script `Step2&4_TreeMix.R`.
This helps identify the optimum number of m.

## 3. Final runs with optimum number of migration edges

Run `Step3_TreeMix.sh` by providing an input file, maximum number of cores, block size, outgroup (alternatively 'noRoot' for unrooted trees), number of bootstrap replicates, number of migrations, output file name, number of independent runs (N), name of consensus tree built in Step 1, and path to consense program.

`sh Step3_TreeMix.sh input.treemix.gz 10 100 Nipponicus 500 3 3spine 30 3spine_constree.newick /appl/soft/phylip-3.697/exe/consense`

Returns trees from chosen number of independent runs with optimum number of m. 
Final runs of trees will be stored in the *final_runs* folder.

## 4. Tree visualization + Migration stats and support 

For this step you will need to have saved the file `TreeMix_functions.R`
Set working directory to the *final_runs* folder, run steps B and C from the `Step2&4_TreeMix.R` script.
From the final runs, compares tree likelihoods, plots ML tree with bootstrap values and migration weights. Returns Migration Support (MS), exact MS (MSE) and statistics such as least significant p-value from all runs, standard error and migrations weight for each migration event averaged over N runs.

### References

Milanesi, M., Capomaccio, S., Vajana, E., Bomba, L., Garcia, J.F., Ajmone-Marsan, P., Colli, L., 2017. BITE: an R package for biodiversity analyses. bioRxiv 181610. doi:10.1101/181610

Pickrell, J., & Pritchard, J. (2012). Inference of population splits and mixtures from genome-wide allele frequency data. *Nature Precedings*, 1-1. 

Zecca, G., Labra, M., & Grassi, F. (2020). Untangling the Evolution of American Wild Grapes: Admixed Species and How to Find Them. *Frontiers in Plant Science*, 10, 1814.

