library(dplyr)
library(data.table)
library(BITE)
library(OptM)
library(plyr)
source("C:/TreeMix_functions.R") #path to required functions for this analysis

########################################################
############# (A) Test migration events ################
########################################################

folder <- file.path(path="~/test_migrations")                     #path to files of TreeMix replicates with different migration edges (m) to test
test.linear = optM(folder, method = "linear", tsv="linear.txt")   #test m: produces a list of which the $out dataframe suggests optimum m based on multiple linear models
plot_optM(test.linear, method = "linear")                         #shows changes in log likelihood for different m, and suggests optimum m as 'change points'

test.optM = optM(folder, tsv ="Evanno.variance.txt")              #another option is the Evanno method - see optM package description for detailed information on output
#if data is robust and all runs have the same likelihoods, SD will be 0 and this function will give an error as it can't produce the ad hoc statistic. 
#in this case you might want to increase variance by varying -k (SNP block size), change permutation methods etc.
plot_optM(test.optM, method = "Evanno")                           #plot the proportion of variation explained by each migration event. Calculates deltaM, which is a second-order rate of change in likelihood weighted by the standard deviation

#Choose optimum number of m and continue with step 3 in the TreeMix pipeline

########################################################
################## (B) Plot Tree #######################
########################################################

## 1. From the final runs, compare tree likelihoods, select tree with highest likelihood, remove duplicates and retain tree(s) with unique topology. 
#Adapted from R functions written by Zecca, Labra and Grassi (2019).

setwd("~/final_runs")                                             #folder with all TreeMix outputs from the final runs
maxLL("finalrun_", nt=30)                                         #first argument is stem of TreeMix output files, nt = number of runs
                                                                  #shows ML trees and highest likelihood, as well as tree(s) with unique topology. Outputs "TreeLLs.txt" into workign directory

#If n of unique trees = 1, continue with step #2, if n > 1, you might want to create a consensus tree. 
#Note that bootstrap and migration values will not be available for consensus tree, thus you could also choose one ML tree 

cfTrees("finalrun_", nt=30, p=1, m='PH85')                        #m is the method to calculate pairwise distances between unique trees (default is Robinson-Foulds distance)
                                                                  #p (number from 0.5 to 1)-proportion for a clade to be represented in the consensus tree. Default (1) is strict tree, 0.5 for majority-rule
                                                                  #plots consensus tree and saves "Consensus.newick" into working directory

## 2. Now plot and save unique tree with highest likelihood:

pdf("TreeMix_output.pdf")                                          
treemix.bootstrap("finalrun_1", out.file = "tmp",                 #stem of TreeMix files (as above) + number of run with highest likelihood and unique topology
                  phylip.file = "tree_constree.newick",           #consensus tree in newick format (from the bootstrap procedure generated with PHYLIP)    
                  nboot = 500, fill = TRUE,                       #nboot is the number of bootstraps used
                  pop.color.file = "col.txt",                     #specify colours with a tab delimited pop.color.file - first column is pop name, second the colour
                  boot.legend.location = "topright")

treemix.drift(in.file = "finalrun_1",                             #pairwise matrix for drift estimates with specified order of pops 
              pop.order.color.file = "poporder.txt") + 
              title("Drift")     

plot_resid("finalrun_1",                                          #pairwise matrix for residuals
           pop_order = "poporder.txt") +
           title("Residuals")                         
dev.off()


########################################################
######### (C) Weights, Std. Err and p-values ###########
########################################################

#Output reports mean weight of edge, the jackknife estimate of the weight and standard error (averaged over N independent runs), 
#the least significant p-value recovered  over N runs for each migration event. Set working directory to the final_runs folder.

GetMigrStats(input_stem="finalrun_", nt=30)                       #arguments as above, writes file "MS_and_stats.txt" into current directory

########################################################
#### Migration support and corrected MS (Optional) #####
########################################################

#From bootstrap replicates, few other support statistic might be calculated.
#The MS is the percentage of times each pair of label sets is present among n bootstrap replicates.
#Calculated as: (number of matches / number of bootstrap replicates)*100 from all independent runs in the current working directory.
#For the Extended MS (MSE), the number of counts is corrected for multiple matches to avoid over-counting.
#Based on R funcions written by Zecca, Labra and Grassi, 2019.

GetPairsOfSets(skipL=1)                                           #create pairs of sets of populations/taxa from TreeMix output (with treeout.gz extension) in /final_runs folder, writes "PairsOfSets.txt" file
                                                                  #if you used the flag -noss, set skipL=2 (default is 1) - the number of lines to skip before reading the tree

#Now set working directory to folder with all bootstrap replicates generated with optimum number of m in Step 3.
setwd("~/final_runs/bootstrap")

#Copy PairsOfSets.txt into directory
GetMigrSupp(skipL=1)                                              #calculates MS over all bootstrap replicates, writes file "MigrSupp.txt" into current directory

GetMS_MSe(nmigr=2, min_n=2, fixed="To", skipL=1)                  #default input file is "MigrSupp.txt" created with GetMigrSupp(), writes file "MS_MSE.txt" into working directory
                                                                  #nmigr = number of migrations, fixed = specifies which taxa/species label set of each pair is kept fixed
                                                                  #fixed = "From" fixes the origin of m; fixed = "To" (default) fixes the destination of the same m 
                                                                  #min_n = minimum number of taxa/species labels to be included within the unfixed set(s)

#Ouputs table with columns 'From' (subset of species below the origin of migration edges),
#'To' (the subset of species below the destination of migration edges), Migration Support (MS) and corrected MS with respect to bootstraps (MSE).




