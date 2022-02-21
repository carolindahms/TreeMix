#####################################################
###### Step 1: OptM prep with consensus tree ########
#####################################################
# This script builds a consensus tree from a specified number of bootstraps and uses it to run multiple replicates of trees with a specified range of migrations. 
# Treemix replicates will be used to test for optimum number of migrations in the next step.
# You will need to have installed TreeMix, Parallel and PHYLIP Consense.

infile=$1 		# TreeMix input file
ncore=$2 		# maximum number of cores to use
blockk=$3 		# SNP block size
outgroup=$4 	        # name of the selected outgroup population (for an unrooted ML tree put here 'noRoot' (without quotes))
nboot=$5		# number of bootstrap replicates of tree
pathP=$6		# path to Phylip consense program. Example: /biosoftware/phylip/phylip-3.696/exe/consense
outname=$7		# name for output file
minM=$8                 # minimum number of migrations to be tested with OptM (minM >= 1 as Treemix automatically tests for 0 migration events) 
maxM=$9                 # maximum number of migrations to be tested
migrep=${10}            # number of replicates to test each number of migrations with OptM (no less than 2)

############################
###### Settings file #######
############################

echo "Input file name = $1" > $outname"_Step1_Settings.txt"
echo "Output file name = $7" >> $outname"_Step1_Settings.txt"
echo "Number of SNPS per block = $3" >> $outname"_Step1_Settings.txt"
if [ $outgroup = "noRoot" ]; then
	echo "Unrooted ML tree" >> $outname"_Step1_Settings.txt"
else
	echo "Outgroup = $4" >> $outname"_Step1_Settings.txt"
fi
echo "Consensus tree built from $5 bootstrap replicates" >> $outname"_Step1_Settings.txt"
echo "$8 to $9 migration events tested with ${10} replicates each" >> $outname"_Step1_Settings.txt"

mkdir bootstrap

############################
### Bootstrap procedure ####
############################

echo "**** Running TreeMix ****"
if [ $outgroup = "noRoot" ]; then
	dowork() { 
		a=$RANDOM
		b=1`date +%N`
		let "c=$a+$b"
		treemix -i $2 -bootstrap -k $3 -se -seed $c -o "bootstrap/"$4"_constree_bootrep_"$1
	}
	export -f dowork
	seq 1 $nboot | parallel -j$ncore dowork {} $infile $blockk $outname > $outname"_logfile_constree_bootrep.log"
else
	dowork() { 
		a=$RANDOM
		b=1`date +%N`
		let "c=$a+$b"
		treemix -i $2 -bootstrap -k $3 -se -root $5 -seed $c -o "bootstrap/"$4"_constree_bootrep_"$1
	}
	export -f dowork
	seq 1 $nboot | parallel -j$ncore dowork {} $infile $blockk $outname $outgroup > $outname"_logfile_constree_bootrep.log"
fi
echo "**** Running TreeMix: DONE ****"


### Create a file with all the bootstrapped trees
for a in `seq 1 $nboot`
do
 bootf="bootstrap/"$outname"_constree_bootrep_"$a".treeout.gz"
 gunzip -c $bootf | head -1 >> $outname"_bootconstree.tre"
done

echo "***** Bootstrap procedure: DONE *****"

#########################################################################
#### Run PHYLIP on the bootstrapped trees to obtain a consensus tree ####
#########################################################################
echo "***** Phylip - consensus tree construction: START *****"
### Clean the environment
rm -rf outfile outtree screanout

# Create parameters file
if [ $outgroup = "noRoot" ]; then
	echo $outname"_bootconstree.tre" > $outname"_PhylipInputFile"
	echo "Y" >> $outname"_PhylipInputFile"
else
	# Find the position of Outgroup population
	posOutgroup=`head -1 $outname"_bootconstree.tre" | tr "," "\n" | grep $outgroup -n | cut -d":" -f1`
	# echo $posOutgroup
	echo $outname"_bootconstree.tre" > $outname"_PhylipInputFile"
	echo "O" >> $outname"_PhylipInputFile"
	echo $posOutgroup >> $outname"_PhylipInputFile"
	echo "Y" >> $outname"_PhylipInputFile"
fi

# Run Phylip
$pathP < $outname"_PhylipInputFile" > screanout

### The output from Phylip will be modified because:
### 1) TreeMix accepts trees with only one line
### 2) TreeMix accepts newick format file 

cat outtree | tr -d "\n" > $outname"_constree.newick"
echo >> $outname"_constree.newick"
echo "***** Phylip - consensus tree construction: DONE *****"

##################################################
###### Run TreeMix with the consensus tree #######
##################################################

echo "**** Running TreeMix with consensus tree adding $8 to $9 migration edges****"
mkdir test_migrations
if [ $outgroup = "noRoot" ]; then
                                 for ((m=minM; m<=maxM; m++))
                                    do
                                    for ((i=1; i<=migrep; i++))
                                       do
                                       treemix \
                                            -i $infile \
                                            -o ./test_migrations/treemix.${m}.${i} \
                                            -global \
                                            -m ${m} \
                                            -k $blockk \
				            -tf $outname"_constree.newick" 
                                       done 
                                 done

else
for ((m=minM; m<=maxM; m++))
   do
   for ((i=1; i<=migrep; i++))
      do
      treemix \
           -i $infile \
           -o ./test_migrations/treemix.${m}.${i} \
           -global \
           -m ${m} \
           -k $blockk \
	   -root $outgroup \
           -tf $outname"_constree.newick" 
      done 
done
fi

echo "**** TreeMix - Bootstrap Analysis with migrations: DONE ****"
