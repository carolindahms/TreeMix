#####################################################
######### Step 3: Final runs with optimum m #########
#####################################################
# This script builds a consensus tree with the optimum number of migration events from a specified number of bootstraps and uses it to run multiple final independent runs. 
# You will need to have installed TreeMix, Parallel and PHYLIP Consense.

infile=$1           # TreeMix input file
ncore=$2 	    # maximum number of cores to use
blockk=$3 	    # SNP block size
outgroup=$4         # set outgroup, for an unrooted ML tree put here 'noRoot' (without quotes)
nboot=$5	    # number of bootstrap replicates of tree with migration
mig=$6              # number of migration events
outname=$7	    # name for output file
runs=$8             # number of independent runs (N)
tree=$9             # name of consensus tree build without migration events (in newick format)
pathP=${10}	    # path to Phylip consense program

############################
### Bootstrap procedure ####
############################

mkdir final_runs
mkdir final_runs/bootstrap

echo "**** Running TreeMix ****"
if [ $outgroup = "noRoot" ]; then
	dowork() { 
		a=$RANDOM
		b=1`date +%N`
		let "c=$a+$b"
		treemix -i $2 -bootstrap -k $3 -m $5 -se -seed $c -tf $6 -o "final_runs/bootstrap/"$4"_constree_bootrep_"$1
	}
	export -f dowork
	seq 1 $nboot | parallel -j$ncore dowork {} $infile $blockk $outname $mig $tree > $outname"_logfile_constree_bootrep.log"
else
	dowork() { 
		a=$RANDOM
		b=1`date +%N`
		let "c=$a+$b"
		treemix -i $2 -bootstrap -k $3 -m $6 -se -root $5 -seed $c -tf $7 -o "final_runs/bootstrap/"$4"_constree_bootrep_"$1
	}
	export -f dowork
	seq 1 $nboot | parallel -j$ncore dowork {} $infile $blockk $outname $outgroup $mig $tree > $outname"_logfile_constree_bootrep.log"
fi
echo "**** Running TreeMix with $6 migration events: DONE ****"


### Create a file with all the bootstrapped trees
for a in `seq 1 $nboot`
do
 bootf="final_runs/bootstrap/"$outname"_constree_bootrep_"$a".treeout.gz"
 gunzip -c $bootf | head -1 >> $outname"_boottree.tre"
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
	echo $outname"_boottree.tre" > $outname"_PhylipInputFile"
	echo "Y" >> $outname"_PhylipInputFile"
else
	# Find the position of Outgroup population
	posOutgroup=`head -1 $outname"_boottree.tre" | tr "," "\n" | grep $outgroup -n | cut -d":" -f1`
	# echo $posOutgroup
	echo $outname"_boottree.tre" > $outname"_PhylipInputFile"
	echo "O" >> $outname"_PhylipInputFile"
	echo $posOutgroup >> $outname"_PhylipInputFile"
	echo "Y" >> $outname"_PhylipInputFile"
fi

# Run Phylip
$pathP < $outname"_PhylipInputFile" > screanout

### The output from Phylip will be modified because:
### 1) TreeMix accepts trees with only one line
### 2) TreeMix accepts newick format file 

cat outtree | tr -d "\n" > $outname"_finalconstree.newick"
echo >> $outname"_finalconstree.newick"
echo "***** Phylip - consensus tree construction: DONE *****"

##################################################
#### Run TreeMix with the new consensus tree #####
##################################################
echo "***** Final $8 runs of TreeMix with $6 migration events: START *****"

if [ $outgroup = "noRoot" ]; then
	dowork() { 
		a=$RANDOM
		b=1`date +%N`
		let "c=$a+$b"
		treemix -i $2 -k $3 -global -m $5 -se -seed $c -tf $4"_finalconstree.newick" -o "final_runs/"$4_$5"m_finalrun_"$1
	}
	export -f dowork
	seq 1 $runs | parallel -j$ncore dowork {} $infile $blockk $outname $mig > $outname_$mig"m_logfile_finalruns.log"
else
	dowork() { 
		a=$RANDOM
		b=1`date +%N`
		let "c=$a+$b"
		treemix -i $2 -k $3 -global -m $6 -se -root $5 -seed $c -tf $4"_finalconstree.newick" -o "final_runs/"$4_$6"m_finalrun_"$1
	}
	export -f dowork
	seq 1 $runs | parallel -j$ncore dowork {} $infile $blockk $outname $outgroup $mig > $outname_$mig"m_logfile_finalruns.log"
fi
echo "**** Running final runs of TreeMix: DONE ****"


