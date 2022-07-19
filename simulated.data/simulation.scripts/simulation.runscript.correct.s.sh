# runscript - calls amended slim script ensures correct selection coeff is being modelled
# resource requests
# status updates
# define all necessary arguments for use in both msprime & slim

#$ -S /bin/bash
#$ -cwd
#$ -N chimp_ancestral_selection
#$ -M hpawar@bchuckle.cs.ucl.ac.uk
#$ -m besa
#$ -l h_vmem=6G,tmem=6G
#$ -l h_rt=6:00:0
#$ -l h='burns*'
#$ -e /dev/null
#$ -o /dev/null

#Special variables
rep=${SGE_TASK_ID}
job=${JOB_ID}
host=${HOSTNAME}

#ARGV
selcoeff=${1}

#PATHs
base="/SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps"

scripts="/SAN/ugi/chimp_internalbranches/scripts/scripts"

python="/home/hpawar/conda_envs/andres-lab-msprime-env" 

slim="/share/apps/SLiM-3.2/bin/slim"


outcoal=${base}/coalescent_trees

outslim=${base}/forward_trees

outmts=${base}/mutations_out

logdir=${base}/log

# for msprime outfile 
filename="coalescent_${rep}_${selcoeff}"

# make logfiles
log="${logdir}/simulationlog_${rep}_${selcoeff}.txt"
echo "beginning simulation run" > ${log}

qstat -f -j ${host} &>> ${log};

# path to python env
PATH=$PATH:/share/apps/anaconda3/bin 
## actiuvate the python environment
source activate /home/hpawar/conda_envs/andres-lab-msprime-env &>> ${log};


## paths for step3
input_for_mts=${base}/forward_trees
inputfile="slimoutput_gen19305_selection_${rep}_${selcoeff}.trees"
outfile="3pclrinput_${rep}_${selcoeff}.txt" 


## run msprime to generate the coalescent history - works!
python ${scripts}/pan_demography_truncated_before_chimp_split.py \
-outdir ${outcoal} \
-sim_file_name ${filename} \
-size 1.2 &>> ${log};
## script from josh
echo "step 1 coalescent tree made" &>> ${log};


## take that prehistory into SLim....
## call slim -works!!
${slim} \
-d "indir='/SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/coalescent_trees/'" \
-d "outdir='/SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/forward_trees/'" \
-d simID=${rep} \
-d selcoeff=${selcoeff} \
-d "infile='/SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/coalescent_trees/coalescent_${rep}_${selcoeff}.trees'" \
${scripts}/updated_pandemog_contextdependentselection.txt  \
&>> ${log};

echo "step 2 forward tree made" &>> ${log};

## take the slim output into msprime, adds mutations, downsamples and produces basic 3pclr input.

# path to script - works
python ${scripts}/SLiM.TreeSeq_to_3pclr.input.py \
    -indir ${input_for_mts} \
    -slim_ts ${inputfile} \
    -outdir ${outmts} \
    -output ${outfile} \
    &>> ${log};

echo "step 3 done, ending simulation run" &>> ${log};

# works
## remove step1 python's .trees file - keep the slim tree & the output .txt file from the 2nd python (suitable for 3pclr input)
if [ -f /SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/mutations_out/3pclrinput_${rep}_${selcoeff}.txt ]
    then
        rm /SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/coalescent_trees/coalescent_${rep}_${selcoeff}.trees
fi;

## remove initial.forward.trees file once final.forward.trees file for given rep is produced
if [ -f /SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/forward_trees/slimoutput_gen19305_selection_${rep}_${selcoeff}.trees ]
    then
        rm /SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/forward_trees/gen_before_mt_${rep}_${selcoeff}.trees
fi;

qstat -f -j ${job} &>> ${log};

exit 0