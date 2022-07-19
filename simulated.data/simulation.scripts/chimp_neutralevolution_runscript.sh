
# script to call the prgms
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
#$ -e /dev/null
#$ -o /dev/null

#Special variables
rep=${SGE_TASK_ID}
job=${JOB_ID}

#PATHs
base="/home/hpawar"

scripts=${base}/scripts

python="/home/hpawar/conda_envs/andres-lab-msprime-env" 

slim="/share/apps/SLiM-3.2/bin/slim"

# PATH=$PATH:/share/apps/SLiM-3.2/bin/slim



outcoal=${base}/data/neutrality/coalescent_trees

outslim=${base}/data/neutrality/forward_trees

outmts=${base}/data/neutrality/mutations_out

logdir=${base}/data/neutrality/log

# for msprime outfile 
filename="coalescent_${rep}"

# make logfiles
log="${logdir}/simulationlog_${rep}.txt"
echo "begining simulation run" > ${log}

# path to python env
PATH=$PATH:/share/apps/anaconda3/bin 
## actiuvate the python environment
source activate /home/hpawar/conda_envs/andres-lab-msprime-env &>> ${log};



## paths for step3
input_for_mts=${base}/data/neutrality/forward_trees
inputfile="slimoutput_gen19305_neutrality_${rep}.trees"
outfile="3pclrinput_${rep}.txt" 


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
-d "indir='/home/hpawar/data/neutrality/coalescent_trees/'" \
-d "outdir='/home/hpawar/data/neutrality/forward_trees/'" \
-d simID=${rep} \
-d "infile='/home/hpawar/data/neutrality/coalescent_trees/coalescent_${rep}.trees'" \
${scripts}/final_slim_pandemog_neutralevolution.txt  \
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
if [ -f /home/hpawar/data/neutrality/mutations_out/3pclrinput_${rep}.txt ]
    then
        rm /home/hpawar/data/neutrality/coalescent_trees/coalescent_${rep}.trees
fi;

qstat -f -j ${job} &>> ${log};

exit 0