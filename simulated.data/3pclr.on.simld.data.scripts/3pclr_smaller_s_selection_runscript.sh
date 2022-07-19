## apply 3P-CLR to smaller s selection simulated data
# script to call 3pclr
# resource requests
# status updates
# define all necessary arguments for use in 3pclr

#$ -S /bin/bash
#$ -cwd
#$ -N chimp_3pclr_selection
#$ -M hpawar@bchuckle.cs.ucl.ac.uk
#$ -m besa
#$ -l h_vmem=4G,tmem=4G
#$ -l h_rt=4:00:0
#$ -e /dev/null
#$ -o /dev/null

#Special variables
rep=${SGE_TASK_ID}
job=${JOB_ID}

#ARGV
selcoeff=${1}

#PATHs
base="/SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/mutations_out/subset_nc"

scripts="/SAN/ugi/chimp_internalbranches/scripts/scripts"

threepclr="/SAN/ugi/chimp_internalbranches/bin/3P-CLR/src/threepclr"

indir=${base}/genpos_morgans

outdir=${base}/genpos_morgans/3pclr_output

logdir=${base}/genpos_morgans/3pclr_log

# make logfiles
log="${logdir}/3pclrlog_${rep}_${selcoeff}.txt"
echo "beginning 3pclr on newly simulated selection simulation" > ${log}


## call 3pclr
/SAN/ugi/chimp_internalbranches/bin/3P-CLR/src/threepclr \
"/SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/mutations_out/subset_nc/genpos_morgans/m_polymorphic_subset_nc_3pclrinput_${rep}_${selcoeff}.txt" \
"/SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/mutations_out/subset_nc/genpos_morgans/3pclr_output/m_polymorphic_subset_nc_3pclroutput_${rep}_${selcoeff}.txt" \
10 100 0.0025 0.05969224496824,0.126500174969109,0.314340752714501 NA \
&>> ${log};


echo "3pclr done, ending job" &>> ${log};

qstat -f -j ${job} &>> ${log};

exit 0