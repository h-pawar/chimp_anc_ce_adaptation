

# script to call 3pclr
# resource requests
# status updates
# define all necessary arguments for use in 3pclr

#$ -S /bin/bash
#$ -cwd
#$ -N chimp_3pclr
#$ -M hpawar@bchuckle.cs.ucl.ac.uk
#$ -m besa
#$ -l h_vmem=4G,tmem=4G
#$ -l h_rt=4:00:0
#$ -e /dev/null
#$ -o /dev/null

#Special variables
rep=${SGE_TASK_ID}
job=${JOB_ID}


#PATHs
base="/SAN/ugi/chimp_internalbranches"

scripts="/home/hpawar/scripts"

threepclr="/SAN/ugi/chimp_internalbranches/bin/3P-CLR/src/threepclr"

indir=${base}/data/neutrality/3pclr_data/genpos_morgans
outdir=${base}/data/neutrality/3pclr_data/genpos_morgans/3pclr_output

logdir=${base}/data/neutrality/3pclr_data/genpos_morgans/log

# make logfiles
log="${logdir}/3pclrlog_${rep}.txt"
echo "beginning 3pclr on neutral simulation runs" > ${log}


## call 3pclr
/SAN/ugi/chimp_internalbranches/bin/3P-CLR/src/threepclr \
"/SAN/ugi/chimp_internalbranches/data/neutrality/3pclr_data/genpos_morgans/m_aligned_polymorphic_subset_nc_3pclrinput_${rep}.txt" \
"/SAN/ugi/chimp_internalbranches/data/neutrality/3pclr_data/genpos_morgans/3pclr_output/m_aligned_polymorphic_subset_nc_3pclroutput_${rep}.txt" \
10 100 0.0025 0.05969224496824,0.126500174969109,0.314340752714501 NA \
&>> ${log};


echo "3pclr done, ending job" &>> ${log};

qstat -f -j ${job} &>> ${log};

exit 0