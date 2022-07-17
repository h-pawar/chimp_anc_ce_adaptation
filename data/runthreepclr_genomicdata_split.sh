#!/bin/bash
#$ -S /bin/bash #defines the used shell
#$ -cwd
#$ -V
#$ -m e
#$ -l h_vmem=1.0G,s_vmem=1.0G,mem_total=1.0G,h_rt=7:00:00
# Job Array with jobs corresponding to chromosome chunks

#$ -e /dev/null 
#$ -o /dev/null
## ARGV
chr=${1}

## directories
src="/SAN/ugi/chimp_internalbranches/bin/3P-CLR/src"
files="/SAN/ugi/chimp_internalbranches/data/genomicdata_unzipped/genomes_subset_perchr_nigeria/3pclr_input"
out="/SAN/ugi/chimp_internalbranches/data/genomicdata_unzipped/genomes_subset_perchr_nigeria/3pclr_output/split_3pclroutput"
logdir="/SAN/ugi/chimp_internalbranches/data/genomicdata_unzipped/genomes_subset_perchr_nigeria/log/split_log"

## make log files
log=${logdir}/chr${chr}_chunk${SGE_TASK_ID}_log.txt;

## focal range based on task id
if [ ${SGE_TASK_ID} -eq 1 ]
then firstfocal=0
     lastfocal=$(($firstfocal+5000));
fi
if [ ${SGE_TASK_ID} -gt 1 ]
then idx=$(($SGE_TASK_ID-1));
    firstfocal=$(($idx*5000));
    offset=$(($idx*10));
    firstfocal=$(($firstfocal+offset));
    lastfocal=$(($firstfocal+5000));
fi
sleep ${SGE_TASK_ID}; # pause for time in seconds equal to task id. as the c program gets a seed from the time, this helps each job have an indpendent seed.
echo "start ${chr} $firstfocal $lastfocal" ${log};
${src}/threepclr_idxrange "${files}/aligned_chr${chr}.derived.allele.counts.central.eastern.nigeria.epoAA_withGenPosM_NCpolymorphic.txt" "${out}/chr${chr}.3pclr.out.central.eastern.nigeria.epoAA_chunck${SGE_TASK_ID}.txt" 10 100 0.0025 0.0905256105285065,0.137699478522873,0.535909737822827 NA ${firstfocal} ${lastfocal} &>> ${log}

exit 0

