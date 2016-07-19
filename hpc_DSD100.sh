#!/bin/bash

#$ -N dsd100_bss
#$ -l h_vmem=15G
#$ -cwd
#$ -q default.q
#$ -t 1-100:1
# ------------------------------------------------

# Start script
# --------------------------------
#
#printf "Starting execution of job $JOB_ID from user $SGE_O_LOGNAME\n"
#printf "Starting at `date`\n"
#printf "Calling Matlab now\n"
#printf "---------------------------\n"
echo ${SGE_TASK_ID}
/soft/MATLAB/R2013b/bin/matlab -nodisplay -nosplash -nodesktop  -singleCompThread -r "DSD100_eval_only_bss('/homedtic/mmiron/data/DSD100/','/homedtic/mmiron/data/DSD100/outputs/',$SGE_TASK_ID);"
# Copy data back, if any
#printf "---------------------------\n"
#printf "Matlab processing done.\n"
#printf "Job done. Ending at `date`\n"
