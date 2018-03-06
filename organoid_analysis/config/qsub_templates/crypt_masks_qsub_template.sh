# ============================
# QSUB TEMPLATE: CRYPT MASKS
# ============================

#!/bin/bash                         

# SGE Commands                        
#$ -S /bin/bash                    
#$ -o {{JOB_OUTPUT_PATH}}                     
#$ -e {{JOB_OUTPUT_PATH}}                      
#$ -cwd                            
#$ -r y                            
#$ -j y                            
#$ -l mem_free={{MEM_FREE}}                 
#$ -l arch=linux-x64               
#$ -l netapp=1G,scratch=1G         
#$ -l h_rt={{RUNTIME}}              

# Job information
date
hostname

# Setup paths
PYTHON_PATH="{{PYTHON_PATH}}"
SCRIPT_PATH="{{SCRIPT_PATH}}"
ARGS="{{ARGS}}"

JOB_DONE_FILE="{{JOB_DONE_FILE}}"

echo "CRYPT FINDER"
echo "=========================================================================================="

# Run program
"$PYTHON_PATH" "$SCRIPT_PATH" "$ARGS"

echo "=========================================================================================="

# Job info
date
echo "hostname: $HOSTNAME"
echo "job id: $JOB_ID"
qstat -j $JOB_ID | grep "usage\ *"

# Save done file
echo "Done!" > "${OUTPUT_PATH}/${JOB_DONE_FILE}"