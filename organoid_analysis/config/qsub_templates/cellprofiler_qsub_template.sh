# ============================
# QSUB TEMPLATE
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
#$ -l netapp=1G,scratch=8G         
#$ -l h_rt={{RUNTIME}}    
#$ -t {{START_JOB_NUM}}-{{END_JOB_NUM}}     

# Job information
date
hostname

export TMPDIR=/scratch
export MYTMP=`mktemp -d`
echo $MYTMP

echo "Task number: $SGE_TASK_ID"

# Setup jobs
IMG_STARTS={{IM_START_IDS}}
IMG_ENDS={{IM_END_IDS}}

FIRST_SET=${IMG_STARTS[$SGE_TASK_ID-1]}
LAST_SET=${IMG_ENDS[$SGE_TASK_ID-1]}

echo "From image number: $FIRST_SET"
echo -e "To image number: $LAST_SET\n"

echo "CellProfiler Output:"
echo "=========================================================================================="

# Setup Cellprofiler paths
CP_PATH="{{CP_PATH}}"
CP_PYTHON="{{CP_PYTHON}}"
CP_SCRIPT="{{CP_SCRIPT}}"

# Setup output paths
PREFIX="{{OUTPUT_PREFIX}}"
OUTPUT_PATH="$PREFIX/output/${FIRST_SET}_${LAST_SET}"
BATCH_PATH="$OUTPUT_PATH/Batch_data.h5"
CP_PROJ_PATH="{{CP_PROJ_PATH}}"
CP_INPUT_FILE="{{CP_INPUT_FILE}}"

JOB_DONE_FILE="{{JOB_DONE_FILE}}"

# Create task output directory
if [ ! -d "$OUTPUT_PATH" ]; then
  mkdir "$OUTPUT_PATH"
fi

# Set environment
export PATH="$CP_PATH/bin:${PATH}"
export LD_LIBRARY_PATH="$CP_PATH/lib:$CP_PATH/lib/mysql:$CP_PATH/lib64:${LD_LIBRARY_PATH}"
export PATH="$JAVA_HOME/bin:${PATH}"
export LD_LIBRARY_PATH="${JAVA_HOME}/lib/amd64:${JAVA_HOME}/lib/amd64/server:${LD_LIBRARY_PATH}"

# Run program
"$CP_PYTHON" "$CP_SCRIPT" -p "$CP_PROJ_PATH" -c -r -b --do-not-fetch -o "$OUTPUT_PATH" --file-list "$CP_INPUT_FILE"

"$CP_PYTHON" "$CP_SCRIPT" -p "$BATCH_PATH" -c -r -b --do-not-fetch -f $FIRST_SET -l $LAST_SET -t "$MYTMP"

echo "=========================================================================================="

# Job info
date

du -h $MYTMP
echo "hostname: $HOSTNAME"
echo "job id: $JOB_ID"
qstat -j $JOB_ID | grep "usage\ *$SGE_TASK_ID\:"

# Save done file
echo "Done!" > "${OUTPUT_PATH}/${JOB_DONE_FILE}"