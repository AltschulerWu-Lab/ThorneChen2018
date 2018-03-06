"""
Set up merged-segmentation of organoid images
=============================================

Run after cluster jobs from pre-segmentation are finished
1. Check all jobs complete: check each task (subjob) folder for done file 
   (e.g. 'done.txt')
2. Combine Cellprofiler output from separate jobs for each experiment
   - Copy images into a single directory
   - Append CSV files into one file
3. Generate segmentation masks for merged-segmentation (dense regions with 
   cytoplasmic borders)
4. Setup merged-segmentation folders
5. Generates cpproj, input file list, and qsub script for each job
6. Save qsub commands in a text file

Bash command to execute qsub jobs is printed at completion
(e.g. 'bash home/organoid_analysis/scripts/20160428-005229_qsub_commands.txt'

"""

import os
import sys
os.chdir('/awlab/projects/2015_08_gut_organoid_analysis/organoid_analysis/organoid_analysis')
sys.path.append('/awlab/projects/2015_08_gut_organoid_analysis/organoid_analysis/organoid_analysis')

import cluster
import cryptops
import fileops

def main():

  # combine preseg job output
  fileops.combine_jobs('preseg')
  print('Combined all job output')

  # generate segmentation masks for merged seg
  cryptops.make_seg_masks()
  print('Generated segmentation masks')

  # setup merged seg jobs
  bash_commands_path = cluster.setup_jobs('merged')

  # terminal command for running qsub files
  print('bash {path:s}'.format(path=bash_commands_path))

  return

if __name__ == "__main__":
    main()
