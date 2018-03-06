"""
Set up pre-segmentation of organoid images
==========================================


1. Setup input and output folders for each experiment
2. Generate cpproj template, raw file list, and qsub script in 'input' folder
3. Save qsub commands in a text file

Bash command to execute qsub jobs is printed at completion
(e.g. 'bash home/organoid_analysis/scripts/20160428-005229_qsub_commands.txt'
"""

import os
os.chdir('/awlab/projects/2015_08_gut_organoid_analysis/organoid_analysis/organoid_analysis')

import cluster

def main():

  # setup preseg jobs
  bash_commands_path = cluster.setup_jobs('preseg')

  # terminal command for running qsub files
  print('bash {path:s}'.format(path=bash_commands_path))

  return

if __name__ == "__main__":
    main()