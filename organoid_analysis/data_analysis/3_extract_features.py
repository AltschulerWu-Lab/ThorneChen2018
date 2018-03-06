"""
Setup crypt masks
==================

Run after cluster jobs from merged segmentation are finished
1. Check all jobs complete: check each task (subjob) folder for done file 
   (e.g. 'done.txt')
2. Combine Cellprofiler output from separate jobs for each experiment
   - Copy images into a single directory
   - Append CSV files into one file
3. Generate crypt masks
4. Extract features

"""

import os
os.chdir('/awlab/projects/2015_08_gut_organoid_analysis/organoid_analysis/organoid_analysis')

import cryptops
import fileops
import sys

def main():

  # combine merged seg job output
  fileops.combine_jobs('merged')
  print('Combined all job output')

  # extract csv
  fileops.extract_csv('merged')
  print('Extracted CSV files')

  # make crypt objs and labels
  cryptops.make_crypt_masks('merged', 'crypt_masks')
  print('Generated crypt objects and labels')

  # measure properties
  cryptops.measure_exps('merged', 'celltype')
  print('Measured crypt properties')

if __name__ == "__main__":
  main()