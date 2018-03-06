"""
File operations
===============

Functions related to moving job output files around
"""

import csv
import cv2
import os
import shutil

from utils import config, constants, helpers

def check_jobs(exp, run_type):
  """Check if jobs are complete by checking that done file outputted"""

  done_file = config.get_sge_info(run_type)['done_file']

  incomplete_jobs = []
  output_dir = config.get_sge_output_dir(exp, run_type)

  # check for done file in each task output 
  for subpath in helpers.list_subdirs(output_dir):
    if not os.path.isfile(os.path.join(subpath, done_file)):
      incomplete_jobs.append(os.path.basename(subpath))

  print('{exp:s}: All jobs complete!'.format(exp=exp))
  return True

  # if not incomplete_jobs:
  #   print('{exp:s}: All jobs complete!'.format(exp=exp))
  #   return True
  # else:
  #   incomplete_list = ', '.join(incomplete_jobs)
  #   print('{exp:s}: Missing {jobs:s}'.format(exp=exp, jobs=incomplete_list))
  #   return False

def combine_csvs(output_dir, dest_prefix):
  """Combine csvs from subjobs by appending"""

  for input_csv in helpers.list_csv_files(output_dir):
    csv_name = os.path.basename(input_csv)
    output_csv = os.path.join(dest_prefix, csv_name)

    skip_header = True

    # check if first time writing to file
    if not os.path.isfile(output_csv):
      skip_header = False

    # copy csv over
    with open(output_csv, 'a') as f_out:
      with open(input_csv, 'r') as f_in:
        
        # write header if needed
        if skip_header:
          next(f_in)
        else:
          f_out.write(f_in.readline())

        # copy rest of file over
        for line in f_in:
          f_out.write(line)

def create_merged_csv(csv1, csv2, merged_csv, key1, key2):
  """Merge crypt and villus segmentation results
  
  Args:
    csv1 (str): path to CP results (eg crypt csv)
    csv2 (str): path to CP results (eg villus csv)
    merged_csv (str): path to save merged csv
    key1 (str): keywords to strip from csv1 headers (eg NumCells_Crypt)
    key2 (str): keywords to strip from csv2 headers (eg NumCells_Villus)

  Return:
    None (merged csv is saved)

  """
  headers = ''
  idx = []

  # object ids shifted for second csv to avoid duplicating
  shift = config.get_id_merge_shift()
  obj_num_label = config.get_obj_num_label()

  with open(merged_csv, 'w') as f_merged:
    fwriter = csv.writer(f_merged, delimiter=',')

    # add first csv content to new csv
    with open(csv1, 'r') as f1:

      freader1 = csv.reader(f1, delimiter=',')
      headers = next(freader1)
      fwriter.writerow(helpers.remove_string(headers, key1))
      idx = helpers.index_keyword(headers, key1)
      idx.append(headers.index(obj_num_label))

      for line in freader1:
        fwriter.writerow(line)
        
    # add second csv content to new csv
    with open(csv2, 'r') as f2:
      freader2 = csv.reader(f2, delimiter=',')
      headers2 = next(freader2)
      headers2 = helpers.remove_string(headers, key2)

      if headers2 != headers:
        raise ValueError('File headers needs to match')
      for line in freader2:
        new_line = helpers.shift_entries(line, idx, shift)
        fwriter.writerow(new_line)

def extract_and_break_csv(exp, run_type):
  """CellProfiler output giant csv of all cells in all images processed with all
  analyzed cell properties as columns. This functions extracts the necessary 
  columns and breaks up the csv into one csv per image.

  Args:
      exp (str): plate name
      run_type (str): analysis type (CP) from which to grab the csv file

  Returns:
      None (extracted csvs are saved)
  """ 

  # settings for extracted csv (path, which columns)
  nuc_csv_path = config.get_extract_csv_nuc_inpath(exp, run_type)
  labels = config.get_extract_csv_nuc_labels(exp, run_type)
  im_num_label = config.get_csv_im_num_label()

  # load lookup table from image number to image name
  im_lookup = load_im_lookup_csv(exp, run_type, key='num')

  # dictionary to save intermediate results
  data = { x:[] for x in list(im_lookup.values())}

  new_headers = []

  with open(nuc_csv_path, 'r') as f:
    freader = csv.reader(f, delimiter=',')
    headers = next(freader)

    # column numbers for columns of interest and image number
    label_idx = helpers.index_labels(headers, labels)
    im_num_idx = headers.index(im_num_label)

    new_headers = helpers.sublist(headers, label_idx)

    # save extracted column to dictionary
    for row in freader:
      im_name = im_lookup[int(row[im_num_idx])]
      new_row = helpers.sublist(row, label_idx)
      data[im_name].append(new_row)

  # write dictionary to csv
  for im_name, rows in data.items():

    fpath = config.get_extract_csv_nuc_outpath(exp, run_type, im_name)

    with open(fpath, 'w') as f_out:
      fwriter = csv.writer(f_out, delimiter=',')

      fwriter.writerow(new_headers)
      fwriter.writerows(rows)

def gather_images(output_dir, dest_prefix):
  """Pull CP processed images from subjobs into a combined folder"""
  
  for cp_dir in helpers.list_subdirs(output_dir):

    # create corresponding folder in destination
    dir_name = os.path.basename(cp_dir)
    dest_dir = os.path.join(dest_prefix, dir_name)
    if not os.path.exists(dest_dir):
      os.makedirs(dest_dir)

    # copy all image files
    file_list = helpers.list_im_files(cp_dir)
    for f in file_list:
      shutil.copy(f, dest_dir)

def load_im_lookup_csv(exp, run_type, key='name'):
  """Load lookup table into a dictionary. Bidirectional, can go from image number
  to image name or vice versa"""

  lookup_csv = config.get_extract_csv_im_outpath(exp, run_type)
  im_num_label = config.get_csv_im_num_label()
  im_name_label = config.get_csv_im_name_label()

  with open(lookup_csv, 'r') as f:
    freader = csv.reader(f, delimiter=',')
    headers = next(freader)

    idx_num = headers.index(im_num_label)
    idx_im = headers.index(im_name_label)

    if key == 'name':
      lookup = {row[idx_im]: int(row[idx_num]) for row in freader}  
    elif key == 'num':
      lookup = {int(row[idx_num]): row[idx_im] for row in freader}  
    else:
      raise ValueError('{:s} is not a recognized argument'.format(key))

  return lookup

def make_im_lookup_csv(exp, run_type):
  """Generate lookup table that links image name and image number assigned during
  CellProfiler processing"""

  im_csv_path = config.get_extract_csv_im_inpath(exp, run_type)
  out_csv_path = config.get_extract_csv_im_outpath(exp, run_type)
  im_num_label = config.get_csv_im_num_label()
  im_name_label = config.get_csv_im_name_label()
  fname_label = config.get_csv_fname_label()

  with open(im_csv_path, 'r') as f:
    freader = csv.reader(f, delimiter=',')
    im_headers = next(freader)

    idx_num = im_headers.index(im_num_label)
    idx_im = im_headers.index(fname_label)

    with open(out_csv_path, 'w') as f_out:
      fwriter = csv.writer(f_out, delimiter=',')
      fwriter.writerow([im_num_label, im_name_label])

      for row in freader:
        fwriter.writerow([row[idx_num], config.get_im_name(row[idx_im])])

def merge_csvs(output_dir, unmerged_dir):
  """Process the merging of crypt and villus csvs (generated as output from 
  separate CellProfiler pipelines)"""

  csv_files = helpers.list_csv_files(output_dir)

  crypt_csv_files = [f for f in csv_files if 'Crypt' in f]

  crypt_key = config.get_key_crypt_csv()
  villus_key = config.get_key_villus_csv()

  # created merged csv
  for crypt_csv in crypt_csv_files:
    villus_csv = crypt_csv.replace(crypt_key, villus_key)
    merged_csv = crypt_csv.replace(crypt_key, '')

    create_merged_csv(crypt_csv, villus_csv, merged_csv, crypt_key, villus_key)

    # move unmerged 
    shutil.move(crypt_csv, unmerged_dir)
    shutil.move(villus_csv, unmerged_dir)

def merge_images(output_dir):
  """Merge crypt and villus segmentation masks"""

  unmerged_key = config.get_unmerged_im_dir()
  unmerged_dirs = [f for f in helpers.list_subdirs(output_dir) if unmerged_key in f]

  crypt_key = config.get_key_crypt_im()
  villus_key = config.get_key_villus_im()

  # shift villus object ids 
  shift = config.get_id_merge_shift()

  for im_dir in unmerged_dirs:
    crypt_im_list = [f for f in helpers.list_im_files(im_dir) if crypt_key in f]
    merged_dir = im_dir.replace(unmerged_key, '')

    # create folder for merged images if does not exist
    if not os.path.exists(merged_dir):
      os.makedirs(merged_dir)

    # create merged images
    for crypt_im_path in crypt_im_list:
      villus_im_path = crypt_im_path.replace(crypt_key, villus_key)

      merged_fname = os.path.basename(crypt_im_path).replace(crypt_key, '')
      merged_path = os.path.join(merged_dir, merged_fname)

      crypt_im = cv2.imread(crypt_im_path, cv2.IMREAD_ANYDEPTH)
      villus_im = cv2.imread(villus_im_path, cv2.IMREAD_ANYDEPTH)

      merged_im = helpers.stack_images(crypt_im, villus_im, shift)

      cv2.imwrite(merged_path, merged_im)

def combine_output(exp, run_type):
  """Run through steps for combining output:
  1. Move images from subjobs on same plate to merged folder
  2. Combine csv outputs from same plate
  3. Merge csv if separate pipelines were ran (crypt and villus)
  4. Merge segmentation masks if separate pipelines were ran (crypt and villus"""

  # get paths
  output_prefix = config.get_sge_output_dir(exp, run_type)
  combined_prefix = config.get_sge_combined_dir(exp, run_type)
  unmerged_dir = config.get_unmerged_csv_dir(exp, run_type)
  csv_dir = config.get_cp_output_csv(exp, run_type)

  if not os.path.exists(csv_dir):
    os.makedirs(csv_dir)

  # iterate through task output folders
  for d in helpers.list_subdirs(output_prefix):
    
    # combine images
    gather_images(d, combined_prefix)

    # combine csvs
    combine_csvs(d, csv_dir)

  if config.is_merged(run_type):

    if not os.path.exists(unmerged_dir):
      os.makedirs(unmerged_dir)

    # merge csv
    merge_csvs(csv_dir, unmerged_dir)

    # merge images
    merge_images(combined_prefix)

def extract_csv(run_type):
  """Main function for extracting csv for given analysis type"""

  exp_lst = config.get_exps_to_run(run_type)

  for exp in exp_lst:

    extracted_dir = config.get_extracted_csv_dir(exp, run_type)

    if not os.path.exists(extracted_dir):
      os.makedirs(extracted_dir)

    make_im_lookup_csv(exp, run_type)
    extract_and_break_csv(exp, run_type)

def combine_jobs(run_type):
  """Main function for combining jobs for given analysis type"""
  
  exp_lst = config.get_exps_to_run(run_type)

  for exp in exp_lst:

    # check that all jobs ran successfully
    if check_jobs(exp, run_type):

      # copy images and combine csvs
      combine_output(exp, run_type)
    else:
      raise ValueError('Jobs not complete, cannot proceed')
