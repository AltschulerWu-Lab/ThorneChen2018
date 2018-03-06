"""
Cluster utils
=============

Functions related to setting up to run Cellprofiler jobs on the cluster
"""

import math
import os
import shutil

from utils import config, expinfo, helpers

def add_qsub_command(command_file, qsub_path):
  """Adds a qsub command to text file"""
  command = 'qsub {:s}\n'.format(qsub_path)

  with open(command_file, 'a') as f:
    f.write(command)

def gen_cpproj_input(exp, run_type):
  """Generate list of input image paths as a text file"""

  # raw images
  raw_path = config.get_raw_data_path(exp)
  file_list = helpers.list_im_files(raw_path)
  output_filelist_path = config.get_sge_input_rawlist(exp, run_type)

  # crypt masks
  if config.is_merged(run_type):
    seg_mask_path = config.get_cp_output_seg_masks(exp)
    seg_mask_file_list = helpers.list_im_files(seg_mask_path)
    file_list = file_list + seg_mask_file_list

  with open(output_filelist_path, 'w') as f:
    for line in file_list:
      f.write('{:s}\n'.format(line))

def gen_cpproj_template(exp, run_type):
  """Copies Cellprofiler pipeline of given run type into experiment folder"""
  cpproj_template_path = config.get_template_cpproj(exp, run_type)
  cpproj_run_path = config.get_sge_input_cpproj(exp, run_type)

  shutil.copy(cpproj_template_path, cpproj_run_path)

def gen_im_job_array(start_id, end_id, num_jobs):
  """Given list of ids, generate list of starting and ending job numbers.

  Args:
      start_id (TYPE): Description
      end_id (TYPE): Description
      num_jobs (TYPE): Description

  Returns:
      {start_ids: [], end_ids: []}: Description
  """
  num_ids = end_id - start_id + 1
  num_per_job = math.ceil(num_ids / num_jobs)
  
  start_ids = list(range(start_id, end_id, num_per_job))
  end_ids = [x+num_per_job-1 for x in start_ids]

  if end_ids[-1] > end_id:
    end_ids[-1] = end_id

  return {'start_ids': start_ids, 'end_ids': end_ids}     

def gen_cp_qsub_constants(exp, run_type):
  """Generate qsub script settings for cell profiler jobs"""

  # qsub job output directory
  job_output_path = config.get_sge_job_output_dir(exp, run_type)

  # sge info
  sge_info = config.get_sge_info(run_type)
  mem_free = sge_info['mem_free']
  runtime = sge_info['runtime']

  num_jobs = sge_info['num_jobs']
  start_job_num = sge_info['start_job_num']
  end_job_num = num_jobs

  # image sets
  end_im_num = expinfo.get_num_images(exp)
  im_idx =  gen_im_job_array(1, end_im_num, num_jobs)
  im_start_ids = helpers.bash_array(im_idx['start_ids'])
  im_end_ids = helpers.bash_array(im_idx['end_ids'])
  
  # cellprofiler paths
  cp_path = config.get_cp_prefix()
  cp_python = config.get_cp_program()
  cp_script = config.get_cp_script()

  # project output
  output_prefix = config.get_processed_data_path(exp, run_type)
  cp_proj_path = config.get_sge_input_cpproj(exp, run_type)
  cp_input_file = config.get_sge_input_rawlist(exp, run_type)

  job_done_file = sge_info['done_file']

  reps = {
    'JOB_OUTPUT_PATH': job_output_path,
    'MEM_FREE': mem_free,
    'RUNTIME': runtime,    
    'START_JOB_NUM': start_job_num,
    'END_JOB_NUM': end_job_num,
    'IM_START_IDS': im_start_ids,
    'IM_END_IDS': im_end_ids,
    'CP_PATH': cp_path,
    'CP_PYTHON': cp_python,
    'CP_SCRIPT': cp_script,
    'OUTPUT_PREFIX': output_prefix,
    'CP_PROJ_PATH': cp_proj_path,
    'CP_INPUT_FILE': cp_input_file,
    'JOB_DONE_FILE': job_done_file
  }

  # convert dict keys into placeholder format ({{string}})
  reps = helpers.convert_placeholder_reps(reps)

  return reps

def gen_mask_qsub_constants(exp, run_type):
  """Generate qsub script settings for mask jobs"""
  
  # qsub job output directory
  job_output_path = config.get_sge_job_output_dir(exp, run_type)

  # sge info
  sge_info = config.get_sge_info(run_type)
  mem_free = sge_info['mem_free']
  runtime = sge_info['runtime']

  # paths
  python_path = config.get_linux_python_path()
  script_path = config.get_crypt_finder_script(run_type)
  args = exp

  job_done_file = sge_info['done_file']

  reps = {
    'JOB_OUTPUT_PATH': job_output_path,
    'MEM_FREE': mem_free,
    'RUNTIME': runtime,    
    'PYTHON_PATH': python_path,
    'SCRIPT_PATH': script_path,
    'ARGS': args,
    'JOB_DONE_FILE': job_done_file
  }

  reps = helpers.convert_placeholder_reps(reps)

  return reps

def gen_qsub_script(exp, run_type):
  """Populate qsub script with settings"""

  reps = {}
  qsub_path = ''

  if config.is_cp_job(run_type):
    reps = gen_cp_qsub_constants(exp, run_type)
  else:
    reps = gen_mask_qsub_constants(exp, run_type)

  qsub_script = ''
  qsub_template_path = config.get_template_qsub(run_type)

  # read qsub template
  with open(qsub_template_path, 'r') as f:
    qsub_template = f.read()
    qsub_script = helpers.multi_str_replace(qsub_template, reps)

  # write qsub script
  qsub_path = config.get_sge_input_qsub(exp, run_type)

  with open(qsub_path, 'w') as f:
    f.write(qsub_script)

  return qsub_path

def setup_cpproj(exp, run_type):
  """Setup cellprofiler input: cellprofiler pipeline and input images """

  gen_cpproj_template(exp, run_type)
  gen_cpproj_input(exp, run_type)

def setup_folders(exp, run_type):
  """Sets up folders for input, job logs, output, and combined output"""

  # check if directory existing
  if os.path.isdir(config.get_sge_output_dir(exp, run_type)):
    raise IOError('{exp:s} has already been run!'.format(exp=exp))

  job_folders = config.get_job_dirs()
  prefix = config.get_processed_data_path(exp, run_type)

  for d in job_folders:

    path = os.path.join(prefix, d)
    os.makedirs(path)

def setup_jobs(run_type):
  """
  Main function for setting up a cluster job batch

  This function will 
  - take list of experiments from config file
  - setup files and folders for running CellProfiler pipeline (based on run_type)

  Args: 
    run_type: 'preseg' or 'merged'

  Returns: 
    qsub_commands_path: path to text file containing commands for submitting jobs
      - to submit jobs: bash [qsub_commands_path] 

  """

  qsub_commands_path = config.get_qsub_commands_path()

  exp_lst = config.get_exps_to_run(run_type)

  # for each experiment
  for exp in exp_lst:
  
    # set up folders
    setup_folders(exp, run_type)

    # set up qsub scripts
    qsub_path = gen_qsub_script(exp, run_type)

    if config.is_cp_job(run_type):
      # put template cpproj and cpproj in input folder for cp jobs
      setup_cpproj(exp, run_type)

    # write commands to a bash script
    add_qsub_command(qsub_commands_path, qsub_path)

  return qsub_commands_path


