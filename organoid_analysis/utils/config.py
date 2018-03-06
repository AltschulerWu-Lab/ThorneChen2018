"""
Configuration settings
======================
Loads config file and extracts properties from config file
"""

import datetime
import os
import re
import yaml

from utils import constants
from utils import expinfo
from utils import helpers

# ==============
# INITIALIZATION
# ==============

# store configuration
CONF = {}

# load config file
with open('config/config.yml') as f:
  CONF = yaml.safe_load(f.read())

# ==============
# PATH FUNCTIONS
# ==============

def format_path(template, exp='none', run_type='preseg', vals={}):
  """Format file path given template, exp, run_type """

  # check valid experiment
  if exp in expinfo.get_exp_list() or exp == 'none' :
    
    names = {'exp': exp, 'run_type': get_runtype_val(run_type), 'root': get_root_path()}
    return template.format(**names, **vals)

  else: 
    raise ValueError('{:s} is not a listed experiment'.format(exp)) 

def join_path_parts(dic, base):
  """Extract file path based on dictionary structure """
  
  prefix = dic['prefix']
  subdir = dic['subdir']

  if base == 'subdir':
    return os.path.join(prefix, subdir)    
  else:
    filename = dic[base]
    return os.path.join(prefix, subdir, filename)

# ===================
# PATHS AND FILENAMES
# ===================

# Paths

def get_processed_data_path(exp, run_type):

  template = CONF['paths']['processed_data']
  return format_path(template, exp=exp, run_type=run_type)

def get_raw_data_path(exp):
  template = CONF['paths']['raw_data']
  return format_path(template, exp=exp)

def get_root_path():
  return CONF['paths']['root']

# Filenames

def build_im_name(exp, row, col, fld, ch_name, suffix=False):
  ch = expinfo.get_channels(exp)[ch_name]
  base = CONF['name_format']['base']
  if suffix:
    ftype = CONF['name_format']['suffixes']['raw']
    return (base.format(row=row, col=col, fld=fld, ch=ch) + ftype)
  else:
    return base.format(row=row, col=col, fld=fld, ch=ch)

def build_well_name(row, col):
  return '{row:s}{col:02d}'.format(row=row, col=col)

def get_im_name(filename):
  """Gives basename image name of a string (filename)"""

  p = re.compile(CONF['name_format']['regex'])
  im_names = p.findall(filename)
  im_names = list(set(im_names))

  if len(im_names) == 1:
    first_name = p.search(filename)
    return first_name.group(0)
  elif len(im_names) > 1:
    raise ValueError('More than one image name in string')
  elif len(im_names) == 0:
    return ''

def get_well_name(filename):
  """Generates well basename based on filename"""

  p = re.compile(CONF['name_format']['regex'])
  im_name = p.match(filename)

  row = im_name.group('row')
  col = im_name.group('col')
  
  return build_well_name(row, col)

def get_row(im_name):
  return get_name_part(im_name, 'row')

def get_col(im_name):
  return get_name_part(im_name, 'col')

def get_fld(im_name):
  return get_name_part(im_name, 'fld')

def get_ch(im_name):
  return get_name_part(im_name, 'ch')

def get_name_part(im_name, part):
  """Parse part of image name"""

  p = re.compile(CONF['name_format']['regex'])
  im_name = p.search(filename)

  return im_name.group(part)

# =====
# PROGS
# =====

# Cellprofiler

def get_cp_prefix():
  return get_cp_path('prefix')

def get_cp_program():
  return get_cp_path('program')

def get_cp_script():
  return get_cp_path('script')

def get_cp_path(dest):

  os_type = helpers.get_os_type()

  prefix = CONF['progs'][os_type]['cellprofiler']['prefix']

  if dest == 'prefix':
    return format_path(prefix)
  else:
    subpath = CONF['progs'][os_type]['cellprofiler'][dest]
    return format_path(os.path.join(prefix, subpath))

# Anaconda
def get_linux_python_path():
  return format_path(CONF['progs']['linux']['python'])
  
def get_python_path():
  os_type = helpers.get_os_type()
  return format_path(CONF['progs'][os_type]['python'])

# ============
# CELLPROFILER
# ============

# CP Job Directories 

def get_job_dirs():

  job_dirs = CONF['cellprofiler']['job_dirs'].values()
  return list(job_dirs)

# CP Output Directories

def get_cp_output_csv(exp, run_type):
  return get_cp_output_dir(exp, run_type, 'csv')

def get_cp_output_crypt_labels(exp, run_type):
  return get_cp_output_dir(exp, run_type, 'crypt_labels')

def get_cp_output_crypt_measure(exp, run_type):
  return get_cp_output_dir(exp, run_type, 'crypt_measure')

def get_cp_output_crypt_objs(exp, run_type):
  return get_cp_output_dir(exp, run_type, 'crypt_objs')

def get_cp_output_cyto_objs(exp, run_type):
  return get_cp_output_dir(exp, run_type, 'cyto_objs')

def get_cp_output_frame_measure(exp, run_type):
  return get_cp_output_dir(exp, run_type, 'frame_measure')

def get_cp_output_nuc_objs(exp, run_type):
  return get_cp_output_dir(exp, run_type, 'nuc_objs')

def get_cp_output_nuc_labels(exp, run_type):
  return get_cp_output_dir(exp, run_type, 'nuc_labels')

def get_cp_output_seg_masks(exp):
  run_type = get_runtype_val('preseg')
  return get_cp_output_dir(exp, run_type, 'seg_masks')

def get_cp_output_well_measure(exp, run_type):
  return get_cp_output_dir(exp, run_type, 'well_measure')

def get_cp_output_dir(exp, run_type, file_type):
  template = join_path_parts(CONF['cellprofiler']['output_dirs'], file_type)
  return format_path(template, exp=exp, run_type=run_type)

# CP Output Files

def get_fpath_crypt_labels(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'crypt_labels')

def get_fpath_crypt_measure(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'crypt_measure')

def get_fpath_crypt_dispersion(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'crypt_dispersion')

def get_fpath_crypt_objs(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'crypt_objs')

def get_fpath_cyto_objs(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'cyto_objs')

def get_fpath_frame_measure(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'frame_measure')

def get_fpath_nuc_objs(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'nuc_objs')

def get_fpath_paneth_objs(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'paneth_objs')

def get_fpath_seg_mask(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'seg_masks')

def get_fpath_well_measure(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'well_measure')

def get_fpath_drop_well_measure(exp, run_type, im_name):
  return get_fpath_file(exp, run_type, im_name, 'drop_well_measure')

def get_fpath_file(exp, run_type, im_name, file_type):
  prefix = get_cp_output_dir(exp, run_type, file_type)
  filename = im_name + CONF['name_format']['suffixes'][file_type]
  return os.path.join(prefix, filename)


# ============
# CRYPT FINDER
# ============

def get_crypt_finder_dilate_kn_size(run_type):
  return get_crypt_finder_prop(run_type, 'dilate_kn_size')

def get_crypt_finder_edu_thresh(run_type):
  return get_crypt_finder_prop(run_type, 'edu_thresh')

def get_crypt_finder_min_area(run_type):
  return get_crypt_finder_prop(run_type, 'min_area')

def get_crypt_finder_min_edu_num(run_type):
  return get_crypt_finder_prop(run_type, 'min_edu_num')

def get_crypt_finder_paneth_thresh(run_type):
  return get_crypt_finder_prop(run_type, 'paneth_thresh')

def get_crypt_finder_prop(run_type, prop):
  return CONF['crypt_finder'][run_type][prop]

def get_crypt_finder_script(run_type):
  path_dict = get_crypt_finder_prop(run_type, 'script')
  template = join_path_parts(path_dict, 'py_script')
  return format_path(template)

# ========
# SGE INFO
# ========

def get_exps_to_run(run_type):
  return CONF['sge_info'][run_type]['exp_list']


# ========
# FILE OPS
# ========

def get_id_merge_shift():
  return CONF['fileops']['id_shift']

def get_extract_csv_nuc_labels(exp, run_type):
  return CONF['fileops']['extracted_csv_nuc']['labels']

def get_extract_csv_nuc_inpath(exp, run_type):
  prefix = get_csv_dir(exp, run_type, 'subdir')
  fname = CONF['fileops']['extracted_csv_nuc']['fname']

  return os.path.join(prefix, fname)

def get_csv_drop_label():
  return get_csv_label('drop_mark')

def get_csv_im_num_label():
  return get_csv_label('im_num')

def get_csv_im_name_label():
  return get_csv_label('im_name')

def get_csv_cellstain_label():
  return get_csv_label('mean_cellstain')

def get_csv_edu_label():
  return get_csv_label('mean_edu')

def get_csv_fname_label():
  return get_csv_label('fname')

def get_csv_obj_num_label():
  return get_csv_label('obj_num')

def get_csv_paneth_label():
  return get_csv_label('mean_paneth')

def get_csv_paneth_fname_label():
  return get_csv_label('paneth_fname')

def get_csv_paneth_count_label():
  return get_csv_label('paneth_count')

def get_csv_label(label):
  return CONF['fileops']['csv_labels'][label]

def get_extract_csv_nuc_outpath(exp, run_type, im_name):
  prefix = get_extracted_csv_dir(exp, run_type)
  im_name = get_im_name(im_name)
  fname = im_name + CONF['fileops']['extracted_csv_nuc']['out_suffix']

  return os.path.join(prefix, fname)

def get_extract_csv_im_outpath(exp, run_type):
  
  fname = CONF['fileops']['extracted_csv_im']['out_fname']
  output_prefix = get_csv_dir(exp, run_type, 'subdir')
  return os.path.join(output_prefix, fname)

def get_extract_csv_im_inpath(exp, run_type):
  
  fname = CONF['fileops']['extracted_csv_im']['in_fname']
  input_prefix = get_csv_dir(exp, run_type, 'subdir')
  return os.path.join(input_prefix, fname)

def get_extract_csv_im_drop_inpath(exp, run_type):
  
  fname = CONF['fileops']['extracted_csv_im']['drop_fname']
  input_prefix = get_csv_dir(exp, run_type, 'subdir')
  return os.path.join(input_prefix, fname)

# Filename Regex

def get_key_crypt_csv():
  return get_keyword('crypt_csv')

def get_key_crypt_im():
  return get_keyword('crypt_im')

def get_key_villus_csv():
  return get_keyword('villus_csv')

def get_key_villus_im():
  return get_keyword('villus_im')

def get_keyword(key):
  return CONF['fileops']['keywords'][key]

def get_obj_num_label():
  return CONF['fileops']['obj_num_label']

# Directories

def get_extracted_csv_dir(exp, run_type):
  return get_csv_dir(exp, run_type, 'extracted_csv')

def get_unmerged_csv_dir(exp, run_type):
  return get_csv_dir(exp, run_type, 'unmerged_csv')

def get_csv_dir(exp, run_type, csv_type):
  template = join_path_parts(CONF['fileops']['csv'], csv_type)
  return format_path(template, exp=exp, run_type=run_type)


def get_unmerged_im_dir():
  return CONF['fileops']['im']['dir_key']

# ========
# SGE INFO
# ========

def get_qsub_commands_path():
  vals = {'time_now': datetime.datetime.now()}

  return format_path(CONF['cellprofiler']['qsub_commands_path'], vals=vals)

def get_sge_info(run_type):
  return CONF['sge_info'][run_type]

# SGE Input Files

def get_sge_input_cpproj(exp, run_type):
  return get_sge_input_path(exp, run_type, 'cpproj')

def get_sge_input_qsub(exp, run_type):
  return get_sge_input_path(exp, run_type, 'qsub')

def get_sge_input_rawlist(exp, run_type):
  return get_sge_input_path(exp, run_type, 'raw_filelist')

def get_sge_input_path(exp, run_type, input_type):
  template = join_path_parts(CONF['sge_info']['input_files'], input_type)
  return format_path(template, exp=exp, run_type=run_type)

# Job type

def is_merged(run_type):
  return run_type == get_runtype_val('merged')

def is_preseg(run_type):
  return run_type == get_runtype_val('preseg')

def is_paneth(run_type):
  return run_type == get_runtype_val('paneth')

def is_cp_job(run_type):
  cp_jobs = ['merged', 'preseg', 'paneth']
  return run_type in cp_jobs

def is_mask_job(run_type):
  mask_jobs = ['crypt_masks', 'seg_masks']
  return run_type in mask_jobs

def get_runtype_val(run_type):
  return CONF['sge_info']['run_types'][run_type]

# SGE directories

def get_sge_combined_dir(exp, run_type):
  return get_sge_dir(exp, run_type, 'combined')

def get_sge_input_dir(exp, run_type):
  return get_sge_dir(exp, run_type, 'input')

def get_sge_job_output_dir(exp, run_type):
  return get_sge_dir(exp, run_type, 'job_output')

def get_sge_output_dir(exp, run_type):
  return get_sge_dir(exp, run_type, 'output')

def get_sge_dir(exp, run_type, subdir):
  prefix = get_processed_data_path(exp, run_type)
  base = CONF['cellprofiler']['job_dirs'][subdir]

  return os.path.join(prefix, base)

# SGE inputs

def get_template_cpproj(exp, run_type):
  template_name = expinfo.get_cpproj_template_name(exp, run_type)
  return get_template_path(template_name, 'cpproj_templates')

def get_template_path(field, template_type):
  return join_path_parts(CONF['cellprofiler'][template_type], field)

def get_template_qsub(run_type):
  return get_template_path(run_type, 'qsub_templates')

