"""
Experiment info
===============
Loads exp_info.yml
"""


import yaml

from utils import helpers

# store exp info
EXP_INFO = {}

def get_channel_name(exp, ch_type):
  return get_channels(exp)[ch_type]

def get_channel_dna(exp):
  return get_channel_name(exp, 'dna')

def get_channel_edu(exp):
  return get_channel_name(exp, 'edu')

def get_channel_paneth(exp):
  return get_channel_name(exp, 'paneth')

def get_channels(exp):
  return get_exp_prop(exp, 'channels')

def get_field_start(exp):
  return get_exp_prop(exp, 'field_start')

def get_format(exp):
  return get_exp_prop(exp, 'format')

def get_name(exp):
  return get_exp_prop(exp, 'name')

def get_num_fields(exp):
  return get_exp_prop(exp, 'num_fields')

def get_num_images(exp):
  return get_exp_prop(exp, 'num_images')

def get_thresh(exp):
  return get_exp_prop(exp, 'thresh')

def get_thresh_val(exp, stain):
  return get_thresh(exp)[stain]

def get_thresh_cellstain(exp, well=None):
  if well is None:
    return get_thresh_val(exp, 'cellstain')
  else: 
    stains = EXP_INFO[exp]['stain_name']
    for k, v in stains.items():
      if list(well) in v:
        return get_thresh_val(exp, k)
    raise ValueError('Stain name not found for well {row:s}{col:d}'.format(row=well[0], col=well[1]))

def get_thresh_edu(exp):
  return get_thresh_val(exp, 'edu')

def get_thresh_paneth(exp):
  return get_thresh_val(exp, 'paneth')

def get_wells(exp):
  return get_exp_prop(exp, 'wells')

def get_exp_list():
  return list(EXP_INFO.keys())

def get_exp_prop(exp, field):
  return EXP_INFO[exp][field]

def get_cpproj_template_name(exp, run_type):
  return EXP_INFO[exp]['cpproj_template_names'][run_type]

# load exp info file
with open('config/exp_info.yml') as f:
  EXP_INFO = yaml.safe_load(f.read())