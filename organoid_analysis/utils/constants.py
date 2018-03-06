"""
Constant settings
==================
Loads constants from constants.yml
"""

import yaml

# store constants
CONSTS = {}

def get_csv_filetypes():
  """csv extension """
  return [x.lower() for x in CONSTS['csv_filetypes']]

def get_im_filetypes():
  """image file extensions """
  return [x.lower() for x in CONSTS['im_filetypes']]

def get_os_name(os_type):
  """convert os name """
  return CONSTS['os'][os_type.lower()].lower()

def get_96_rows():
  """list of rows in 96-plate """
  return get_rows('96_well')

def get_96_cols():
  """list of columns in 96-plate """
  return get_cols('96_well')

def get_rows(plate_type):
  """list of rows given plate name """
  return CONSTS['plate_info'][plate_type]['rows']

def get_cols(plate_type):
  """list of columns given plate name """
  return CONSTS['plate_info'][plate_type]['cols']

def get_96_wells():
  """get list of wells (A01, etc) in 96-plate """
  return get_wells('96_well')

def get_wells(plate_type):
  """get list of wells given plate name """
  rows = get_rows(plate_type)
  cols = get_cols(plate_type)

  wells = []
  for r in rows:
    for c in cols:
      wells.append('{row:s}{col:02d}'.format(row=r, col=c))

  return wells

def is_drop(s):
  return s == CONSTS['drop_mark']

# load constants file
with open('config/constants.yml') as f:
  CONSTS = yaml.safe_load(f.read())