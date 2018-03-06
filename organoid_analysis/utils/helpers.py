import collections
import csv
import math
import numpy as np
import os
import sys

from utils import constants

def check_positives(lst, idx):
  """Binary array of positive values given value and index list """
  return not [x for i,x in enumerate(lst) if (i in idx) and int(x)<1]

def convert_placeholder_reps(reps):
  """Converts dictionary keys to placeholder format and dictionary values to strings"""
  
  for old_key in list(reps.keys()):
    new_key = make_placeholder(old_key)
    reps[new_key] = str(reps.pop(old_key))

  return reps

def dict_to_csv(fpath, dic):
  """Write dictionary of values to csv file """
  with open(fpath, 'w') as f:
    fwriter = csv.writer(f)
    keys = list(dic.keys())
    vals = list(zip(*list(dic.values())))
    fwriter.writerow(keys)
    fwriter.writerows(vals)

def get_os_type():
  """Gives os type ('mac', 'linux') as output"""
  
  os_name = sys.platform.lower()

  if os_name == constants.get_os_name('linux'):
    return 'linux'
  
  elif os_name == constants.get_os_name('mac'):
    return 'mac'

def bash_array(lst):
  """Converts python array [a, b, c] to bash array (a b c)"""

  contents = ' '.join(str(x) for x in lst)
  return '({:s})'.format(contents)

def index_keyword(lst, keyword):
  """Index of all list members that include keyword """
  return [i for i, x in enumerate(lst) if keyword in x]

def index_labels(lst, labels):
  """Index of all list members contained in second given list"""
  return [i for i, x in enumerate(lst) if x in labels]

def is_csv_file(fpath, filename):
  """Check if file is csv file """
  csv_filetypes = constants.get_csv_filetypes()
  return is_file_type(fpath, filename, csv_filetypes)

def is_im_file(fpath, filename):
  """Check if file is image type file """
  im_filetypes = constants.get_im_filetypes()
  return is_file_type(fpath, filename, im_filetypes)

def is_file_type(fpath, filename, ext_list):
  """Returns true if file is valid, not hidden, and has extension of given type"""

  file_parts = filename.split('.')

  # invalid file
  if not os.path.isfile(os.path.join(fpath, filename)):
    return False

  # hidden file
  elif filename.startswith('.'):
    return False

  # no extension
  elif len(file_parts) < 2:
    return False
  
  # check file type
  extension = file_parts[-1].lower()

  if extension in ext_list:
    return True  
  else:
    return False
  
def list_csv_files(fpath):
  """List all csv files in given directory path """
  file_list = [os.path.join(fpath, f) for f in os.listdir(fpath) if is_csv_file(fpath, f)]
  return file_list

def list_im_files(fpath):
  """List all image files in given directory path """
  file_list = [os.path.join(fpath, f) for f in os.listdir(fpath) if is_im_file(fpath, f)]
  return file_list

def list_subdirs(fpath):
  """Returns list of full paths to subdirectories in a given path"""

  return [os.path.join(fpath, f) for f in os.listdir(fpath) if os.path.isdir(os.path.join(fpath, f))]

def logical_array(ar):
  """Convert ndarray (int, float, bool) to array of 1 and 0's"""
  
  out = ar.copy()
  out[out!=0] = 1

  return out

def negate_entries(lst, idx):
  """Negate subset of list members indicated by indices list """
  return [-int(x) if i in idx else x for i, x in enumerate(lst)]

def make_placeholder(s):
  """Converts a string to placeholder format {{string}}"""

  return '{{'+s+'}}'

def multi_str_replace(text, reps):
  """Makes multiple string substitutions

  Args:
      text (str): Original string
      reps (dict): Dictionary with keys as substrings to be replaced and values are replacement substrings

  Returns:
      str: new string
  """

  for old, new in reps.items():
    text = text.replace(old, new)

  return text

def remove_string(str_list, substr):
  """Given list of strings, remove substring from each string"""
  return [s.replace(substr, '') for s in str_list]

def shift_entries(lst, idx, shift):
  """Shift select entries by [shift] """
  return [int(x)+shift if i in idx else x for i, x in enumerate(lst)]

def stack_images(im1, im2, shift):
  """Merge two images"""
  out = im2.copy()
  out[im1 != 0] = 0
  out[out != 0] += shift
  return np.add(out, im1)

def sublist(lst, idx):
  """Return sublist defined by indices"""
  return [lst[i] for i in idx]
