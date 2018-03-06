"""
Crypt Operations
================

Functions related to identifying and analyzing crypt regions
"""

import cv2
import os
import numpy as np
import pandas as pd
import scipy
from scipy import ndimage as ndi
from scipy.spatial import distance
from skimage import img_as_ubyte, morphology, measure

from utils import config, expinfo, helpers, constants


def count_masked_objs(objs, mask, label=1):
  """Count objects given particular label

  Args:
      objs (array): labeled objects (e.g. nuclei)
      mask (array): labeled regions (e.g. crypts) 
      label: mask label to filter for

  Returns: 
      int: number of objects in a mask with given label
  """
  
  masked_objs = mask_objects(objs, mask, mask_val=label)
  return len(nonzero_unique(masked_objs))

def count_stained_objs(crypt_objs, stained_objs):
  """Count number of overlap in two lists of interest

  Args:
      crypt_objs (array): first set of objects of interest
        - e.g. cells in a crypt (could be mask array or list)
      stained_objs (list): second set of objects of interest
        - e.g. cells positive for a stain

  Returns:
      int: number of overlap
  """  
  
  obj_labels = nonzero_unique(crypt_objs)
  return len(set(obj_labels) & set(stained_objs))

def measure_objs_dispersion(crypt_objs, stained_objs):
  """Count number of overlap in two lists of interest

  Args:
      crypt_objs (array): first set of objects of interest
        - e.g. cells in a crypt (could be mask array or list)
      stained_objs (list): second set of objects of interest
        - e.g. cells positive for a stain

  Returns:
      int: number of overlap
  """  

  # get edu+ objects in this crypt
  edu_obj_mask = np.isin(crypt_objs, stained_objs)
  crypt_objs[~edu_obj_mask] = 0

  # convert objects matrix to location indices

  crypt_props = measure.regionprops(crypt_objs)

  locs = [el.centroid for el in crypt_props]

  try:

    # compute distance matrix between all object pairs
    d_matrix = distance.cdist(locs, locs)
    d_matrix[d_matrix==0] = np.nan

    # nearest distance
    return np.nanmin(d_matrix, axis=1)
    
  except:
    return None

def measure_objs_test(crypt_objs, stained_objs, d_cutoff=5):
  """Count number of overlap in two lists of interest

  Args:
      crypt_objs (array): first set of objects of interest
        - e.g. cells in a crypt (could be mask array or list)
      stained_objs (list): second set of objects of interest
        - e.g. cells positive for a stain

  Returns:
      int: number of overlap
  """  

  # get edu+ objects in this crypt
  edu_obj_mask = np.isin(crypt_objs, stained_objs)
  crypt_objs[~edu_obj_mask] = 0

  # convert objects matrix to location indices

  crypt_props = measure.regionprops(crypt_objs)

  locs = [el.centroid for el in crypt_props]

  try:

    # compute distance matrix between all object pairs
    d_matrix = distance.cdist(locs, locs)
    d_matrix[d_matrix==0] = np.nan

    # nearest distance
    return np.nanmin(d_matrix, axis=1)
    
  except:
    return None

def dense_objects(objs, kernel_size=5, min_area=6000):
  """Filter for dense objects

  Args:
      objs (array): object segmentation mask (e.g. nuclei)
      kernel_size (int, optional): size of dilation kernel
      min_area (int, optional): minimum pixel area of a "dense region"

  Returns:
      ndarrary: ndarray retaining objects that are close together
  """
  
  mask_bool = helpers.logical_array(objs)

  # dilate objects
  kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernel_size, kernel_size))
  mask_dilated = cv2.dilate(mask_bool, kernel, iterations = 1)

  # drop small objects
  mask_large = morphology.remove_small_objects(mask_dilated.astype(bool), 
    min_size=min_area, connectivity=4)

  # mask over original objects, discarding objects not completely masked
  dense_objs = mask_objects(objs, mask_large)

  struct = np.ones((3,3))
  label_mask, num_feats = ndi.measurements.label(mask_large, structure=struct)

  return dense_objs, label_mask

def filter_crypts(objs, label_mask, objs_list, min_count):
  """Count objs of interest in each group"""

  crypt_labels = nonzero_unique(label_mask)
  filtered_mask = label_mask.copy()


  for l in crypt_labels:
    crypt_mask = mask_objects(objs, label_mask, mask_val=l)
    num_pos = count_stained_objs(crypt_mask, objs_list)

    if num_pos < min_count:
      remove_object(filtered_mask, l)

  filtered_crypts = objs.copy()
  filtered_crypts[filtered_mask==0] = 0

  return filtered_crypts

def filter_objects(objs, labels, none_val=0):
  """Keep objects specified by label list"""
  
  out = objs.copy()

  all_labels = set(nonzero_unique(out))
  labels = set(labels)
  remove_labels = all_labels - labels

  for l in remove_labels:
    remove_object(out, l)

  return out

def thresh_edu_objects(im_df, thresh):
  """Threshold image dataframe for EdU+ objects based on given thershold"""
  
  edu_label = config.get_csv_edu_label()
  return thresh_objects(im_df, edu_label, thresh)

def thresh_cellstain_objects(im_df, thresh):
  """Threshold image dataframe for EdU+ objects based on given thershold"""
  
  cellstain_label = config.get_csv_cellstain_label()
  return thresh_objects(im_df, cellstain_label, thresh)

def thresh_paneth_objects(im_df, thresh):
  """Threshold image dataframe for Paneth objects based on given thershold"""
  
  label = config.get_csv_paneth_label()
  return thresh_objects(im_df, label, thresh)

def thresh_objects(im_df, stain_label, thresh):
  """Filter table of cells for cells positively staining for given stain

  Args:
      im_df (dataframe): rows are cells, columns are object label, stain avg
        - output from CellProfiler
      stain_label (str): stain name (edu, paneth, etc)
      thresh (float): threshold for positive stain

  Returns:
      list: list of object labels that are positively stained
  """
  obj_label = config.get_csv_obj_num_label()
  return im_df[(im_df[stain_label]>thresh)][obj_label].tolist()

def gen_crypt_masks(exp, input_run_type, run_type):
  """Generate crypt masks based on nuclei density and EdU count

  Args:
      exp (string): name of experiment 
      input_run_type (str): cellprofiler pipeline output to use (merged, etc)
      run_type (str): name for current process ('crypt_mask')

  Returns:
      None
  """

  # get file paths
  nuc_objs_dir = config.get_cp_output_nuc_objs(exp, input_run_type)

  for nuc_path in helpers.list_im_files(nuc_objs_dir):

    # paths
    im_name = config.get_im_name(nuc_path)
    csv_path = config.get_extract_csv_nuc_outpath(exp, input_run_type, im_name)

    # read in image
    nuc_im = cv2.imread(nuc_path, cv2.IMREAD_ANYDEPTH)

    # dense nuclei
    kn_size = config.get_crypt_finder_dilate_kn_size('crypt_masks')
    min_area = config.get_crypt_finder_min_area('crypt_masks')
    nuc_dense, label_mask = dense_objects(nuc_im, kernel_size=kn_size, min_area=min_area)

    # filter by edu
    edu_thresh = expinfo.get_thresh_edu(exp)
    min_edu_num = config.get_crypt_finder_min_edu_num(run_type)
    im_data = pd.read_csv(csv_path)

    edu_pos_objs = thresh_edu_objects(im_data, edu_thresh)

    nuc_crypt = filter_crypts(nuc_dense, label_mask, edu_pos_objs, min_edu_num)
    
    # give nuclei in same crypt the same label
    nuc_labeled = label_mask.copy()
    nuc_labeled[nuc_crypt==0] = 0

    # convex hull on crypts
    # crypt_mask = chull_groups(nuc_labeled)

    # save images
    objs_fname = config.get_fpath_crypt_objs(exp, input_run_type, im_name)
    cv2.imwrite(objs_fname, nuc_crypt)

    label_fname = config.get_fpath_crypt_labels(exp, input_run_type, im_name)
    cv2.imwrite(label_fname, nuc_labeled)

def gen_seg_masks(exp, input_run_type, run_type):
  """Generate dense region masks for second segmentation (based on nuclei density)

  Args:
      exp (string): name of experiment 
      input_run_type (str): cellprofiler pipeline output to use (preseg, etc)
      run_type (str): name for current process ('preseg')

  Returns:
      None
  """

  # get file paths
  nuc_objs_dir = config.get_cp_output_nuc_objs(exp, run_type)

  # create path for directory
  out_dir = config.get_cp_output_seg_masks(exp)

  if not os.path.exists(out_dir):
    os.makedirs(out_dir)

  for nuc_path in helpers.list_im_files(nuc_objs_dir):

    # paths
    im_name = config.get_im_name(nuc_path)
    cyto_path = config.get_fpath_cyto_objs(exp, run_type, im_name)

    # read in image
    nuc_im = cv2.imread(nuc_path, cv2.IMREAD_ANYDEPTH)
    cyto_im = cv2.imread(cyto_path, cv2.IMREAD_ANYDEPTH)

    # dense nuclei
    kn_size = config.get_crypt_finder_dilate_kn_size('seg_masks')
    min_area = config.get_crypt_finder_min_area('seg_masks')
    nuc_dense, label_mask = dense_objects(nuc_im, kernel_size=kn_size, min_area=min_area)

    # cyto mask
    nuc_dense_mask = helpers.logical_array(nuc_dense)
    cyto_dense = mask_objects(cyto_im, nuc_dense_mask, complete=False)
    cyto_mask = ndi.binary_fill_holes(helpers.logical_array(cyto_dense)).astype(int)

    # save image
    filename = config.get_fpath_seg_mask(exp, run_type, im_name)
    cv2.imwrite(filename, cyto_mask)

def get_object(objs, label):
  """Return mask with only objects of given label"""
  out = objs.copy()
  out[out!=label] = 0

  return out

def mask_objects(objs, mask, mask_val=1, complete=True):
  """Remove objects in labeled array that are not entirely masked

  Args:
      objs (int ndarray): ndarray where each object is identified by integers
      mask (ndarray): mask
      mask_val (int, optional): value in mask array considered to be masked regions
      complete (bool, optional): if True, objects not completely in mask are filtered out

  Returns:
      ndarray: new array with unmasked objects removed
  """

  outside = objs.copy()
  outside[mask==mask_val] = 0
  out_ids = nonzero_unique(outside)

  not_ids = []

  if complete:
    not_ids = out_ids

  else:
    inside = objs.copy()
    inside[mask!=mask_val] = 0
    in_ids = nonzero_unique(inside)
    not_ids = set(out_ids) - set(in_ids)

  out = objs.copy()
  for i in not_ids:
    remove_object(out, i)

  return out

def measure_crypts_props(crypt_objs, label_mask, edu_objs, paneth_objs, df, row, col, fld):
  """Measure crypt level properties for all crypts in image

  Args:
      crypt_objs (array): labeled cell objects (e.g. nuclei segmentation)
      label_mask (array): labeled crypt objects 
      edu_objs (list): ids of cell objects positive for EdU
      paneth_objs (list): ids of cell objects that are paneth cells
      df (dataframe): dataframe of crypt measurements in this well 
        - add results to dataframe
      row (char): row of current well
      col (int): column of current well
      fld (int): field of current frame

  Returns:
      dataframe: dataframe with measurements from this field added
  """

  # list of crypt labels
  crypt_labels = nonzero_unique(label_mask)

  for l in crypt_labels:

    # measure properties for one crypt
    crypt_mask = get_object(label_mask, l)
    objs = mask_objects(crypt_objs, crypt_mask, mask_val=l)
    crypt_props = measure.regionprops(crypt_mask)[0]

    # add properties to dataframe
    df['num_cells'].append(len(nonzero_unique(objs)))
    df['num_edu'].append(count_stained_objs(objs, edu_objs))
    df['num_paneth'].append(count_masked_objs(paneth_objs, crypt_mask, l))
    df['nuc_area'].append(crypt_props.area)
    df['eccentricity'].append(crypt_props.eccentricity)
    df['solidity'].append(crypt_props.solidity)
    df['row'].append(row)
    df['col'].append(col)
    df['fld'].append(fld)

  return df

def measure_crypts_props_cellstain(crypt_objs, label_mask, edu_objs, cellstain_objs, df, row, col, fld):
  """Measure crypt level properties for all crypts in image

  Args:
      crypt_objs (array): labeled cell objects (e.g. nuclei segmentation)
      label_mask (array): labeled crypt objects 
      edu_objs (list): ids of cell objects positive for EdU
      cellstain_objs (list): ids of cell objects that are cellstain cells
      df (dataframe): dataframe of crypt measurements in this well 
        - add results to dataframe
      row (char): row of current well
      col (int): column of current well
      fld (int): field of current frame

  Returns:
      dataframe: dataframe with measurements from this field added
  """

  # list of crypt labels
  crypt_labels = nonzero_unique(label_mask)

  for l in crypt_labels:

    # measure properties for one crypt
    crypt_mask = get_object(label_mask, l)
    objs = mask_objects(crypt_objs, crypt_mask, mask_val=l)
    crypt_props = measure.regionprops(crypt_mask)[0]

    # add properties to dataframe
    df['num_cells'].append(len(nonzero_unique(objs)))
    df['num_edu'].append(count_stained_objs(objs, edu_objs))
    df['num_cellstain'].append(count_stained_objs(objs, cellstain_objs))
    df['nuc_area'].append(crypt_props.area)
    df['eccentricity'].append(crypt_props.eccentricity)
    df['solidity'].append(crypt_props.solidity)
    df['row'].append(row)
    df['col'].append(col)
    df['fld'].append(fld)

  return df

def measure_crypts_props_no_paneth(crypt_objs, label_mask, edu_objs, df, row, col, fld):
  """Measure crypt level properties for all crypts in image

  Args:
      crypt_objs (array): labeled cell objects (e.g. nuclei segmentation)
      label_mask (array): labeled crypt objects 
      edu_objs (list): ids of cell objects positive for EdU
      df (dataframe): dataframe of crypt measurements in this well 
        - add results to dataframe
      row (char): row of current well
      col (int): column of current well
      fld (int): field of current frame

  Returns:
      dataframe: dataframe with measurements from this field added
  """

  # list of crypt labels
  crypt_labels = nonzero_unique(label_mask)

  for l in crypt_labels:

    # measure properties for one crypt
    crypt_mask = get_object(label_mask, l)
    objs = mask_objects(crypt_objs, crypt_mask, mask_val=l)
    crypt_props = measure.regionprops(crypt_mask)[0]

    # add properties to dataframe
    df['num_cells'].append(len(nonzero_unique(objs)))
    df['num_edu'].append(count_stained_objs(objs, edu_objs))
    df['nuc_area'].append(crypt_props.area)
    df['eccentricity'].append(crypt_props.eccentricity)
    df['solidity'].append(crypt_props.solidity)
    df['row'].append(row)
    df['col'].append(col)
    df['fld'].append(fld)

  return df

def measure_crypts_props_pulsechase(crypt_objs, label_mask, edu_objs, cellstain_objs, double_objs, df, row, col, fld):
  """Measure crypt level properties for all crypts in image

  Args:
      crypt_objs (array): labeled cell objects (e.g. nuclei segmentation)
      label_mask (array): labeled crypt objects 
      edu_objs (list): ids of cell objects positive for EdU
      cellstain_objs (list): ids of cell objects that are cellstain cells
      df (dataframe): dataframe of crypt measurements in this well 
        - add results to dataframe
      row (char): row of current well
      col (int): column of current well
      fld (int): field of current frame

  Returns:
      dataframe: dataframe with measurements from this field added
  """

  # list of crypt labels
  crypt_labels = nonzero_unique(label_mask)

  for l in crypt_labels:

    # measure properties for one crypt
    crypt_mask = get_object(label_mask, l)
    objs = mask_objects(crypt_objs, crypt_mask, mask_val=l)
    crypt_props = measure.regionprops(crypt_mask)[0]

    # add properties to dataframe
    df['num_cells'].append(len(nonzero_unique(objs)))
    df['num_edu'].append(count_stained_objs(objs, edu_objs))
    df['num_cellstain'].append(count_stained_objs(objs, cellstain_objs))
    df['num_double'].append(count_stained_objs(objs, double_objs))
    df['nuc_area'].append(crypt_props.area)
    df['eccentricity'].append(crypt_props.eccentricity)
    df['solidity'].append(crypt_props.solidity)
    df['row'].append(row)
    df['col'].append(col)
    df['fld'].append(fld)

  return df

def measure_celltype(exp, run_type, drop=False):
  """Measure properties for all wells in a plate

  Args:
      exp (str): plate name
      run_type (str): analysis output to use for measurements
      drop (bool): whether to drop frames

  Returns:
      None
  """
  
  # measurement names
  crypt_measures = ['num_cells', 'num_edu', 'num_cellstain', 'nuc_area', 'eccentricity', 
    'solidity', 'row', 'col', 'fld']
  
  well_num_measures = ['num_cells', 'num_edu', 'num_cellstain', 'num_crypt_cells', 'num_crypt_edu', 
    'num_crypt_cellstain', 'num_villus_cells', 'num_villus_edu', 'num_villus_cellstain']
  well_measures = well_num_measures + ['avg_eccentricity', 'avg_solidity', 'num_crypts']

  # output directories
  crypt_dir = config.get_cp_output_crypt_measure(exp, run_type)
  well_dir = config.get_cp_output_well_measure(exp, run_type)

  # create output directories
  if not os.path.exists(crypt_dir):
    os.makedirs(crypt_dir)
  if not os.path.exists(well_dir):
    os.makedirs(well_dir)

  rows = constants.get_96_rows()
  cols = constants.get_96_cols()

  if drop:
    im_csv_path = config.get_extract_csv_im_drop_inpath(exp, run_type)

  # label
  drop_label = config.get_csv_drop_label()

  # matrices of well-level measurements      
  well_mats = {k: np.zeros((len(rows), len(cols))) for k in well_measures}

  wells = expinfo.get_wells(exp)

  for w in wells:

    row = w[0]
    col = w[1]

    c = col - 1 
    r = ord(row) - 65

    # store crypt-level measurements
    crypt_props = {k: [] for k in crypt_measures}

    # store well-level measurements
    well_nums = {k: [] for k in well_num_measures}

    # iterate over fields
    num_flds = expinfo.get_num_fields(exp)
    start = expinfo.get_field_start(exp)

    for fld in list(range(start, num_flds+start, 1)):

      im_name = config.build_im_name(exp, row, col, fld, 'dna')

      objs_path = config.get_fpath_crypt_objs(exp, run_type, im_name)
      label_path = config.get_fpath_crypt_labels(exp, run_type, im_name)
      im_data_path = config.get_extract_csv_nuc_outpath(exp, run_type, im_name)

      crypt_objs = cv2.imread(objs_path, cv2.IMREAD_ANYDEPTH)
      label_mask = cv2.imread(label_path, cv2.IMREAD_ANYDEPTH)

      im_data = pd.read_csv(im_data_path)

      edu_thresh = expinfo.get_thresh_edu(exp)
      cellstain_thresh = expinfo.get_thresh_cellstain(exp, well=[row, col])
      
      edu_objs = thresh_edu_objects(im_data, edu_thresh)
      cellstain_objs = thresh_cellstain_objects(im_data, cellstain_thresh)

      # crypt-level measurements
      crypt_props = measure_crypts_props_cellstain(crypt_objs, label_mask, edu_objs, cellstain_objs, crypt_props, row, col, fld)

      # well-level measurements
      well_props = {
        'num_cells': len(im_data),
        'num_edu': len(edu_objs),
        'num_cellstain': len(cellstain_objs),
        'num_crypt_cells': len(nonzero_unique(crypt_objs)),
        'num_crypt_edu': count_stained_objs(crypt_objs, edu_objs),
        'num_crypt_cellstain': count_stained_objs(crypt_objs, cellstain_objs),
        'row': row,
        'col': col,
        'fld': fld
        }

      well_props.update({
        'num_villus_cells': well_props['num_cells'] - well_props['num_crypt_cells'],
        'num_villus_edu': well_props['num_edu'] - well_props['num_crypt_edu'],
        'num_villus_cellstain': well_props['num_cellstain'] - well_props['num_crypt_cellstain']
        })

      for k in list(well_nums.keys()):
        well_nums[k].append(well_props[k])

      for k in list(well_nums.keys()):
        well_nums[k].append(well_props[k])

    # save well-level measurements
    for k in well_num_measures:
      well_mats[k][r][c] = np.sum(well_nums[k])

    well_mats['avg_eccentricity'][r][c] = np.mean(crypt_props['eccentricity'])
    well_mats['avg_solidity'][r][c] = np.mean(crypt_props['solidity'])
    well_mats['num_crypts'][r][c] = len(crypt_props['num_cells'])

    # write crypt-level measurements to file
    well_name = config.build_well_name(row, col)
    out_path = config.get_fpath_crypt_measure(exp, run_type, well_name)

    helpers.dict_to_csv(out_path, crypt_props)

  # write well-level measurements to file
  for k in list(well_mats.keys()):
    out_path = config.get_fpath_well_measure(exp, run_type, k)
    np.savetxt(out_path, well_mats[k], delimiter=',')

def measure_dispersion(exp, run_type, drop=False):
  """Measure properties for all wells in a plate

  Args:
      exp (str): plate name
      run_type (str): analysis output to use for measurements
      drop (bool): whether to drop frames

  Returns:
      None
  """

  # output directories
  well_dir = config.get_cp_output_well_measure(exp, run_type)

  # create output directories
  if not os.path.exists(well_dir):
    os.makedirs(well_dir)

  rows = constants.get_96_rows()
  cols = constants.get_96_cols()

  if drop:
    im_csv_path = config.get_extract_csv_im_drop_inpath(exp, run_type)

  # label
  drop_label = config.get_csv_drop_label()

  rows = ['A', 'B', 'C', 'D']
  cols = [1,2, 5]

  # measurement names
  well_disp = pd.DataFrame()

  for r, row in enumerate(rows):
    for c, col in enumerate(cols):

      # wellname
      well_name = config.build_well_name(row, col)

      # store crypt-level measurements
      crypt_disp = []

      # iterate over fields
      num_flds = expinfo.get_num_fields(exp)
      start = expinfo.get_field_start(exp)

      for fld in list(range(start, num_flds+start, 1)):

        im_name = config.build_im_name(exp, row, col, fld, 'dna')
        

        objs_path = config.get_fpath_crypt_objs(exp, run_type, im_name)
        label_path = config.get_fpath_crypt_labels(exp, run_type, im_name)
        im_data_path = config.get_extract_csv_nuc_outpath(exp, run_type, im_name)

        crypt_objs = cv2.imread(objs_path, cv2.IMREAD_ANYDEPTH)
        label_mask = cv2.imread(label_path, cv2.IMREAD_ANYDEPTH)

        im_data = pd.read_csv(im_data_path)

        edu_thresh = expinfo.get_thresh_edu(exp)
        
        edu_objs = thresh_edu_objects(im_data, edu_thresh)

        # list of crypt labels
        crypt_labels = nonzero_unique(label_mask)

        for l in crypt_labels:

          # measure properties for one crypt
          crypt_mask = get_object(label_mask, l)
          objs = mask_objects(crypt_objs, crypt_mask, mask_val=l)

          # add properties to dataframe
          edu_dist = measure_objs_dispersion(objs, edu_objs)

          # add data to dataframe
          if edu_dist is not None:
            crypt_disp = pd.DataFrame([{'d_nn': d, 'crypt_num': l, 'well': well_name} for d in edu_dist])
            well_disp = well_disp.append(crypt_disp)
            
  out_path = config.get_fpath_well_measure(exp, run_type, 'edu_dispersion')
  well_disp.to_csv(out_path, index=False)
  print(out_path)


def measure_props(exp, run_type, drop=False):
  """Measure properties for all wells in a plate

  Args:
      exp (str): plate name
      run_type (str): analysis output to use for measurements
      drop (bool): whether to drop frames

  Returns:
      None
  """
  
  # measurement names
  crypt_measures = ['num_cells', 'num_edu', 'num_paneth', 'nuc_area', 'eccentricity', 
    'solidity', 'row', 'col', 'fld']

  frame_measures = ['num_cells', 'num_edu', 'num_paneth', 'num_crypt_cells', 'num_crypt_edu', 
    'num_crypt_paneth', 'num_villus_cells', 'num_villus_edu', 'num_villus_paneth', 'num_crypts', 'row', 'col', 'fld']
  
  well_num_measures = ['num_cells', 'num_edu', 'num_paneth', 'num_crypt_cells', 'num_crypt_edu', 
    'num_crypt_paneth', 'num_villus_cells', 'num_villus_edu', 'num_villus_paneth']
  well_measures = well_num_measures + ['avg_eccentricity', 'avg_solidity', 'num_crypts']

  # output directories
  crypt_dir = config.get_cp_output_crypt_measure(exp, run_type)
  frame_dir = config.get_cp_output_frame_measure(exp, run_type)
  well_dir = config.get_cp_output_well_measure(exp, run_type)

  # create output directories
  if not os.path.exists(crypt_dir):
    os.makedirs(crypt_dir)
  if not os.path.exists(frame_dir):
    os.makedirs(frame_dir)
  if not os.path.exists(well_dir):
    os.makedirs(well_dir)

  rows = constants.get_96_rows()
  cols = constants.get_96_cols()

  if drop:
    im_csv_path = config.get_extract_csv_im_drop_inpath(exp, run_type)

  # label
  drop_label = config.get_csv_drop_label()
  count_label = config.get_csv_paneth_count_label()

  # paneth info
  paneth_csv_path = config.get_extract_csv_im_inpath(exp, 'paneth')
  paneth_index = config.get_csv_paneth_fname_label()
  paneth_data = pd.read_csv(paneth_csv_path).set_index(paneth_index)

  # matrices of well-level measurements      
  well_mats = {k: np.zeros((len(rows), len(cols))) for k in well_measures}

  for r, row in enumerate(rows):
    for c, col in enumerate(cols):

      # store crypt-level measurements
      crypt_props = {k: [] for k in crypt_measures}

      # store well-level measurements
      well_nums = {k: [] for k in well_num_measures}

      # iterate over fields
      num_flds = expinfo.get_num_fields(exp)
      start = expinfo.get_field_start(exp)

      for fld in list(range(start, num_flds+start, 1)):

        im_name = config.build_im_name(exp, row, col, fld, 'dna')
        paneth_im_name = config.build_im_name(exp, row, col, fld, 'paneth')
        paneth_im_file = config.build_im_name(exp, row, col, fld, 'paneth', suffix=True)

        objs_path = config.get_fpath_crypt_objs(exp, run_type, im_name)
        paneth_objs_path = config.get_fpath_paneth_objs(exp, 'paneth', paneth_im_name)
        label_path = config.get_fpath_crypt_labels(exp, run_type, im_name)
        im_data_path = config.get_extract_csv_nuc_outpath(exp, run_type, im_name)

        crypt_objs = cv2.imread(objs_path, cv2.IMREAD_ANYDEPTH)
        paneth_objs = cv2.imread(paneth_objs_path, cv2.IMREAD_ANYDEPTH)
        label_mask = cv2.imread(label_path, cv2.IMREAD_ANYDEPTH)

        im_data = pd.read_csv(im_data_path)

        edu_thresh = expinfo.get_thresh_edu(exp)
        # paneth_thresh = expinfo.get_thresh_paneth(exp)
        
        edu_objs = thresh_edu_objects(im_data, edu_thresh)
        # paneth_objs = thresh_paneth_objects(im_data, paneth_thresh)

        # crypt-level measurements
        crypt_props = measure_crypts_props(crypt_objs, label_mask, edu_objs, paneth_objs, crypt_props, row, col, fld)

        # well-level measurements
        well_props = {
          'num_cells': len(im_data),
          'num_edu': len(edu_objs),
          'num_paneth': paneth_data.loc[paneth_im_file][count_label],
          'num_crypt_cells': len(nonzero_unique(crypt_objs)),
          'num_crypt_edu': count_stained_objs(crypt_objs, edu_objs),
          'num_crypt_paneth': count_masked_objs(paneth_objs, scipy.sign(crypt_objs)),
          'row': row,
          'col': col,
          'fld': fld
          }

        well_props.update({
          'num_villus_cells': well_props['num_cells'] - well_props['num_crypt_cells'],
          'num_villus_edu': well_props['num_edu'] - well_props['num_crypt_edu'],
          'num_villus_paneth': well_props['num_paneth'] - well_props['num_crypt_paneth']
          })

        for k in list(well_nums.keys()):
          well_nums[k].append(well_props[k])

        for k in list(well_nums.keys()):
          well_nums[k].append(well_props[k])

      # save well-level measurements
      for k in well_num_measures:
        well_mats[k][r][c] = np.sum(well_nums[k])

      well_mats['avg_eccentricity'][r][c] = np.mean(crypt_props['eccentricity'])
      well_mats['avg_solidity'][r][c] = np.mean(crypt_props['solidity'])
      well_mats['num_crypts'][r][c] = len(crypt_props['num_cells'])

      # write crypt-level measurements to file
      well_name = config.build_well_name(row, col)
      out_path = config.get_fpath_crypt_measure(exp, run_type, well_name)

      helpers.dict_to_csv(out_path, crypt_props)

  # write well-level measurements to file
  for k in list(well_mats.keys()):
    out_path = config.get_fpath_well_measure(exp, run_type, k)
    np.savetxt(out_path, well_mats[k], delimiter=',')

def measure_props_no_paneth(exp, run_type, drop=False):
  """Measure properties for all wells in a plate

  Args:
      exp (str): plate name
      run_type (str): analysis output to use for measurements
      drop (bool): whether to drop frames

  Returns:
      None
  """
  
  # measurement names
  crypt_measures = ['num_cells', 'num_edu', 'nuc_area', 'eccentricity', 
    'solidity', 'row', 'col', 'fld']
  
  well_num_measures = ['num_cells', 'num_edu', 'num_crypt_cells', 'num_crypt_edu', 
    'num_villus_cells', 'num_villus_edu']
  well_measures = well_num_measures + ['avg_eccentricity', 'avg_solidity', 'num_crypts']

  # output directories
  crypt_dir = config.get_cp_output_crypt_measure(exp, run_type)
  well_dir = config.get_cp_output_well_measure(exp, run_type)

  # create output directories
  if not os.path.exists(crypt_dir):
    os.makedirs(crypt_dir)
  if not os.path.exists(well_dir):
    os.makedirs(well_dir)

  rows = constants.get_96_rows()
  cols = constants.get_96_cols()

  if drop:
    im_csv_path = config.get_extract_csv_im_drop_inpath(exp, run_type)

  # label
  drop_label = config.get_csv_drop_label()

  # matrices of well-level measurements      
  well_mats = {k: np.zeros((len(rows), len(cols))) for k in well_measures}

  for r, row in enumerate(rows):
    for c, col in enumerate(cols):

      # store crypt-level measurements
      crypt_props = {k: [] for k in crypt_measures}

      # store well-level measurements
      well_nums = {k: [] for k in well_num_measures}

      # iterate over fields
      num_flds = expinfo.get_num_fields(exp)
      start = expinfo.get_field_start(exp)

      for fld in list(range(start, num_flds+start, 1)):

        im_name = config.build_im_name(exp, row, col, fld, 'dna')

        objs_path = config.get_fpath_crypt_objs(exp, run_type, im_name)
        label_path = config.get_fpath_crypt_labels(exp, run_type, im_name)
        im_data_path = config.get_extract_csv_nuc_outpath(exp, run_type, im_name)

        crypt_objs = cv2.imread(objs_path, cv2.IMREAD_ANYDEPTH)
        label_mask = cv2.imread(label_path, cv2.IMREAD_ANYDEPTH)

        im_data = pd.read_csv(im_data_path)

        edu_thresh = expinfo.get_thresh_edu(exp)
        
        edu_objs = thresh_edu_objects(im_data, edu_thresh)

        # crypt-level measurements
        crypt_props = measure_crypts_props_no_paneth(crypt_objs, label_mask, edu_objs, crypt_props, row, col, fld)

        # well-level measurements
        well_props = {
          'num_cells': len(im_data),
          'num_edu': len(edu_objs),
          'num_crypt_cells': len(nonzero_unique(crypt_objs)),
          'num_crypt_edu': count_stained_objs(crypt_objs, edu_objs),
          'row': row,
          'col': col,
          'fld': fld
          }

        well_props.update({
          'num_villus_cells': well_props['num_cells'] - well_props['num_crypt_cells'],
          'num_villus_edu': well_props['num_edu'] - well_props['num_crypt_edu']
          })

        for k in list(well_nums.keys()):
          well_nums[k].append(well_props[k])

        for k in list(well_nums.keys()):
          well_nums[k].append(well_props[k])

      # save well-level measurements
      for k in well_num_measures:
        well_mats[k][r][c] = np.sum(well_nums[k])

      well_mats['avg_eccentricity'][r][c] = np.mean(crypt_props['eccentricity'])
      well_mats['avg_solidity'][r][c] = np.mean(crypt_props['solidity'])
      well_mats['num_crypts'][r][c] = len(crypt_props['num_cells'])

      # write crypt-level measurements to file
      well_name = config.build_well_name(row, col)
      out_path = config.get_fpath_crypt_measure(exp, run_type, well_name)

      helpers.dict_to_csv(out_path, crypt_props)

  # write well-level measurements to file
  for k in list(well_mats.keys()):
    out_path = config.get_fpath_well_measure(exp, run_type, k)
    np.savetxt(out_path, well_mats[k], delimiter=',')

def measure_paneth(exp, run_type, drop=False):
  """Measure properties for all wells in a plate (analysis on Paneth channel only)

  Args:
      exp (str): plate name
      run_type (str): analysis output to use for measurements
      drop (bool): whether to drop frames

  Returns:
      None
  """

  # input measurements
  im_csv_path = config.get_extract_csv_im_inpath(exp, run_type)

  if drop:
    im_csv_path = config.get_extract_csv_im_drop_inpath(exp, run_type)

  im_index = config.get_csv_paneth_fname_label()
  im_data = pd.read_csv(im_csv_path).set_index(im_index)

  # label
  count_label = config.get_csv_paneth_count_label()
  drop_label = config.get_csv_drop_label()

  # measurements
  well_num_measures = ['num_paneth']

  # output directories
  well_dir = config.get_cp_output_well_measure(exp, run_type)

  if not os.path.exists(well_dir):
    os.makedirs(well_dir)

  rows = constants.get_96_rows()
  cols = constants.get_96_cols()

  # matrices of well-level measurements      
  well_mats = {k: np.zeros((len(rows), len(cols))) for k in well_num_measures}

  for r, row in enumerate(rows):
    for c, col in enumerate(cols):
      
      # count dropped frames
      drop_count = 0

      # store well-level measurements
      well_nums = {k: [] for k in well_num_measures}

      # iterate over fields
      num_flds = expinfo.get_num_fields(exp)
      start = expinfo.get_field_start(exp)

      for fld in list(range(start, num_flds+start, 1)):

        im_name = config.build_im_name(exp, row, col, fld, 'paneth', suffix=True)

        # dropping
        if drop and constants.is_drop(im_data.loc[im_name][drop_label]):
          drop_count += 1
        else:
          well_props = {
              'num_paneth': im_data.loc[im_name][count_label]
            }

          for k in list(well_nums.keys()):
            well_nums[k].append(well_props[k])

      # save well-level measurements
      for k in well_num_measures:
        if drop:
          total = np.sum(well_nums[k])
          num = len(well_nums[k])
          
          adjusted = 0
          if num != 0:
            adjusted = total*(1+drop_count/num)
          
          well_mats[k][r][c] = adjusted
        else:
          well_mats[k][r][c] = np.sum(well_nums[k])

  # write well-level measurements to file
  for k in list(well_mats.keys()):
    out_path = config.get_fpath_well_measure(exp, run_type, k)

    if drop:
      out_path = config.get_fpath_drop_well_measure(exp, run_type, k)

    np.savetxt(out_path, well_mats[k], delimiter=',')

def measure_pulsechase(exp, run_type, drop=False):
  """Measure properties for all wells in a plate

  Args:
      exp (str): plate name
      run_type (str): analysis output to use for measurements
      drop (bool): whether to drop frames

  Returns:
      None
  """
  
  # measurement names
  crypt_measures = ['num_cells', 'num_edu', 'num_cellstain', 'num_double', 'nuc_area', 'eccentricity', 
    'solidity', 'row', 'col', 'fld']

  well_num_measures = ['num_cells', 'num_edu', 'num_cellstain', 'num_double', 'num_single_edu', 'num_single_ki67', 'num_crypt_cells', 'num_crypt_edu', 
    'num_crypt_cellstain', 'num_crypt_double', 'num_villus_cells', 'num_villus_edu', 'num_villus_cellstain', 'num_villus_double']
  well_measures = well_num_measures + ['avg_eccentricity', 'avg_solidity', 'num_crypts']

  # output directories
  crypt_dir = config.get_cp_output_crypt_measure(exp, run_type)
  well_dir = config.get_cp_output_well_measure(exp, run_type)

  # create output directories
  if not os.path.exists(crypt_dir):
    os.makedirs(crypt_dir)
  if not os.path.exists(well_dir):
    os.makedirs(well_dir)

  rows = constants.get_96_rows()
  cols = constants.get_96_cols()

  if drop:
    im_csv_path = config.get_extract_csv_im_drop_inpath(exp, run_type)

  # label
  drop_label = config.get_csv_drop_label()

  # matrices of well-level measurements      
  well_mats = {k: np.zeros((len(rows), len(cols))) for k in well_measures}

  wells = expinfo.get_wells(exp)

  for w in wells:

    row = w[0]
    col = w[1]

    c = col - 1 
    r = ord(row) - 65

    # store crypt-level measurements
    crypt_props = {k: [] for k in crypt_measures}

    # store well-level measurements
    well_nums = {k: [] for k in well_num_measures}

    # iterate over fields
    num_flds = expinfo.get_num_fields(exp)
    start = expinfo.get_field_start(exp)

    for fld in list(range(start, num_flds+start, 1)):

      im_name = config.build_im_name(exp, row, col, fld, 'dna')

      objs_path = config.get_fpath_crypt_objs(exp, run_type, im_name)
      label_path = config.get_fpath_crypt_labels(exp, run_type, im_name)
      im_data_path = config.get_extract_csv_nuc_outpath(exp, run_type, im_name)

      crypt_objs = cv2.imread(objs_path, cv2.IMREAD_ANYDEPTH)
      label_mask = cv2.imread(label_path, cv2.IMREAD_ANYDEPTH)

      im_data = pd.read_csv(im_data_path)

      edu_thresh = expinfo.get_thresh_edu(exp)
      cellstain_thresh = expinfo.get_thresh_cellstain(exp)
      
      edu_objs = thresh_edu_objects(im_data, edu_thresh)
      cellstain_objs = thresh_cellstain_objects(im_data, cellstain_thresh)
      double_objs = list(set(edu_objs).intersection(cellstain_objs))
      single_edu_objs = list(set(edu_objs).difference(cellstain_objs))
      single_ki67_objs = list(set(cellstain_objs).difference(edu_objs))

      # crypt-level measurements
      crypt_props = measure_crypts_props_pulsechase(crypt_objs, label_mask, edu_objs, cellstain_objs, double_objs, crypt_props, row, col, fld)

      # well-level measurements
      well_props = {
        'num_cells': len(im_data),
        'num_edu': len(edu_objs),
        'num_cellstain': len(cellstain_objs),
        'num_double': len(double_objs),
        'num_single_edu': len(single_edu_objs),
        'num_single_ki67': len(single_ki67_objs),
        'num_crypt_cells': len(nonzero_unique(crypt_objs)),
        'num_crypt_edu': count_stained_objs(crypt_objs, edu_objs),
        'num_crypt_cellstain': count_stained_objs(crypt_objs, cellstain_objs),
        'num_crypt_double': count_stained_objs(crypt_objs, double_objs),
        'row': row,
        'col': col,
        'fld': fld
        }

      well_props.update({
        'num_villus_cells': well_props['num_cells'] - well_props['num_crypt_cells'],
        'num_villus_edu': well_props['num_edu'] - well_props['num_crypt_edu'],
        'num_villus_cellstain': well_props['num_cellstain'] - well_props['num_crypt_cellstain'],
        'num_villus_double': well_props['num_double'] - well_props['num_crypt_double']
        })

      for k in list(well_nums.keys()):
        well_nums[k].append(well_props[k])

      for k in list(well_nums.keys()):
        well_nums[k].append(well_props[k])

    # save well-level measurements
    for k in well_num_measures:
      well_mats[k][r][c] = np.sum(well_nums[k])

    well_mats['avg_eccentricity'][r][c] = np.mean(crypt_props['eccentricity'])
    well_mats['avg_solidity'][r][c] = np.mean(crypt_props['solidity'])
    well_mats['num_crypts'][r][c] = len(crypt_props['num_cells'])

    # write crypt-level measurements to file
    well_name = config.build_well_name(row, col)
    out_path = config.get_fpath_crypt_measure(exp, run_type, well_name)

    helpers.dict_to_csv(out_path, crypt_props)

  # write well-level measurements to file
  for k in list(well_mats.keys()):
    out_path = config.get_fpath_well_measure(exp, run_type, k)
    np.savetxt(out_path, well_mats[k], delimiter=',')

def nonzero_unique(ar):
  """Return list of unique values excluding zero"""

  uniques = np.unique(ar).tolist()
  uniques.pop(0)

  return uniques

def remove_object(objs, label, none_val=0):
  """Remove object specified by id value"""
  
  objs[objs==label]=none_val
  return True

def make_crypt_masks(input_run_type, run_type):
  """Generate crypt masks for all experiments listed in run type"""
  exp_list = config.get_exps_to_run(run_type)

  for exp in exp_list:
    # create path for directory
    objs_dir = config.get_cp_output_crypt_objs(exp, input_run_type)
    labels_dir = config.get_cp_output_crypt_labels(exp, input_run_type)

    if not os.path.exists(objs_dir):
      os.makedirs(objs_dir)

    if not os.path.exists(labels_dir):
      os.makedirs(labels_dir)

    gen_crypt_masks(exp, input_run_type, run_type)

def make_seg_masks():
  """Generate segmentation masks (dense regions) for all exps in run type"""
  run_type = 'preseg'
  input_run_type = 'preseg'
  exp_list = config.get_exps_to_run(run_type)

  for exp in exp_list:
    gen_seg_masks(exp, input_run_type, run_type)

def measure_exps(run_type, measure_type, drop=False):
  """Analyze crypt and cells for all exps in run type"""
  exp_list = config.get_exps_to_run(run_type)

  for exp in exp_list:
    if measure_type == 'crypt_measure':
      measure_props(exp, run_type)
    elif measure_type == 'no_paneth':
      measure_props_no_paneth(exp, run_type)
    elif measure_type == 'paneth':
      measure_paneth(exp, run_type, drop)
    elif measure_type == 'pulsechase':
      measure_pulsechase(exp, run_type)
    elif measure_type == 'celltype':
      measure_celltype(exp, run_type)
    elif measure_type == 'dispersion':
      measure_dispersion(exp, run_type)



