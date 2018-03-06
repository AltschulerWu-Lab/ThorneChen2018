import glob
import os
import sys

os.chdir('/awlab/projects/2015_08_gut_organoid_analysis/organoid_analysis/organoid_analysis')
sys.path.append('/awlab/projects/2015_08_gut_organoid_analysis/organoid_analysis/organoid_analysis')

from utils import config, expinfo

def check_jobs(exp_name):
  runtype = 'merged'   # CHECK
  im_dir = 'seg_overlay'
  ext = '.png'

  sge_info = config.get_sge_info(runtype)
  num_jobs = sge_info['num_jobs']
  filenum = expinfo.get_num_images(exp_name) / num_jobs

  prefix_template = '/awlab/projects/2015_08_gut_organoid_analysis/processed_data/{exp_name:s}/{runtype:s}/output'
  prefix = prefix_template.format(exp_name=exp_name, runtype=runtype)

  subdirs = [os.path.join(prefix, d, im_dir) for d in os.listdir(prefix) if os.path.isdir(os.path.join(prefix, d))]

  print('Checking experiment: ' + exp_name + ', runtype: '+runtype)

  for s in subdirs:
    if filenum != len(glob.glob(os.path.join(s, '*'+ext))):
      print(s)

  print('Check done')

if __name__ == "__main__":
  exp_lst = ['ct_rspo_drc1', 'ct_rspo_drc2']  # CHECK
  for exp_name in exp_lst:
    check_jobs(exp_name)