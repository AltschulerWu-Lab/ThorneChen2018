# 2D Enteroid Analysis 

Data analysis code for [Thorne, Chen, et al (2018 Dev Cell)](https://www.ncbi.nlm.nih.gov/pubmed/29503158)

## Setup

This code was developed to run on the UCSF QB3 cluster. Running the code on a different system will require installation of the required software and adapting the code accordingly. This project was built using [Python](https://www.python.org/downloads/release/python-352/) and certain analysis steps requires [CellProfiler](http://cellprofiler.org/).

Further, the code expects the input data to be organized in particular folder structures (more details below). If copying the data to a new location, keep the internal folder structure unchanged.


## How to run analysis

This section details how to use the code to do the following
1. Generate processed and feature data from raw data (raw data -> processed/feature data:)
2. Generate figures from feature data (feature data -> paper figures)


### Processed/Feature Data

**Data required**: raw image data ("plates") 

We use image analysis to produce the cell/crypt segmentation results (processed data) and extract various image features. The analysis is done in three steps:

1. Pre-segmentation
    - Sets up processed data folder 
    - Sets up preseg jobs on the cluster
    - Rough segmentation using CellProfiler
2. Segmentation
    - Combines job array output from preseg step 
    - Identifies dense nuclear regions and generates masks for these regions 
    - Sets up merged segmentation jobs on the cluster
    - Final segmentation (mixed-strategy) using CellProfiler
3. Feature extraction
    - Combines job array output from merged step
     - Extracts "number" features from CellProfiler output 
    - Identifies crypts  
    - Measures crypt related features 
    - Saves features in csv file as both well- and crypt- level features

To run the analysis:

1. Setup paths
    
    Update the paths in `config/config.yml`, in the paths (top) section. 

    - root: path to top-level organoid_analysis folder
    - config: path to config folder
    - raw_data: path to raw data 
    - processed_data: path to folder for storing processed data (this is an empty new directory before running the code)
    - analysis_output: path to store any output files

    In general, if all folders are in the same place, simply change `/awlab/projects/2015_08_gut_organoid_analysis` to `path/to/project`, the top-level folder that contains `organoid_analysis` (code folder) and `plates` (raw data folder). Remember to create the `processed_data` and `analysis_output` folders.

2. Update more paths
    
    There are three scripts located in data_analysis:
    - 1_setup_preseg_jobs.py
    - 2_setup_mergedseg_jobs.py
    - 3_gen_crypt_masks.py

    Change the path(s) at the top of each file such that they correctly point to the code folder. As above, change `/awlab/projects/2015_08_gut_organoid_analysis` to `path/to/project`, the top-level folder that contains `organoid_analysis` (code folder)

3. Specify experiments
    
    Specify which experiments / plates to run in the `config.yml` file, in `exp_list` under `cp_sge_job`. For example, to run **ct20a** and **ct24_control1**:

    ```
    cp_sge_job: &cp_sge_job
      mem_free: 12G
      runtime: '12:00:00'
      start_job_num: 1
      num_jobs: 16
      exp_list: [ct20a, ct24_control1]
      done_file: 'done.txt
    ```

4. Run pre-segmentation 

    ```
    python 1_setup_preseg_jobs.py
    ```

    The output will be printed which is a string 'bash xxx'
    
    Paste the output string 'bash xxx' on the cluster. This will submit the pre-segmentation jobs. 

5. Run segmentation 

    ```
    python 2_setup_mergedseg_jobs.py
    ```

    The output will be printed which is a string 'bash xxx'
    
    Paste the output string 'bash xxx' on the cluster. This will submit the segmentation jobs. 

6. Extract features

    ```
    python 3_extract_features.py
    ```

    The features are stored in csv files in `processed_data/[plate]/merged/combined/[crypt_measures|well_measures]`

For simplicity, steps 4-6 are written together in `run_analysis.ipynb` 


### Publication Figures

From the feature data, we can generate the figures in the paper. All the plots are done in the notebook `data_analysis/paper_figures.ipynb`.

1. Setup paths

    Update the paths in the **"Initialize settings"** block of the notebook file such that processed_path, output_path, and program_path points to the processed data folder, the desired output folder to save graphs, and the code folder, respectively

    Also update **subdir** to point to the sub-folder in the output folder that graphs should be saved to.

2. Save plots 

    If **save_plots** is True, whenever a graph is generated, the graph and the associated data will be saved to the specified output folder. The graph will be saved as both a png and a svg file; the data will be saved as a csv file.

3. Run code blocks to generate graphs


## License

This project is licensed under the MIT License. See LICENSE.md for details.