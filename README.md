This repository contains the custom code used in 
[Authors] et al., [name] 2025 (Manuscript in preparation).

## Overview
This codebase includes custom code to preprocess and analyze Hi-C data, as well as supplementary files that accompany these analyses. The repository is organized according to different analytic steps, numbered from 0 to N. Each step represents a different aspect of the Hi-C data processing pipeline.

## Dependencies
To run the scripts provided in this repository, create a conda environment from the yml file. This ensures you have the correct Python packages and system libraries installed.

1. Install Miniconda or Anaconda if you have not done so already.
2. Clone this repository to your local machine.
3. Create the conda environment:
   ```
   conda env create -f environment_hic_1.yml
   ```
4. Activate the environment:
   ```
   conda activate <environment-name-from-yml>
   ```
5. Use `environment_hic_1.yml` for:
    **0_prepare_maps**
    **1_general_chromatin_properties**  
    **2_tads**  
    **3_chromatin_loops** 


## Scripts Overview
Below is a brief description of each script folder in this repository:

1. **0_prepare_maps**  
   Contains scripts and notebooks for preparing raw Hi-C contact maps (e.g., normalizing contacts, generating matrix files).

2. **1_general_chromatin_properties**  
   Scripts for computing genome-wide and chromosome-level properties (e.g., contact probability decay, compartment analysis).

3. **2_tads**  
   Scripts to identify and analyze Topologically Associating Domains (TADs) or other domain-like architectures.

4. **3_chromatin_loops**  
   Scripts to detect and annotate loops or point interactions in Hi-C data (e.g., using peak-calling algorithms).

5. **4_additional_files**  
   Contains supporting files necessary for the analyses (such as annotation files, reference genome files, or any custom data used by the scripts).

## Usage Instructions
1. Ensure that you have prepared your Hi-C data in the required format (e.g., contact matrix files, or raw data that can be handled by the scripts in `0_prepare_maps`).
2. Run the scripts in the logical order:
   1. Navigate to `0_prepare_maps` and run the corresponding scripts to generate or normalize Hi-C maps.
   2. Use `1_general_chromatin_properties` scripts to analyze broader-scale chromatin features.
   3. Use `2_tads` scripts to detect and analyze TADs.
   4. Use `3_chromatin_loops` scripts to detect and analyze loops.
   5. Access additional or supporting data in `4_additional_files`.

