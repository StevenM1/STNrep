### Analysis code for 
# `7T functional MRI finds no evidence for distinct functional subregions in the subthalamic nucleus during a speeded decision-making task`

###
The code is organized in ipython notebooks, and combined produce all results reported in the paper:

1. 0a_coregister_flash_to_t1w.ipynb
Coregisters the FLASH images to T1w space for dataset 1. In this dataset, a 0.5mm isotropic FLASH image was collected in a separate session from the task and T1w (MP2RAGE) data. The manually delineated STN masks were in FLASH space. In this notebook, we compute the warp between the FLASH to T1w space images.

2. 0b_warp_stn_masks.ipynb
Here, we warp the STN masks to T1w images of both datasets. For the Accolla masks, we warp to T1w space via the warp computed by fmriprep.

3. 0c_PCA_masks.ipynb
Here, we compute the STN segments based on a PCA that identifies the main dorsolateral-ventromedial axis of the STN. We also make plots to ensure the STN masks are in the correct location for each subject.

4. 0d_extract_signals.ipynb
Here, we extract signals from each STN subregion.

5. 1_Combine_events_confounds_signal.ipynb
This is a simple notebook that combines all relevant events, confounds, and signal to individual pandas DataFrames, which is easier for further ROI analysis.

6. 2_deconvolve_stn_responses.ipynb
This notebook does the signal deconvolutions based on Fourier basis set.

7. 3_Fit_wholebrain_GLM-FEAT.ipynb
This notebook does whole-brain GLM fitting. It creates .fsf-files for FEAT per individual run, then combines runs with flame (fixed effects), and then combines subjects with flame1.

8. 4_Analyze_wholebrain_GLM.ipynb
This notebook makes plots of the whole-brain GLM results.

9. 6a_ROI_GLM-PCA.ipynb
This notebook does ROI-wise GLMs based on the PCA-generated STN subdivisions (ie from notebook 0c)

10. 6b_ROI_GLM-Accolla.ipynb
This notebook does ROI-wise GLMs based on the Accolla et al. (2014) STN subdivisions

11. 7a_Interindividual_correlations-PCA.ipynb
This notebook does the interindividual correlation analyses between DDM parameters and BOLD signal changes in each STN subregion, based on the PCA-defined STN subregion masks

12. 7b_Interindividual_correlations-Accolla.ipynb
This notebook does the interindividual correlation analyses between DDM parameters and BOLD signal changes in each STN subregion, based on the STN subregion masks defined by Accolla et al. (2014)
