Representation of different spatial scales in the human brain:

Despite a large body of research on spatial representations, a fundamental question remains unresolved - whether the brain uses similar 
mechanisms to represent environments at different scales. In our study (Michael Peer, Yorai Ron, Rotem Monsa and Shahar Arzy, eLife 2019) 
we describe functional gradients of spatial processing organized according to spatial scales.
As part of this publication, our statistical maps and results as well as the full analysis scripts are available below. 
for questions please contact michael.peer (at) mail.huji.ac.il . 




Experimental results:

The paper's statistical maps (in MNI space), representing cortical preference to processing of different spatial scales, can be found 
here: 

1. spatial_scales_gaussian_fit.nii - preference to spatial scale, as determined by the peak of a Gaussian fit to each voxel's beta graph 
(Figure 2in the paper).

2. spatial_scales_maximal_activity.nii - preference to spatial scale, as determined by the maximal beta in each voxel (Figure 2=-Figure 
Supplement 3 in the paper).

3. spatial_scales_results.nii - separate maps for the previous two measures as well as intermediate analysis stages. This file contains 
the following statistical maps:
- map 1 - results of ANOVA analysis within each voxel, to identify voxels sensitive to changes in spatial scale
- maps 2-7 - the beta values for each condition (room, building, neighborhood, city, country, continent) in each voxel, within the voxels identified by the ANOVA analysis
- map 8 - winner map indicating the scale with maximal activity (within the voxels identified in the ANOVA analysis)
- map 9 - R^2 of Gaussian fit to the voxels' beta graphs (within the voxels identified in the ANOVA analysis)
- map 10 - amplitude of Gaussian fit to the voxels' beta graphs (within the voxels identified in the ANOVA analysis)
- map 11 - location of peak of Gaussian fit to the voxels' beta graphs (within the voxels identified in the ANOVA analysis)
- map 12 - width of Gaussian fit to the voxels' beta graphs (within the voxels identified in the ANOVA analysis)
- map 13 - corrected location of peak of Gaussian fit - peaks below 1 corrected to 1, peaks above 6 corrected to 6
- map 14 - corrected location of peak of Gaussian fit (map 13), thresholded to show voxels with r^2 >= 0.7
- maps 15-20 - preferred scale of each voxel that was identified using the ANOVA analysis, calculated by the maximal beta value within the voxel (maps: voxels preferentially active for room, building, neighborhood, city, country and continent).
- maps 21-26 - preferred scale of each voxel that was identified using the ANOVA analysis, calculated by the location of the Gaussian fit peak (maps: voxels preferentially active for room, building, neighborhood, city, country and continent).



Analysis software:

Our full analysis codes, used to identify brain activity related to spatial representation at different scales, are provided here. 
These scripts use Matlab (tested on Matlab 2017B), BrainVoyager (tested on BrainVoyager 20.6) and Neuroelf (tested on NeuroElf v1.1).
The zip file contains the following scripts: 
- 1_Preprocessing_pipeline.m - pre-processing pipeline for the fMRI data using BrainVoyager and NeuroElf.
- 2_GLM_estimation.m - estimation of BOLD response to each experimental condition using BrainVoyager and NeuroElf.
- 3_Analyze_GLM_results.m - analysis of the beta weights from the GLM to identify spatial scale selectivity in each voxel, analysis 
of the order of activations along the posterior-anterior axis, and getting the coordinates of each activity cluster.
- 4_analyzing_effects_of_correlated_factors.m - using parametric modulation to measure the effects of possible factors on the results, 
and comparison of activations to a lexical control condition.
- create_sdm_with_other_factors.m - helper function to create design matrices with the parametrically modulated conditions according to 
the different contributing factors, used inside the previous script.

