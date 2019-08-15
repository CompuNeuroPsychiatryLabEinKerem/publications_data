%% Before running this script - manual stages:
% 1. Compute ANOVA across all subjects with the six spatial scales conditions, to identify regions with scale sensitivity ("ANCOVA random effects analysis" dialog); Change threshold to q=0.01 (FDR-corrected); Save the results to a VMP file - "anova_results.vmp"

subjects_group_analysis_dir = 'F:\Scales_of_space_project\Analysis_group\';
vmp_filename = 'anova_results.vmp';
glm_filename = 'all_subjs_MNI_combined_no_control.glm';

% Reading the VMP and GLM using Neuroelf
vmp = xff(fullfile(subjects_group_analysis_dir, vmp_filename));
glm = xff(fullfile(subjects_group_analysis_dir, glm_filename));

% Getting the beta values from the GLM, and computing an average of all subjects' beta maps for each condition 
s = size(glm.GLMData.Subject(1).BetaMaps);
all_beta_maps = zeros(s(1), s(2), s(3), 6);
for s = 1:length(glm.GLMData.Subject)
    for i = 1:6
        all_beta_maps(:,:,:,i) = all_beta_maps(:,:,:,i) + glm.GLMData.Subject(s).BetaMaps(:,:,:,i);
    end
end
all_beta_maps = all_beta_maps / 6;

% masking with the anova results, where FDR threshold passes 0.01
for i = 1:6
    all_beta_maps(:,:,:,i) = all_beta_maps(:,:,:,i) .* (vmp.Map(1).VMPData >= vmp.Map(1).FDRThresholds(6,2));
end


% Computing maximum beta, and gaussian fit to betas graph
max_map = zeros(size(vmp.Map(1).VMPData));
gaussian_fit_rsquared_map = zeros(size(vmp.Map(1).VMPData));
gaussian_fit_a_amplitude_map = zeros(size(vmp.Map(1).VMPData));
gaussian_fit_b_center_map = zeros(size(vmp.Map(1).VMPData));
gaussian_fit_c_width_map = zeros(size(vmp.Map(1).VMPData));
% going over all voxels
for i = 1:size(all_beta_maps, 1)
    disp(i)
    for j = 1:size(all_beta_maps, 2)
        for k = 1:size(all_beta_maps, 3)
            % getting the betas for all conditions for this voxel
            curr_beta_vec = squeeze(all_beta_maps(i,j,k,:));
            % normalizing the betas so that they'll all have positive values
            curr_beta_vec_normalized = curr_beta_vec - min(curr_beta_vec);
            
            if sum(curr_beta_vec)~= 0   % making sure the cell is not empty
                % calculating the maximum beta
                [~, max_map(i,j,k)] = max(curr_beta_vec);
                
                % calculating a Gaussian fit to the normalized betas graph
                try
                    [f, gof] = fit([1:6]', curr_beta_vec_normalized, 'gauss1', 'Lower', [0 -100 0], 'Upper', [100 100 100]);
                    gaussian_fit_rsquared_map(i,j,k) = gof.rsquare;
                    gaussian_fit_a_amplitude_map(i,j,k) = f.a1;
                    gaussian_fit_b_center_map(i,j,k) = f.b1;
                    gaussian_fit_c_width_map(i,j,k) = f.c1;
                catch
                    gaussian_fit_rsquared_map(i,j,k) = -1;
                    gaussian_fit_a_amplitude_map(i,j,k) = nan;
                    gaussian_fit_b_center_map(i,j,k) = nan;
                    gaussian_fit_c_width_map(i,j,k) = nan;
                end
            end
            
            clear curr_beta_vec; clear curr_beta_vec_normalized; clear f; clear gof;
        end
    end
end

% Additional processing of Gaussian fit maps
% removing places with NaN values in R squared (unable to fit Gaussian)
gaussian_fit_a_amplitude_map(isnan(gaussian_fit_rsquared_map)) = nan;
gaussian_fit_b_center_map(isnan(gaussian_fit_rsquared_map)) = nan;
gaussian_fit_c_width_map(isnan(gaussian_fit_rsquared_map)) = nan;
% Changing Gaussian peak locations below minimum and maximum scale to these values
gaussian_fit_center_corrected = gaussian_fit_b_center_map;
gaussian_fit_center_corrected(gaussian_fit_center_corrected < 1 & gaussian_fit_center_corrected ~= 0) = 1; 
gaussian_fit_center_corrected(gaussian_fit_center_corrected > 6) = 6;
% Thresholding by minimum value of fit
gaussian_fit_center_rsquared_thresholded = gaussian_fit_center_corrected; gaussian_fit_center_rsquared_thresholded(gaussian_fit_rsquared_map < 0.7) = nan;

% Adding the results to the VMP
m = 1;     % the number of the last existing map in the VMP
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'betas_scale1'; vmp.Map(m).VMPData = all_beta_maps(:,:,:,1);
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'betas_scale2'; vmp.Map(m).VMPData = all_beta_maps(:,:,:,2);
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'betas_scale3'; vmp.Map(m).VMPData = all_beta_maps(:,:,:,3);
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'betas_scale4'; vmp.Map(m).VMPData = all_beta_maps(:,:,:,4);
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'betas_scale5'; vmp.Map(m).VMPData = all_beta_maps(:,:,:,5);
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'betas_scale6'; vmp.Map(m).VMPData = all_beta_maps(:,:,:,6);

m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'max_map'; vmp.Map(m).VMPData = max_map;

m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_rsquared_map'; vmp.Map(m).VMPData = gaussian_fit_rsquared_map;
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_a_amplitude_map'; vmp.Map(m).VMPData = gaussian_fit_a_amplitude_map;
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_b_center_map'; vmp.Map(m).VMPData = gaussian_fit_b_center_map;
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_c_width_map'; vmp.Map(m).VMPData = gaussian_fit_c_width_map;
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_corrected'; vmp.Map(m).VMPData = gaussian_fit_center_corrected;
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_rsquared_thresholded'; vmp.Map(m).VMPData = gaussian_fit_center_rsquared_thresholded;

m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'max_map_room'; vmp.Map(m).VMPData = (max_map == 1); vmp.Map(m).RGBUpperThreshPos = [255 0 0]; vmp.Map(m).RGBLowerThreshPos = [255 0 0];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'max_map_building'; vmp.Map(m).VMPData = (max_map == 2); vmp.Map(m).RGBUpperThreshPos = [255 167 25]; vmp.Map(m).RGBLowerThreshPos = [255 167 25];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'max_map_neighborhood'; vmp.Map(m).VMPData = (max_map == 3); vmp.Map(m).RGBUpperThreshPos = [255 255 0]; vmp.Map(m).RGBLowerThreshPos = [255 255 0];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'max_map_city'; vmp.Map(m).VMPData = (max_map == 4); vmp.Map(m).RGBUpperThreshPos = [85 255 0]; vmp.Map(m).RGBLowerThreshPos = [85 255 0];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'max_map_country'; vmp.Map(m).VMPData = (max_map == 5); vmp.Map(m).RGBUpperThreshPos = [0 150 255]; vmp.Map(m).RGBLowerThreshPos = [0 150 255];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'max_map_continent'; vmp.Map(m).VMPData = (max_map == 6); vmp.Map(m).RGBUpperThreshPos = [170 0 255]; vmp.Map(m).RGBLowerThreshPos = [170 0 255];

m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_rsquared_thresholded_room'; vmp.Map(m).VMPData = (round(gaussian_fit_center_rsquared_thresholded) == 1); vmp.Map(m).RGBUpperThreshPos = [255 0 0]; vmp.Map(m).RGBLowerThreshPos = [255 0 0];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_rsquared_thresholded_building'; vmp.Map(m).VMPData = (round(gaussian_fit_center_rsquared_thresholded) == 2); vmp.Map(m).RGBUpperThreshPos = [255 167 25]; vmp.Map(m).RGBLowerThreshPos = [255 167 25];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_rsquared_thresholded_neighborhood'; vmp.Map(m).VMPData = (round(gaussian_fit_center_rsquared_thresholded) == 3); vmp.Map(m).RGBUpperThreshPos = [255 255 0]; vmp.Map(m).RGBLowerThreshPos = [255 255 0];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_rsquared_thresholded_city'; vmp.Map(m).VMPData = (round(gaussian_fit_center_rsquared_thresholded) == 4); vmp.Map(m).RGBUpperThreshPos = [85 255 0]; vmp.Map(m).RGBLowerThreshPos = [85 255 0];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_rsquared_thresholded_country'; vmp.Map(m).VMPData = (round(gaussian_fit_center_rsquared_thresholded) == 5); vmp.Map(m).RGBUpperThreshPos = [0 150 255]; vmp.Map(m).RGBLowerThreshPos = [0 150 255];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_rsquared_thresholded_continent'; vmp.Map(m).VMPData = (round(gaussian_fit_center_rsquared_thresholded) == 6); vmp.Map(m).RGBUpperThreshPos = [170 0 255]; vmp.Map(m).RGBLowerThreshPos = [170 0 255];

% Setting additional map visualization parameters
for m = 2:7
    vmp.Map(m).LowerThreshold = 0;
    vmp.Map(m).UpperThreshold = 8;
end
for m = 8:length(vmp.Map)
    vmp.Map(m).LowerThreshold = 1;
    vmp.Map(m).UpperThreshold = 6;
end
for m = 9:length(vmp.Map)
    vmp.Map(m).UseRGBColor = 1;
end
vmp.NrOfMaps = length(vmp.Map);

% Saving the new VMP
vmp.SaveAs(fullfile(subjects_group_analysis_dir, 'all_subjs_betas_ANOVA_masked_with_parameters.vmp'));



%% Manual stage - save the VOIs of the required regions from the anova GLM results, by using the "convert map clusters to VOIs" function




%% Create table of activation clusters coordinates (at their center-of-mass) and their corresponding locations in the AICHA brain atlas 
% 
% Clusters were defined using the "convert map clusters to VOIs" option in brain voyager. 
% Each VOI file should contain all of the VOIs corresponding to a peak in that scale.
% (Clusters created from the Gaussian fit map for each scale separately)

% Reading the VOI files
voi_dir = 'F:\Scales_of_space_project\Analysis_group\';
voi_files = dir(fullfile(voi_dir, 'vois_gaussianfit_anova_*.voi'));


% identifying the VOIs centers of mass and number of voxels
VOIs_table = cell(length(voi_files), 1);
for f = 1:length(voi_files)
    current_voi = xff(fullfile(voi_dir, voi_files(f).name));
    VOIs_table{f} = zeros(current_voi.NrOfVOIs, 4);
    for i=1:current_voi.NrOfVOIs
        VOIs_table{f}(i, 1:3) = mean(current_voi.VOI(i).Voxels);    % Averaging the VOI values for X,Y,Z separately, to get the center of mass coordinate
        VOIs_table{f}(i, 4) = current_voi.VOI(i).NrOfVoxels / 27;   % Getting the number of voxels in the VOI - dividing by 27 since voxels as defined in the VOI are in anatomical resolution (1mm^3), in contrast to functional resolution (3mm^3)
    end
end
% rounding the values
for f = 1:length(voi_files)
    VOIs_table{f} = round(VOIs_table{f}); 
end


% finding each VOI's correspondence to regions in the AICHA atlas
% reading the AICHA atlas from file using SPM
atlas_filename = 'C:\templates\AICHA.nii';
atlas = spm_vol(atlas_filename);
atlas_mat = atlas.mat;
atlas = spm_read_vols(atlas);
% getting the area names
atlas_regions_list = 'C:\templates\AICHA_ROI_MNI_V1.txt';
fid = fopen(atlas_regions_list); 
atlas_area_list = textscan(fid, '%s %s %s'); 
fclose(fid);
atlas_area_list = atlas_area_list{1};
% identifying for each VOI its corresponding AICHA region
region_names = cell(size(VOIs_table));
for f = 1:length(VOIs_table)
    region_names{f} = cell(size(VOIs_table{f}, 1), 1);
    for i=1:size(VOIs_table{f}, 1)
        current_region_coordinates = inv(atlas_mat) * [VOIs_table{f}(i,1:3) 1]';
        current_region_coordinates = round(current_region_coordinates(1:3));
        current_region_number = atlas(current_region_coordinates(1), current_region_coordinates(2), current_region_coordinates(3));
        if current_region_number~=0
            region_names{f}{i} = atlas_area_list{current_region_number};
        else
            region_names{f}{i} = 0;
        end
    end
end


%% checking the order of activations along the Y axis of the brain, inside each region of interest
% Manual - combine VOIs of all scales from each region, and convert them into MSK files via BrainVoyager

vmp = xff(fullfile(subjects_group_analysis_dir, 'all_subjs_betas_ANOVA_masked_with_parameters.vmp'));
glm = xff(fullfile(subjects_group_analysis_dir, 'all_subjs_MNI_combined_no_control.glm'));

% reading the VOI mask (choose among the following options)
current_VOI_mask = xff(fullfile(subjects_group_analysis_dir, 'Hippocampal_mask_HO_atlas.nii')); current_VOI_mask = double(current_VOI_mask.VoxelData);  % Hippocampus
% current_VOI_mask = xff(fullfile(subjects_group_analysis_dir, 'lateral_parietal_new.msk')); current_VOI_mask = current_VOI_mask.Mask;  % Lateral parietal gradient
% current_VOI_mask = xff(fullfile(subjects_group_analysis_dir, 'medial_temporal_new.msk')); current_VOI_mask = current_VOI_mask.Mask;   % Medial temporal gradient
% current_VOI_mask = xff(fullfile(subjects_group_analysis_dir, 'medial_parietal_new.msk')); current_VOI_mask = current_VOI_mask.Mask;   % Medial parietal gradient

num_scales = 6;
num_subjects = length(glm.GLMData.Subject);
num_coordinates = size(current_VOI_mask, 1);    % The coordinates along the posterior-anterior axis
sa = size(current_VOI_mask);

% Reading each subject's GLM data
beta_maps_all_subjects = cell(1, num_scales);
for i = 1:num_scales, beta_maps_all_subjects{i} = zeros(sa(1),sa(2),sa(3),num_subjects); end
for s = 1:num_subjects
    for i = 1:num_scales
        beta_maps_all_subjects{i}(:,:,:,s) = glm.GLMData.Subject(s).BetaMaps(:,:,:,i);
    end
end


% Calculating the average beta per subject and scale
beta_vecs_coordinates_all_subjects = nan(num_subjects, num_scales, num_coordinates);
for m = 1:num_coordinates
    temp_mask = zeros(size(current_VOI_mask));
    temp_mask(m, :, :) = 1;
    if length(find(current_VOI_mask==1 & temp_mask==1))>10
        for s = 1:num_subjects
            curr_beta_vec = nan(1,num_scales);
            for i = 1:num_scales
                temp_beta = beta_maps_all_subjects{i}(:,:,:,s);
                curr_beta_vec(i) = nanmean(temp_beta(find(current_VOI_mask==1 & temp_mask==1)));
            end
            beta_vecs_coordinates_all_subjects(s, :, m) = curr_beta_vec;
        end
    end
end

% Calculating the Gaussian peak and max value for each coordinate in each subject
max_values_per_subject = nan(num_subjects, num_coordinates);
gaussian_peak_per_subject = nan(num_subjects, num_coordinates);
gaussian_fit_per_subject = nan(num_subjects, num_coordinates);
for m = 1:num_coordinates
    % Per subject
    for s = 1:num_subjects
        % getting the betas for all conditions for this voxel
        curr_beta_vec = beta_vecs_coordinates_all_subjects(s, :, m);
        % normalizing the betas so that they'll all have positive values
        curr_beta_vec_normalized = curr_beta_vec - min(curr_beta_vec);
        
        if sum(isnan(curr_beta_vec))~=6
            % Maximally active scale for each subject
            [~, max_values_per_subject(s, m)] = max(curr_beta_vec);
            
            % calculating a Gaussian fit to the normalized betas graph
            try
                [f, gof] = fit([1:6]', curr_beta_vec_normalized', 'gauss1', 'Lower', [0 -100 0], 'Upper', [100 100 100]);
                gaussian_peak_per_subject(s, m) = f.b1;
                gaussian_fit_per_subject(s, m) = gof.rsquare;
            catch
                gaussian_peak_per_subject(s, m) = nan;
                gaussian_fit_per_subject(s, m) = nan;
            end
        end
    end
end

% Changing Gaussian peak locations below minimum and maximum scale to these values
gaussian_peak_per_subject_corrected = gaussian_peak_per_subject;
gaussian_peak_per_subject_corrected(gaussian_peak_per_subject_corrected < 1 & gaussian_peak_per_subject_corrected ~= 0) = 1; 
gaussian_peak_per_subject_corrected(gaussian_peak_per_subject_corrected > 6) = 6;

% Plotting the result
figure; plot(mean(gaussian_peak_per_subject_corrected))
hold on; plot(mean(max_values_per_subject),'r')
% changing Y axis values to scales
Y_labels = {'Room','Building','Neighborhood','City','Country','Continent'};
ylim([1 6]), yticks(1:6), yticklabels(Y_labels)
% changing minimum and maximum displayed X values
locations_with_data = find(~isnan(a)); 
xlim([locations_with_data(1) locations_with_data(end)])
% changing X axis values to MNI Y coordinates
x_start = 128 - vmp.XStart; x_end = 128 - vmp.XEnd;
mni_coordinates = (x_start:-3:x_end) - 1;   % the coordinate at the center of each voxel
xticks(1:60), xticklabels(mni_coordinates), set(gca, 'xdir', 'reverse')




scale_selectivity_values_gaussian = mean(gaussian_peak_per_subject_corrected);
scale_selectivity_values_gaussian = scale_selectivity_values_gaussian(~isnan(scale_selectivity_values_gaussian));
true_slope_gaussian = fit((1:length(scale_selectivity_values_gaussian))', scale_selectivity_values_gaussian', 'poly1');  % linear fit to the actual data
true_slope_gaussian = true_slope_gaussian.p1;     % The slope of the linear fit
% Shuffling the original data 1000 time and calculating linear fit for each shuffle
num_iters = 1000;
random_slopes = zeros(1, length(num_iters));
for i=1:1000
    temp_fit = fit((1:length(scale_selectivity_values_gaussian))', scale_selectivity_values_gaussian(randperm(length(scale_selectivity_values_gaussian)))', 'poly1');
    random_slopes(i) = temp_fit.p1;
end
permutation_test_p_value_gaussian = sum(random_slopes <= true_slope_gaussian) / num_iters;

scale_selectivity_values_max = mean(max_values_per_subject);
scale_selectivity_values_max = scale_selectivity_values_max(~isnan(scale_selectivity_values_max));
true_slope_max = fit((1:length(scale_selectivity_values_max))', scale_selectivity_values_max', 'poly1');  % linear fit to the actual data
true_slope_max = true_slope_max.p1;     % The slope of the linear fit
% Shuffling the original data 1000 time and calculating linear fit for each shuffle
num_iters = 1000;
random_slopes = zeros(1, length(num_iters));
for i=1:1000
    temp_fit = fit((1:length(scale_selectivity_values_max))', scale_selectivity_values_max(randperm(length(scale_selectivity_values_max)))', 'poly1');
    random_slopes(i) = temp_fit.p1;
end
permutation_test_p_value_max = sum(random_slopes <= true_slope_max) / num_iters;


% Individual subject level analysis
for s = 1:num_subjects
    disp(s);
    
    scale_selectivity_values_gaussian = gaussian_peak_per_subject_corrected(s, :);
    scale_selectivity_values_gaussian = scale_selectivity_values_gaussian(~isnan(scale_selectivity_values_gaussian));
    true_slope_gaussian = fit((1:length(scale_selectivity_values_gaussian))', scale_selectivity_values_gaussian', 'poly1');  % linear fit to the actual data
    true_slope_gaussian = true_slope_gaussian.p1;     % The slope of the linear fit
    % Shuffling the original data 1000 time and calculating linear fit for each shuffle
    num_iters = 1000;
    random_slopes = zeros(1, length(num_iters));
    for i=1:1000
        temp_fit = fit((1:length(scale_selectivity_values_gaussian))', scale_selectivity_values_gaussian(randperm(length(scale_selectivity_values_gaussian)))', 'poly1');
        random_slopes(i) = temp_fit.p1;
    end
    permutation_test_p_value_gaussian = sum(random_slopes <= true_slope_gaussian) / num_iters;
    
    scale_selectivity_values_max = max_values_per_subject(s, :);
    scale_selectivity_values_max = scale_selectivity_values_max(~isnan(scale_selectivity_values_max));
    true_slope_max = fit((1:length(scale_selectivity_values_max))', scale_selectivity_values_max', 'poly1');  % linear fit to the actual data
    true_slope_max = true_slope_max.p1;     % The slope of the linear fit
    % Shuffling the original data 1000 time and calculating linear fit for each shuffle
    num_iters = 1000;
    random_slopes = zeros(1, length(num_iters));
    for i=1:1000
        temp_fit = fit((1:length(scale_selectivity_values_max))', scale_selectivity_values_max(randperm(length(scale_selectivity_values_max)))', 'poly1');
        random_slopes(i) = temp_fit.p1;
    end
    permutation_test_p_value_max = sum(random_slopes <= true_slope_max) / num_iters;
    
    all_slopes_gaussian(s) = true_slope_gaussian;
    all_slopes_max(s) = true_slope_max;
    all_pvalues_gaussian(s) = permutation_test_p_value_gaussian;
    all_pvalues_max(s) = permutation_test_p_value_max;
    
end



%% Analyze relation to resting-state networks from the Yeo et al., Neuropyhisiology 2011 parcellation
masks_dir = 'F:\Rescaling_analysis\1_Rescaling_analysis_paper_final\Analysis_group\VOI_masks';
files = dir(fullfile(masks_dir, '*.nii'));
% reading and reslicing Yeo's networks
Yeo_nets = reslice_data('C:\Users\michaelpeer1\Dropbox\Michael_templates\rYeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii', fullfile(masks_dir, files(1).name), 0);
% calculating number of voxels from each VOI in each network
all_nets_files = zeros(length(files), 7);
for i = 1:length(files)
    current_file = spm_read_vols(spm_vol(fullfile(masks_dir, files(i).name)));
    for j = 1:7
        all_nets_files(i, j) = sum(current_file(:) == 1 & Yeo_nets(:) == j);
    end
end


