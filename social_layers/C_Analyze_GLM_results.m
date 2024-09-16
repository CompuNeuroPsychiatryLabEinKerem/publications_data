%% Before running this script - manual stages:
% 1. Compute ANOVA across all subjects with the five social scales conditions, 
% to identify regions with scale sensitivity ("ANCOVA random effects analysis" dialog); 
% Change threshold to q=0.01 (FDR-corrected); Save the results to a VMP file - "anova_results.vmp"

subjects_group_analysis_dir = 'social_main\3_group_analysis\';
output_dir = 'social_main\3_group_analysis\ANOVA_q0.01_thr1000';
vmp_filename = 'anova_results_q0.01.vmp';  % it's only for the ANOVA results
glm_filename = 'all_subjs_MNI_combined_no_control.glm';

% Reading the VMP and GLM using Neuroelf
vmp = xff(fullfile(subjects_group_analysis_dir, vmp_filename));
glm = xff(fullfile(subjects_group_analysis_dir, glm_filename));

msk_dir = 'social_main\3_group_analysis\VOIs';
msk_filename = 'anova_results_all_clusters_q0.01_thr1000.msk';
VOI_mask = xff(fullfile(msk_dir, msk_filename)); VOI_mask = VOI_mask.Mask;

% Getting the beta values from the GLM, and computing an average of all subjects' beta maps for each condition
s = size(glm.GLMData.Subject(1).BetaMaps);
space_size = s(1:3);
subject_num = size(glm.GLMData.Subject, 2);
scales_num = 5;
all_beta_maps = zeros([subject_num, scales_num, space_size]);
for s = 1:subject_num  % subject
    for c = 1:scales_num  % scale
        % masking with whole brain regions base on anova results
        all_beta_maps(s,c,:,:,:) = glm.GLMData.Subject(s).BetaMaps(:,:,:,c) .* (VOI_mask == 1);
    end
end

% Compute maximum beta
max_map_mean = squeeze(mean(all_beta_maps));
[max_map_betas, max_map] = max(max_map_mean);
max_map_betas = squeeze(max_map_betas);
max_map = squeeze(max_map);
max_map(max_map_betas == 0) = 0;  % if max_beta == 0, make scale in this voxel = 0

% Compute gaussian fit to betas graph
gf_rsquared_map_all = zeros([subject_num, space_size]);
gf_a_amplitude_map_all = zeros([subject_num, space_size]);
gf_b_center_map_all = zeros([subject_num, space_size]);
gf_c_width_map_all = zeros([subject_num, space_size]);

for s = 1:subject_num  % subject
    % going over all voxels    
    for i = 1:space_size(1)
        for j = 1:space_size(2)
            for k = 1:space_size(3)
                
                % get betas for all scales
                curr_vox_beta_vec = squeeze(all_beta_maps(s,:,i,j,k));
                
                % normalizing the betas so that they'll all have positive values
                curr_vox_vec_normalized = curr_vox_beta_vec - min(curr_vox_beta_vec);
                
                if sum(curr_vox_beta_vec) ~= 0  % making sure the cell is not empty

                    % calculating a Gaussian fit to the normalized betas graph
                    try
                        [f, gof] = fit([1:5]', curr_vox_vec_normalized', 'gauss1', 'Lower', [0 -100 0], 'Upper', [100 100 100]);
                        gf_rsquared_map_all(s,i,j,k) = gof.rsquare;
                        gf_a_amplitude_map_all(s,i,j,k) = f.a1;
                        gf_b_center_map_all(s,i,j,k) = f.b1;
                        gf_c_width_map_all(s,i,j,k) = f.c1;
                    catch
                        gf_rsquared_map_all(s,i,j,k) = -1;
                        gf_a_amplitude_map_all(s,i,j,k) = nan;
                        gf_b_center_map_all(s,i,j,k) = nan;
                        gf_c_width_map_all(s,i,j,k) = nan;
                    end
                end
            
                clear curr_vox_beta_vec; clear curr_vox_vec_normalized; clear f; clear gof;
        
            end
        end
    end
end

% Additional processing of Gaussian fit maps
% removing places with NaN values in R squared (unable to fit Gaussian)
gf_a_amplitude_map_all(isnan(gf_rsquared_map_all)) = nan;
gf_b_center_map_all(isnan(gf_rsquared_map_all)) = nan;
gf_c_width_map_all(isnan(gf_rsquared_map_all)) = nan;

% Thresholding, averaging and removing outliers
gf_center_rsquared_thr08 = gf_b_center_map_all; gf_center_rsquared_thr08(gf_rsquared_map_all < 0.8) = nan;

% correcting center of GF below minimum and maximum scale
gf_center_rsquared_thr08(gf_center_rsquared_thr08 < 1) = 1; gf_center_rsquared_thr08(gf_center_rsquared_thr08 > 5) = 5;

gf_center_mean_rsquared_thr08 = squeeze(mean(gf_center_rsquared_thr08, 1, 'omitnan'));
gf_center_std_rsquared_thr08 = squeeze(std(gf_center_rsquared_thr08, 1, 'omitnan'));

% remove outliers
upper_bounds = gf_center_mean_rsquared_thr08 + 3*gf_center_std_rsquared_thr08;
lower_bounds = gf_center_mean_rsquared_thr08 - 3*gf_center_std_rsquared_thr08;
for i = 1:size(gf_center_rsquared_thr08, 1)
    outliers_idx = gf_center_rsquared_thr08(i, :) > upper_bounds | gf_center_rsquared_thr08(i, :) < lower_bounds;
    gf_center_rsquared_thr08(i, outliers_idx) = nan;
end
gf_center_mean_rsquared_thr08 = squeeze(mean(gf_center_rsquared_thr08, 1, 'omitnan'));

% perform t-test in each voxel with a "max map" result, comparing the chosen scale to the average of others.
mm_selec_hs = zeros(space_size);
mm_selec_pvals = zeros(space_size);

for i = 1:space_size(1)
    for j = 1:space_size(2)
        for k = 1:space_size(3)
            
            % get betas for all scales
            curr_vox_betas_mat = squeeze(all_beta_maps(:,:,i,j,k));
            
            if sum(curr_vox_betas_mat, 'all') ~= 0  % making sure the matrix is not empty

                % get chosen scale and betas
                chosen_scale = max_map(i, j, k);
                betas_chosen = squeeze(curr_vox_betas_mat(:, chosen_scale));
                curr_vox_betas_mat(:,chosen_scale) = []; % remove from matrix
                betas_nonchosen = mean(curr_vox_betas_mat, 2);

                % perform paired-sample t-test
                [h,p] = ttest(betas_chosen,betas_nonchosen);
                mm_selec_hs(i,j,k) = h;
                mm_selec_pvals(i,j,k) = p;
                
            end
            
            clear curr_vox_betas_mat; clear chosen_scale; 
            clear betas_chosen; clear betas_nonchosen; clear h; clear p;
    
        end
    end
end

% correct for multiple comparisons
mm_selec_pvals_vector = squeeze(reshape(mm_selec_pvals, 1, []));
mm_selec_pvals_vector(mm_selec_pvals_vector == 0) = nan;
mm_selec_fdrs = mafdr(mm_selec_pvals_vector, 'BHFDR', true);
mm_selec_fdrs = reshape(mm_selec_fdrs, space_size);

% filter max map according to selectivity
max_map_selthr05 = max_map .* (mm_selec_fdrs < 0.05);

gf_rsquared_map_all_mean = squeeze(mean(gf_rsquared_map_all, 1, 'omitnan'));
gf_a_amplitude_map_all_mean = squeeze(mean(gf_a_amplitude_map_all, 1, 'omitnan'));
gf_b_center_map_all_mean = squeeze(mean(gf_b_center_map_all, 1, 'omitnan'));
gf_c_width_map_all_mean = squeeze(mean(gf_c_width_map_all, 1, 'omitnan'));

% Adding the results to the VMP
m = 1;     % the number of the last existing map in the VMP
n = n + 1; vmp.Map(n) = vmp.Map(1); vmp.Map(n).Name = 'max_map_selthr05'; vmp.Map(n).VMPData = max_map_selthr05;

m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_rsquared_map'; vmp.Map(m).VMPData = gf_rsquared_map_all_mean;
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_a_amplitude_map'; vmp.Map(m).VMPData = gf_a_amplitude_map_all_mean;
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_b_center_map'; vmp.Map(m).VMPData = gf_b_center_map_all_mean;
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_c_width_map'; vmp.Map(m).VMPData = gf_c_width_map_all_mean;

m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_r2_thr0.8'; vmp.Map(m).VMPData = gf_center_mean_rsquared_thr08;

n = n + 1; vmp.Map(n) = vmp.Map(1); vmp.Map(n).Name = 'max_map_selthr05_scale1'; vmp.Map(n).VMPData = (max_map_selthr05 == 1); vmp.Map(n).RGBUpperThreshPos = [255 0 0]; vmp.Map(n).RGBLowerThreshPos = [255 0 0];
n = n + 1; vmp.Map(n) = vmp.Map(1); vmp.Map(n).Name = 'max_map_selthr05_scale2'; vmp.Map(n).VMPData = (max_map_selthr05 == 2); vmp.Map(n).RGBUpperThreshPos = [255 167 25]; vmp.Map(n).RGBLowerThreshPos = [255 167 25];
n = n + 1; vmp.Map(n) = vmp.Map(1); vmp.Map(n).Name = 'max_map_selthr05_scale3'; vmp.Map(n).VMPData = (max_map_selthr05 == 3); vmp.Map(n).RGBUpperThreshPos = [255 255 0]; vmp.Map(n).RGBLowerThreshPos = [255 255 0];
n = n + 1; vmp.Map(n) = vmp.Map(1); vmp.Map(n).Name = 'max_map_selthr05_scale4'; vmp.Map(n).VMPData = (max_map_selthr05 == 4); vmp.Map(n).RGBUpperThreshPos = [85 255 0]; vmp.Map(n).RGBLowerThreshPos = [85 255 0];
n = n + 1; vmp.Map(n) = vmp.Map(1); vmp.Map(n).Name = 'max_map_selthr05_scale5'; vmp.Map(n).VMPData = (max_map_selthr05 == 5); vmp.Map(n).RGBUpperThreshPos = [0 150 255]; vmp.Map(n).RGBLowerThreshPos = [0 150 255];

m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_r2_thr0.8_scale1'; vmp.Map(m).VMPData = (round(gf_center_mean_rsquared_thr08) == 1); vmp.Map(m).RGBUpperThreshPos = [255 0 0]; vmp.Map(m).RGBLowerThreshPos = [255 0 0];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_r2_thr0.8_scale2'; vmp.Map(m).VMPData = (round(gf_center_mean_rsquared_thr08) == 2); vmp.Map(m).RGBUpperThreshPos = [255 167 25]; vmp.Map(m).RGBLowerThreshPos = [255 167 25];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_r2_thr0.8_scale3'; vmp.Map(m).VMPData = (round(gf_center_mean_rsquared_thr08) == 3); vmp.Map(m).RGBUpperThreshPos = [255 255 0]; vmp.Map(m).RGBLowerThreshPos = [255 255 0];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_r2_thr0.8_scale4'; vmp.Map(m).VMPData = (round(gf_center_mean_rsquared_thr08) == 4); vmp.Map(m).RGBUpperThreshPos = [85 255 0]; vmp.Map(m).RGBLowerThreshPos = [85 255 0];
m = m + 1; vmp.Map(m) = vmp.Map(1); vmp.Map(m).Name = 'gaussian_fit_center_r2_thr0.8_scale5'; vmp.Map(m).VMPData = (round(gf_center_mean_rsquared_thr08) == 5); vmp.Map(m).RGBUpperThreshPos = [0 150 255]; vmp.Map(m).RGBLowerThreshPos = [0 150 255];

% Setting additional map visualization parameters
for m = 1:length(vmp.Map)
    vmp.Map(m).LowerThreshold = 1;
    vmp.Map(m).UpperThreshold = 5;
end
for m = 2:length(vmp.Map)
    vmp.Map(m).UseRGBColor = 1;
end
vmp.NrOfMaps = length(vmp.Map);

% Saving vmp & nii
vmp.SaveAs(fullfile(output_dir, 'ANOVA_q0.01_thr1000_GF_MM_results.vmp'));
vmp.ExportNifti(fullfile(output_dir, 'ANOVA_q0.01_thr1000_GF_MM_results.nii'), true);


%% social-spatial

spatial_max_map = load('social_main\spatial_scales_data\group_analysis\01_q0.01_thr1000_scales_corrected\max_map.mat');
spatial_max_map = spatial_max_map.max_map;

% calculating Jaccard index for each scale-pairs
jac_mat = nan(6, 5);
for soc_s = 1:5
    soc_s_mask = max_map_selthr05 == soc_s;
    for spat_s = 1:6
        spat_s_mask = spatial_max_map == spat_s;
        intersection = sum(soc_s_mask .* spat_s_mask, 'all');
        union = sum((soc_s_mask + spat_s_mask) > 0, 'all');
        jac_mat(spat_s, soc_s) = intersection/union;
    end
end

% permutation test: shuffle array and calculate Jaccard indices
num_perm = 100000;
perm_jac_mat = nan(6, 5, num_perm);

for p = (1:num_perm)
    
    % shuffle social scales (only in the social scales region)
    soc_mm_vals = nonzeros(max_map_selthr05(:));
    soc_mm_vals_shuf = soc_mm_vals(randperm(length(soc_mm_vals)));
    
    % put the social scale voxels back into array
    soc_mm_array_shuf = zeros(space_size);
    count_vox = 1;
    
    for i = 1:space_size(1)
        for j = 1:space_size(2)
            for k = 1:space_size(3)
                
                if max_map_selthr05(i, j, k) ~= 0
    
                    soc_mm_array_shuf(i, j, k) = soc_mm_vals_shuf(count_vox);
                    count_vox = count_vox + 1;
                    
                end
                
            end
        end
    end
    
    % calculate Jaccard index for each scale-pairs
    for soc_s = 1:5
        soc_s_mask = soc_mm_array_shuf == soc_s;
        for spat_s = 1:6
            spat_s_mask = spatial_max_map == spat_s;
            intersection = sum(soc_s_mask .* spat_s_mask, 'all');
            union = sum((soc_s_mask + spat_s_mask) > 0, 'all');
            perm_jac_mat(spat_s, soc_s, p) = intersection/union;
        end
    end

end

% calculate p-vals from permutation test
jac_p_vals = nan(6, 5);
for soc_s = 1:5
    for spat_s = 1:6
        jac_ind = jac_mat(spat_s, soc_s);
        perm_jac_inds = sort(squeeze(perm_jac_mat(spat_s, soc_s, :)));
        jac_p_vals(spat_s, soc_s) = sum(perm_jac_inds >= jac_ind, 'all')/num_perm;
    end
end

% correct for multiple comparisons
jac_p_vals_vector = squeeze(reshape(jac_p_vals, 1, []));
jac_fdrs = mafdr(jac_p_vals_vector, 'BHFDR', true);
jac_fdrs = reshape(jac_fdrs, spat_s, soc_s);


%% ROI analysis using Yeo's parcellations

% load Yeo's 7 networks and a combination of them
main_dir = 'social_main\3_group_analysis\VOIs\Yeo_ROI\all_7_nets';
nii_pattern = '*.nii.gz';
nii_files = dir(fullfile(main_dir, nii_pattern));
yeo_7_nets = nii_files(1:7);
yeo_surface = nii_files(8);

% load surface mask and transform into BV coordinates
n = neuroelf;
vmp = n.importvmpfromspms(fullfile(yeo_surface.folder, yeo_surface.name), [], [57, 49, 53; 236, 180, 202]);
surf_mask = round(double(vmp.Map.VMPData));

% filter social scale selective voxels to keep only surface
max_map_selthr05_surf = max_map_selthr05 .* surf_mask;
max_map_selthr05_surf_allscales = max_map_selthr05_surf > 0;

% calculate Jaccard of all scales combined vs. each network
jac_mat_yeo = nan(7, 1);
for f = 1:length(yeo_7_nets)
    cur_roi_fname = yeo_7_nets(f).name;
    vmp = n.importvmpfromspms(fullfile(nii_files(f).folder, cur_roi_fname), [], [57, 49, 53; 236, 180, 202]);
    net_mask = round(double(vmp.Map.VMPData));
    intersection = sum(max_map_selthr05_surf_allscales .* net_mask, 'all');
    union = sum((max_map_selthr05_surf_allscales + net_mask) > 0, 'all');
    jac_mat_yeo(f) = intersection/union;
end

% Permutations: shuffle scales in surface and calculate Jaccards
num_perm = 10000;
perm_jac_mat_yeo = nan(7, num_perm);

% number of scale selective voxels on surface
soc_surf_voxles_num = sum(max_map_selthr05_surf_allscales, "all");

% return locations of surface voxels
[is,js,ks] = ind2sub(size(surf_mask),find(surf_mask));
surf_voxels = [is,js,ks];

for p = (1:num_perm)

    % randomly choose surface voxels
    rand_mm_mask = zeros(size(surf_mask));
    rand_vox_inds = randsample(size(surf_voxels, 1), soc_surf_voxles_num);
    for vox = 1:size(rand_vox_inds, 1)
        rand_vox_ind = rand_vox_inds(vox);
        chosen_vox = surf_voxels(rand_vox_ind, :);
        rand_mm_mask(chosen_vox(1), chosen_vox(2), chosen_vox(3)) = 1;
    end
    
    % calculate Jaccard index vs. each network
    for f = 1:length(yeo_7_nets)
        cur_roi_fname = yeo_7_nets(f).name;
        vmp = n.importvmpfromspms(fullfile(nii_files(f).folder, cur_roi_fname), [], [57, 49, 53; 236, 180, 202]);
        net_mask = round(double(vmp.Map.VMPData));
        intersection = sum(rand_mm_mask .* net_mask, 'all');
        union = sum((rand_mm_mask + net_mask) > 0, 'all');
        perm_jac_mat_yeo(f, p) = intersection/union;
    end

end

% calculate p-vals from permutation test
jac_mat_yeo_p_vals = nan(7, 1);
for f = 1:length(yeo_7_nets)
    jac_ind = jac_mat_yeo(f);
    perm_jac_inds = sort(squeeze(perm_jac_mat_yeo(f, :)));
    jac_mat_yeo_p_vals(f) = sum(perm_jac_inds >= jac_ind, 'all')/num_perm;
end

% correct for multiple comparisons
jac_yeo_p_vals_vector = squeeze(reshape(jac_mat_yeo_p_vals, 1, []));
jac_yeo_fdrs = mafdr(jac_yeo_p_vals_vector, 'BHFDR', true);


%% Hierarchical clustering using clusterwise

output_dir = 'social_main\3_group_analysis\hierarchical_clustering';

msk_dir = 'social_main\3_group_analysis\VOIs\q0.01_thr1000_all_msks';
msk_pattern = '*.msk';
msk_names = dir(fullfile(msk_dir, msk_pattern));

clus_mean_betas = zeros([scales_num, length(msk_names)]);

for msk = 1:length(msk_names)
    current_msk = msk_names(msk).name;
    disp(current_msk);

    current_VOI_mask = xff(fullfile(msk_dir, current_msk)); current_VOI_mask = current_VOI_mask.Mask;
    current_VOI_mask = double(current_VOI_mask);
    current_VOI_mask_4d = permute(repmat(current_VOI_mask, [1, 1, 1, 5]), [4, 1, 2, 3]);

    group_beta_maps_masked = group_beta_maps .* current_VOI_mask_4d;
    betas = squeeze(mean(group_beta_maps_masked, [2, 3, 4]));
    clus_mean_betas(:, msk) = betas;

end

col_labels = {'Scale 1', 'Scale 2', 'Scale 3', 'Scale 4', 'Scale 5'};
row_labels = {'rMTG', 'lHC', 'lPHC', 'lTPJ', 'rTPJ', 'rPCS', 'Prc', 'rSMA', 'rRSC', 'dmPFC', 'vmPFC', 'lRSC'};
row_colors = {[.5 .5 .5], [.1 .1 .1], [.1 .1 .1], [.5 .5 .5], [.9 .9 .9], [.75 .9 .5], [.9 .9 .9], [.5 .5 .5], [.1 .1 .1], [.9 .9 .9], [.9 .9 .9], [.1 .1 .1]};
rowLabelColors = struct('Labels', row_labels, 'Colors', row_colors);

clustergram(clus_mean_betas', 'Colormap', redbluecmap, 'Standardize', 'Row', ...
    'RowLabels', row_labels, 'ColumnLabels', col_labels, ...
    'RowLabelsColor', rowLabelColors, 'LabelsWithMarkers', true)
