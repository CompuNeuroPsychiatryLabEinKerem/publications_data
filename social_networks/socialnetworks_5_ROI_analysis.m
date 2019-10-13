% RSA analysis in ROIs


%% Loading the ROIs

% Getting one of the statistical maps, to resample all of the other images to its resolution
filename_for_resampling = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_facebook_percentcommonfriends_radius3_iters10000.nii';

% Getting the RSC, PPA, OPA and FFA from the masks of Julian et al. 2012
ROIs = {}; ROI_names = {}; counter = 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\lRSC.img', filename_for_resampling, 0); ROI_names{counter} = 'lRSC'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\rRSC.img', filename_for_resampling, 0); ROI_names{counter} = 'rRSC'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\lPPA.img', filename_for_resampling, 0); ROI_names{counter} = 'lPPA'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\rPPA.img', filename_for_resampling, 0); ROI_names{counter} = 'rPPA'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\lTOS.img', filename_for_resampling, 0); ROI_names{counter} = 'lOPA'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\rTOS.img', filename_for_resampling, 0); ROI_names{counter} = 'rOPA'; counter = counter + 1;

% Getting the hippocampus ROI from the AAL atlas
AAL_atlas = reslice_data('C:\Users\michaelpeer1\Dropbox\Michael_templates\AAL_61x73x61_YCG.nii', filename_for_resampling, 0);
ROIs{counter} = single(AAL_atlas == 37); ROI_names{counter} = 'lHC'; counter = counter + 1;
ROIs{counter} = single(AAL_atlas == 38); ROI_names{counter} = 'rHC'; counter = counter + 1;

% Orientation study ROI
ROIs{counter} = reslice_data('C:\temp\Orientation_paper_activations\social_vs_control_and_rest_ROI.nii', filename_for_resampling, 0); ROI_names{counter} = 'orientation_social'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\temp\Orientation_paper_activations\space_vs_control_and_rest_ROI.nii', filename_for_resampling, 0); ROI_names{counter} = 'orientation_space'; counter = counter + 1;

% Write all ROIs to files
ROIs_directory = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\ROIs';
for i = 1:length(ROIs)
    save_mat_to_nifti(filename_for_resampling, ROIs{i}, fullfile(ROIs_directory, [ROI_names{i} '.nii']));
end


%% Loading the experimental data
all_subjects_GLM_files_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\GLM_files';          % The per-run GLM results
all_subjects_RDMs_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\Dissimilarity_mats';      % The dissimilarity matrices
all_subjects_analysis_results_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results';          % Directory for outlut
% Finding the subject names (directories)
subject_names = dir(all_subjects_GLM_files_dir); subject_names = subject_names(3:end); subject_names = subject_names([subject_names.isdir]);
% Defining the number of stimuli (individuals)
groups_num = 4; people_in_group_num = 6;
num_conditions = groups_num * people_in_group_num;



%% Computing RSA within ROI
all_results = cell(1, length(ROIs)); for i = 1:length(ROIs), all_results{i} = nan(4, length(subject_names)); end
all_results_reg = cell(1, length(ROIs)); for i = 1:length(ROIs), all_results_reg{i} = nan(4, length(subject_names)); end

for r = 1:length(ROIs)
    disp(ROI_names{r});
    
    % Getting the ROI mask filename
    current_ROI_mask_filename = fullfile(ROIs_directory, [ROI_names{r} '.nii']);
    
    % Reading the data from each subject and converting to cosmoMVPA dataset format
    for s = 1:length(subject_names)
        subj = subject_names(s).name;
        disp(s)
        
        % Defining subject-specific directories
        current_subject_GLM_dir = fullfile(all_subjects_GLM_files_dir, subj);
        current_subject_RDMs_dir = fullfile(all_subjects_RDMs_dir, subj);
        current_subject_results_dir = fullfile(all_subjects_analysis_results_dir, subj);
        
        % Reading the data for analysis
        data_files = dir(fullfile(current_subject_GLM_dir, '*t_map_individuals*.vmp'));       % GLM files - t values
        
        % Reading GLM data from the specific ROI and creating dataset
        all_ds = cell(1, length(data_files));
        for i = 1:length(data_files)
            ds_current = cosmo_fmri_dataset(fullfile(current_subject_GLM_dir, data_files(i).name), 'mask', current_ROI_mask_filename);  % Reading the current GLM file to CosmoMVPA dataset format
            ds_current = cosmo_slice(ds_current, 1:num_conditions, 1);          % Getting only the first 24 conditions - the rest are not interesting (motion parameters, etc.)
            ds_current.sa.targets = [1:num_conditions]';                        % Setting the targets (condition numbers) to be 1-24
            ds_current.sa.chunks = ones(num_conditions, 1) * i;                 % Setting the chunks (independent data runs) to be the run number
            all_ds{i} = ds_current;
        end
        ds = cosmo_stack(all_ds);                   % Combine all runs to a single dataset
        ds = cosmo_remove_useless_data(ds);         % Remove voxels with no information from the data
        
        % Averaging patterns across runs to get one pattern for each individual
        num_conditions = length(unique(ds.sa.targets));
        ds_average_runs = ds;
        for i = 1:num_conditions
            ds_average_runs.samples(i,:) = mean(ds.samples(i:num_conditions:size(ds.samples,1),:));
        end
        ds_average_runs = cosmo_slice(ds_average_runs, 1:num_conditions);
        
        % Loading the dissimilarity matrices
        temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_percentcommonfriends_RDM.mat'));
        %         temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_numcommonfriends_RDM.mat'));
        %         temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_RDM_directconnections.mat'));
        %         temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_RDM_shortestpath.mat'));
        %         temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_RDM_communicability.mat'));
        temp3 = load(fullfile(current_subject_RDMs_dir, 'subject_responses_RDMs.mat'));
        temp_counter = 1;
        RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_percentcommonfriends_RDM'; temp_counter = temp_counter + 1;
        %         RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_numcommonfriends_RDM'; temp_counter = temp_counter + 1;
        %         RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_RDM_directconnections'; temp_counter = temp_counter + 1;
        %         RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_RDM_shortestpath'; temp_counter = temp_counter + 1;
        %         RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_RDM_communicability'; temp_counter = temp_counter + 1;
        %         RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'social_groups_RDM'; temp_counter = temp_counter + 1;
        RDMs{temp_counter} = temp3.responses_self_proximity_RDM;   RDM_names{temp_counter} = 'responses_self_proximity_RDM'; temp_counter = temp_counter + 1;
        RDMs{temp_counter} = temp3.responses_personality_RDM;      RDM_names{temp_counter} = 'responses_personality_RDM'; temp_counter = temp_counter + 1;
        RDMs{temp_counter} = temp3.responses_appearance_RDM;       RDM_names{temp_counter} = 'responses_appearance_RDM'; temp_counter = temp_counter + 1;
        % Removing places that are NaN values in one of the matrices
        sum_RDMs = zeros(size(RDMs{1})); for i = 1:length(RDMs), sum_RDMs = sum_RDMs + RDMs{i}; end   % Creating a matrix that is a sum of the others, to identify locations with nan values in one of the matrices
        for i = 1:length(RDMs), RDMs{i}(isnan(sum_RDMs)) = nan; end
        
        % Running the RSA for each dissimilarity matrix, with or without regression of other matrices from it
        if sum(RDMs{2}(:))~=0       % Making sure there's no matrix that is all zeros - subject's responses were not recorded
            for rdm = 1:length(RDMs)
                % Defining the current dissimilarity matrix
                current_RDM = RDMs{rdm};
                if nansum(current_RDM(:)) ~= 0     % Ignoring empty matrices - subjects for which no data was collected, or another problem
                    % Defining searchlight parameters
                    measure_args = struct(); measure_args.center_data = true;
                    measure_args.type = 'Spearman';     % Use Spearman's correlation to compare dissimilarity matrices
                    measure_args.target_dsm = current_RDM;
                    % Running the searchlight on the dissimilarity matrix
                    result = cosmo_target_dsm_corr_measure(ds_average_runs, measure_args);
                    all_results{r}(rdm, s) = result.samples;
                    % Running the searchlight with regression of other matrices from the current matrix
                    current_RDMs_to_regress_out = {};
                    for i = 1:length(RDMs)
                        if i ~= rdm
                            current_RDMs_to_regress_out{end+1} = RDMs{i};
                        end
                    end
                    measure_args_reg = measure_args;
                    measure_args_reg.regress_dsm = current_RDMs_to_regress_out;
                    result_reg = cosmo_target_dsm_corr_measure(ds_average_runs, measure_args_reg);
                    all_results_reg{r}(rdm, s) = result_reg.samples;
                end
            end
        end
    end
end







