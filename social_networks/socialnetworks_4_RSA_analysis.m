%% Analysis of the functional data using RSA

% Defining parent directories (each containing a separate directory for each subject)
all_subjects_GLM_files_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\GLM_files';          % The per-run GLM results
all_subjects_RDMs_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\Dissimilarity_mats';      % The dissimilarity matrices
all_subjects_analysis_results_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results';          % Directory for analysis output files
% Identifying the subjects' names (directories)
subject_names = dir(all_subjects_GLM_files_dir); subject_names = subject_names(3:end); subject_names = subject_names([subject_names.isdir]);

% Defining the number of stimuli (individuals)
groups_num = 4; people_in_group_num = 6;
num_conditions = groups_num * people_in_group_num;

% Defining the searchlight radius (voxels)
searchlight_radius = 3;


%% RSA analysis in each subject
all_RSA_results = cell(1, 6); for i = 1:6, all_RSA_results{i} = cell(1, num_conditions); end
for s = 1:length(subject_names)
    %% Reading the data and converting it to CosmoMVPA dataset format
    
    % Getting the current subject's directories
    subj = subject_names(s).name;
    disp(s)
    current_subject_GLM_dir = fullfile(all_subjects_GLM_files_dir, subj);
    current_subject_RDMs_dir = fullfile(all_subjects_RDMs_dir, subj);
    current_subject_results_dir = fullfile(all_subjects_analysis_results_dir, subj);
    
    % Reading the data for analysis
    data_files = dir(fullfile(current_subject_GLM_dir, '*t_map_individuals*.vmp'));       % GLM files - t values
    
    % Create a cosmoMVPA dataset from the whole brain GLM data
    all_ds = cell(1, length(data_files));
    for i = 1:length(data_files)
        % Reading the current GLM file to CosmoMVPA dataset format
        ds_current = cosmo_fmri_dataset(fullfile(current_subject_GLM_dir, data_files(i).name));
        % Getting only the first 24 conditions - the rest are not interesting (motion parameters, etc.)
        ds_current = cosmo_slice(ds_current, 1:num_conditions, 1);
        % Setting the targets (condition numbers) to be 1-24
        ds_current.sa.targets = [1:num_conditions]';
        % Setting the chunk (independent data run) to be the run number
        ds_current.sa.chunks = ones(num_conditions, 1) * i;
        % Adding the run to the overall subject dataset
        all_ds{i} = ds_current;
    end
    % Combine all runs to a single dataset
    ds = cosmo_stack(all_ds);
    % Remove voxels with no information from the data
    ds = cosmo_remove_useless_data(ds);
    
    % Averaging the pattern for each individual across runs
    num_conditions = length(unique(ds.sa.targets));
    ds_average_runs = ds;
    for i = 1:num_conditions
        ds_average_runs.samples(i,:) = mean(ds.samples(i:num_conditions:size(ds.samples,1),:));
    end
    ds_average_runs = cosmo_slice(ds_average_runs, 1:num_conditions);
    
    
    
    %% Running the searchlight analysis using the subject's dissimilarity matrices
    
    % Reading the different dissimilarity matrices
    RDMs = {}; RDM_names = {};
    temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_percentcommonfriends_RDM.mat'));     % Social distance by proportion of mutual friends
    %     temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_numcommonfriends_RDM.mat'));   % Alternative social distance measure - total number of mutual friends
    %     temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_RDM_directconnections.mat'));  % Alternative social distance measure - existence of a direct connection between individuals
    %     temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_RDM_shortestpath.mat'));       % Alternative social distance measure - shortest path between individuals along the social network
    %     temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_RDM_communicability.mat'));    % Alternative social distance measure - graph theoretical communicability
    temp3 = load(fullfile(current_subject_RDMs_dir, 'subject_responses_RDMs.mat'));
    
    % Loading the actual matrices
    temp_counter = 1;
    RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_percentcommonfriends_RDM'; temp_counter = temp_counter + 1;
    %     RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_numcommonfriends_RDM'; temp_counter = temp_counter + 1;
    %     RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_directconnections_RDM'; temp_counter = temp_counter + 1;
    %     RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_shortestpath_RDM'; temp_counter = temp_counter + 1;
    %     RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_communicability_RDM'; temp_counter = temp_counter + 1;
    RDMs{temp_counter} = temp3.responses_self_proximity_RDM;   RDM_names{temp_counter} = 'responses_self_proximity_RDM'; temp_counter = temp_counter + 1;
    RDMs{temp_counter} = temp3.responses_personality_RDM;      RDM_names{temp_counter} = 'responses_personality_RDM'; temp_counter = temp_counter + 1;
    RDMs{temp_counter} = temp3.responses_appearance_RDM;       RDM_names{temp_counter} = 'responses_appearance_RDM'; temp_counter = temp_counter + 1;
    
    % Removing places that are NaN values in one of the matrices from all of the others
    sum_RDMs = zeros(size(RDMs{1})); for i = 1:length(RDMs), sum_RDMs = sum_RDMs + RDMs{i}; end   % Creating a matrix that is a sum of the others, to identify locations with nan values in one of the matrices
    for i = 1:length(RDMs), RDMs{i}(isnan(sum_RDMs)) = nan; end
    
    % Defining variables for the searchlight
    current_ds = ds_average_runs;       % Current cosmoMVPA dataset to run the analysis on
    nbrhood = cosmo_spherical_neighborhood(current_ds, 'radius', searchlight_radius);    % Searchlight radius
    measure_RSA = @cosmo_target_dsm_corr_measure;               % Defining the computation at each searchlight sphere - in this case, correlation to the dissimilarity matrix
    measure_args = struct();
    measure_args.center_data = true;    % Centering the data by removing mean
    measure_args.type = 'Spearman';     % Use the rank-based Spearman's correlation to compare dissimilarity matrices to one another (Pearson's correlation is used to calculate the distances themselves between neural patterns)
    
    % Running the searchlight separately for each of the dissimilarity matrices, with or without regression of the matrices from each other
    for r = 1:length(RDMs)
        % Defining the current dissimilarity matrix
        current_RDM = RDMs{r};
        if nansum(current_RDM(:)) ~= 0     % Ignoring empty matrices - subjects for which no data was collected, or another problem
            measure_args.target_dsm = current_RDM;
            % Running the searchlight on the dissimilarity matrix
            RSA_results = cosmo_searchlight(current_ds, nbrhood, measure_RSA, measure_args);
            % Saving the searchlight results to nifti file
            output_filename = fullfile(current_subject_results_dir, ['RSA_searchlight_' RDM_names{r} '_radius' num2str(searchlight_radius) '.nii']);
            cosmo_map2fmri(RSA_results, output_filename);
            
            % Running the searchlight again, this time regressing out the other matrices to eliminate their contribution
            current_RDMs_to_regress_out = {};
            for i = 1:length(RDMs)
                if i ~= r
                    current_RDMs_to_regress_out{end+1} = RDMs{i};
                end
            end
            measure_args_reg = measure_args;
            measure_args_reg.regress_dsm = current_RDMs_to_regress_out;
            % Running the searchlight
            RSA_results_regressed = cosmo_searchlight(current_ds, nbrhood, measure_RSA, measure_args_reg);
            % Saving the searchlight results to nifti file
            output_filename = fullfile(current_subject_results_dir, ['RSA_searchlight_' RDM_names{r} '_radius' num2str(searchlight_radius) '_othermats_regressed_percent.nii']);
            cosmo_map2fmri(RSA_results_regressed, output_filename);
        end
    end
end



%% Group analysis from individual subjects' results

% MANUALLY REMOVE THE RESULTS WITH REGRESSION FOR ONE SUBJECT FOR WHICH RESPONSES WERE NOT RECORDED (SUBJECT 15)

% Define analysis parameters
num_iterations = 10000;     % Number of iterations for the permutation testing

% Define RSA files to run the group analysis on (separately for each RSA file, across all subjects)
RSA_filenames = {'RSA_searchlight_facebook_distance_percentcommonfriends_RDM_radius3.nii',...       % Social distance by proportion of common friends
    'RSA_searchlight_responses_self_proximity_RDM_radius3.nii',...                                  % Dissimilarity in responses to personal affiliation questions
    'RSA_searchlight_responses_personality_RDM_radius3.nii',...                                     % Dissimilarity in responses to personality traits questions
    'RSA_searchlight_responses_appearance_RDM_radius3.nii',...                                      % Dissimilarity in responses to appearance questions
    'RSA_searchlight_facebook_distance_percentcommonfriends_RDM_radius3_othermats_regressed_percent.nii',...    % Social distance by proportion of common friends, other matrices regressed out
    'RSA_searchlight_responses_self_proximity_RDM_radius3_othermats_regressed_percent.nii',...                  % Dissimilarity in responses to personal affiliation questions, other matrices regressed out
    'RSA_searchlight_responses_personality_RDM_radius3_othermats_regressed_percent.nii',...                     % Dissimilarity in responses to personality traits questions, other matrices regressed out
    'RSA_searchlight_responses_appearance_RDM_radius3_othermats_regressed_percent.nii',...                      % Dissimilarity in responses to appearance questions, other matrices regressed out
    };
%     'RSA_searchlight_facebook_distance_communicability_RDM_radius3.nii',...                       % Social distance by communicability
%     'RSA_searchlight_facebook_distance_directconnections_RDM_radius3.nii',...                     % Social distance by direct connectivity
%     'RSA_searchlight_facebook_distance_numcommonfriends_RDM_radius3.nii',...                      % Social distance by total number of mutual friends
%     'RSA_searchlight_facebook_distance_shortestpath_RDM_radius3.nii',...                          % Social distance by shortest path length


% Defining output filenames for each analysis
RSA_output_filenames = {['group_analysis_RSA_searchlight_facebook_percentcommonfriends_radius' num2str(searchlight_radius) '_iters' num2str(num_iterations) '.nii'],...
    ['group_analysis_RSA_searchlight_self_proximity_radius' num2str(searchlight_radius) '_iters' num2str(num_iterations) '.nii'],...
    ['group_analysis_RSA_searchlight_personality_radius' num2str(searchlight_radius) '_iters' num2str(num_iterations) '.nii'],...
    ['group_analysis_RSA_searchlight_appearance_radius' num2str(searchlight_radius) '_iters' num2str(num_iterations) '.nii'],...
    ['group_analysis_RSA_searchlight_facebook_percentcommonfriends_radius' num2str(searchlight_radius) '_othermats_regressed_percent_iters' num2str(num_iterations) '.nii'],...
    ['group_analysis_RSA_searchlight_self_proximity_radius' num2str(searchlight_radius) '_othermats_regressed_percent_iters' num2str(num_iterations) '.nii'],...
    ['group_analysis_RSA_searchlight_personality_radius' num2str(searchlight_radius) '_othermats_regressed_percent_iters' num2str(num_iterations) '.nii'],...
    ['group_analysis_RSA_searchlight_appearance_radius' num2str(searchlight_radius) '_othermats_regressed_percent_iters' num2str(num_iterations) '.nii'],...
    };
%     ['group_analysis_RSA_searchlight_facebook_distance_communicability_RDM_radius' num2str(searchlight_radius) '_iters' num2str(num_iterations) '.nii'],...
%     ['group_analysis_RSA_searchlight_facebook_distance_directconnections_RDM_radius' num2str(searchlight_radius) '_iters' num2str(num_iterations) '.nii'],...
%     ['group_analysis_RSA_searchlight_facebook_distance_numcommonfriends_RDM_radius' num2str(searchlight_radius) '_iters' num2str(num_iterations) '.nii'],...
%     ['group_analysis_RSA_searchlight_facebook_distance_shortestpath_RDM_radius' num2str(searchlight_radius) '_iters' num2str(num_iterations) '.nii'],...


% Running the group analysis separately for each dissimilarity matrix results
for i = 1:length(RSA_filenames)
    RSA_filename = RSA_filenames{i};
    RSA_output_filename = RSA_output_filenames{i};

    % Reading the RSA analysis results files from all subjects
    all_RSA_results = cell(length(subject_names), 1);
    for s = 1:length(subject_names)
        subj = subject_names(s).name;
        current_subject_results_dir = fullfile(all_subjects_analysis_results_dir, subj);
        RSA_result_file = dir(fullfile(current_subject_results_dir, RSA_filename));
        % Converting to CosmoMVPA dataset format
        if ~isempty(RSA_result_file)
            all_RSA_results{s} = cosmo_fmri_dataset(fullfile(current_subject_results_dir, RSA_result_file(1).name));
        end
    end
    
    % Organizing the across-subjects dataset
    all_RSA_results(cellfun(@isempty, all_RSA_results)) = [];
    [idxs, ds_cell_common] = cosmo_mask_dim_intersect(all_RSA_results);     % Choosing only voxels overlapping across all subjects
    ds_group = cosmo_stack(ds_cell_common, 1, 'drop_nonunique');
    ds_group.sa.targets = ones(size(ds_group.samples, 1),1);
    ds_group.sa.chunks = [1:size(ds_group.samples, 1)]';
    % Removing uninformative locations (NaN values, etc.)
    ds_group = cosmo_remove_useless_data(ds_group);
    % Doing fisher-z transform on the correlation values (using the atanh function)
    ds_group.samples = atanh(ds_group.samples);
    % Running the group analysis
    h0_mean = 0;    % The null hypothesis - no correlation
    cluster_nbrhood = cosmo_cluster_neighborhood(ds_group);
    stat_map = cosmo_montecarlo_cluster_stat(ds_group, cluster_nbrhood, 'niter', num_iterations, 'h0_mean', h0_mean);
    % saving the result as a nifti file
    cosmo_map2fmri(stat_map, fullfile(all_subjects_analysis_results_dir, RSA_output_filename));
end




%% Computing overlap of the group analysis results with the 7 resting-state networks from Yeo et al. 2011
% Loading the Yeo et al. networks
ROIs = {}; ROI_names = {}; counter = 1;
Yeo_nets = reslice_data('C:\Users\michaelpeer1\Dropbox\Michael_templates\rYeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii', filename_for_resampling, 0); 
ROIs{counter} = Yeo_nets; ROIs{counter}(find(ROIs{counter}~=1)) = 0; ROI_names{counter} = 'Yeo_visual'; counter = counter + 1;
ROIs{counter} = Yeo_nets; ROIs{counter}(find(ROIs{counter}~=2)) = 0; ROI_names{counter} = 'Yeo_somatomotor'; counter = counter + 1;
ROIs{counter} = Yeo_nets; ROIs{counter}(find(ROIs{counter}~=3)) = 0; ROI_names{counter} = 'Yeo_DAN'; counter = counter + 1;
ROIs{counter} = Yeo_nets; ROIs{counter}(find(ROIs{counter}~=4)) = 0; ROI_names{counter} = 'Yeo_VAN'; counter = counter + 1;
ROIs{counter} = Yeo_nets; ROIs{counter}(find(ROIs{counter}~=5)) = 0; ROI_names{counter} = 'Yeo_OFC_ATL'; counter = counter + 1;
ROIs{counter} = Yeo_nets; ROIs{counter}(find(ROIs{counter}~=6)) = 0; ROI_names{counter} = 'Yeo_FPCN'; counter = counter + 1;
ROIs{counter} = Yeo_nets; ROIs{counter}(find(ROIs{counter}~=7)) = 0; ROI_names{counter} = 'Yeo_DMN'; counter = counter + 1;

% Loading the group analysis result for the main analysis
group_activation = niftiread('C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_facebook_percentcommonfriends_radius3_iters10000.nii');

% Computing overlap of result with each resting-state network
ROIs_overlap2 = zeros(1, length(ROIs));
for i = 1:length(ROIs)
    % Computing overlap as number of ROI voxels that appear in the RSA results above the z=2.33 (p<0.01) threshold
    num_voxels_in_current_ROI = length(find(ROIs{i} ~= 0));
    num_voxels_significant_in_current_ROI = length(find(ROIs{i} ~= 0 & group_activation>=2.33));
    num_voxels_significant_overall = length(find(group_activation>=2.33 & Yeo_nets~=0));      % For comparisons to RSNs
    ROIs_overlap2(i) = num_voxels_significant_in_current_ROI / num_voxels_significant_overall;
end
disp(' '); for i = 1:length(ROIs), disp([num2str(ROIs_overlap2(i)) ' ' ROI_names{i}]); end

% Permutation test to check significance of overlap compared to overlap obtained by chance
num_perms = 1000;
num_Yeo_nets = 7;
ROIs_overlaps2_perms = nan(num_perms, num_Yeo_nets);
locs_GM = find(Yeo_nets(:) ~= 0);
Yeo_nets_GM = Yeo_nets(locs_GM);
for p = 1:num_perms
    Yeo_nets_current = Yeo_nets; 
    Yeo_nets_current(locs_GM) = Yeo_nets_GM(randperm(length(Yeo_nets_GM)));
    for i = 1:num_Yeo_nets
        num_voxels_significant_in_current_ROI = length(find(Yeo_nets_current == i & group_activation>=2.33));
        num_voxels_significant_overall = length(find(group_activation>=2.33 & Yeo_nets~=0));      % For comparisons to RSNs
        ROIs_overlaps2_perms(p, i) = num_voxels_significant_in_current_ROI / num_voxels_significant_overall;
    end
end
pvalues_Yeo_nets_overlap = nan(1, num_Yeo_nets);
for i = 1:num_Yeo_nets
    pvalues_Yeo_nets_overlap(i) = 1 - sum(ROIs_overlap2(i+10) >= ROIs_overlaps2_perms(:,i)) / num_perms;
end




