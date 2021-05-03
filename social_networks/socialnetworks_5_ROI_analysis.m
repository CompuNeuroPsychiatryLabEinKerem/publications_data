% RSA analysis in ROIs


%% Loading the ROIs

% Getting one of the statistical maps, to resample all of the other images to its resolution
filename_for_resampling = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_facebook_percentcommonfriends_radius3_iters10000.nii';
ROIs_directory = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\ROIs';

% Getting the RSC, PPA, OPA and FFA from the masks of Julian et al. 2012
ROIs = {}; ROI_names = {}; counter = 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\lRSC.img', filename_for_resampling, 0); ROI_names{counter} = 'lRSC'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\rRSC.img', filename_for_resampling, 0); ROI_names{counter} = 'rRSC'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\lPPA.img', filename_for_resampling, 0); ROI_names{counter} = 'lPPA'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\rPPA.img', filename_for_resampling, 0); ROI_names{counter} = 'rPPA'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\lTOS.img', filename_for_resampling, 0); ROI_names{counter} = 'lOPA'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Research\Julian_Kanwisher_2012_ROIs\scene_parcels\rTOS.img', filename_for_resampling, 0); ROI_names{counter} = 'rOPA'; counter = counter + 1;

% Getting the hippocampus ROI from the AAL atlas
AAL_atlas = reslice_data('C:\Users\michaelpeer1\Dropbox\Michael_scripts\Templates\AAL_61x73x61_YCG.nii', filename_for_resampling, 0);
ROIs{counter} = single(AAL_atlas == 37); ROI_names{counter} = 'lHC'; counter = counter + 1;
ROIs{counter} = single(AAL_atlas == 38); ROI_names{counter} = 'rHC'; counter = counter + 1;

% Orientation study ROI
ROIs{counter} = reslice_data('C:\Users\michaelpeer1\Dropbox\Michael_scripts\Templates\Orientation_paper_activations\social_vs_control_and_rest_ROI.nii', filename_for_resampling, 0); ROI_names{counter} = 'orientation_social'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Users\michaelpeer1\Dropbox\Michael_scripts\Templates\Orientation_paper_activations\space_vs_control_and_rest_ROI.nii', filename_for_resampling, 0); ROI_names{counter} = 'orientation_space'; counter = counter + 1;

% % Adding the ROI of the region identified in the searchlight - only for the pattern reliability analysis
% ROIs{counter} = reslice_data('C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_facebook_percentcommonfriends_radius3_othermats_regressed_percent_iters10000.nii', filename_for_resampling, 0); ROI_names{counter} = 'searchlight_facebookdist_regressed'; counter = counter + 1;
% ROIs{end} = ROIs{end}>=2.33;
% ROIs{counter} = reslice_data('C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_self_proximity_radius3_othermats_regressed_percent_iters10000.nii', filename_for_resampling, 0); ROI_names{counter} = 'searchlight_selfproximity_regressed'; counter = counter + 1;
% ROIs{end} = ROIs{end}>=2.33;
% ROIs{counter} = reslice_data('C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_personality_radius3_othermats_regressed_percent_iters10000.nii', filename_for_resampling, 0); ROI_names{counter} = 'searchlight_personality_regressed'; counter = counter + 1;
% ROIs{end} = ROIs{end}>=2.33;

% % Write all ROIs to files
% for i = 1:length(ROIs)
%     save_mat_to_nifti(filename_for_resampling, ROIs{i}, fullfile(ROIs_directory, [ROI_names{i} '.nii']));
% end



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
all_means_diag_nondiag_betweenruns = nan(length(ROIs), length(subject_names));
all_means_diag_nondiag_betweenruns_groups = nan(length(ROIs), length(subject_names));

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
        %         temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_RDM_shortestpath.mat'));
        %         temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_RDM_communicability.mat'));
        %         temp2 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_RDM_directconnections.mat'));
        %         temp3 = load(fullfile(current_subject_RDMs_dir, 'social_groups_RDM.mat'));
        temp3 = load(fullfile(current_subject_RDMs_dir, 'subject_responses_RDMs.mat'));
        RDMs = {};
        temp_counter = 1;
        RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_percentcommonfriends_RDM'; temp_counter = temp_counter + 1;
        %         RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_numcommonfriends_RDM'; temp_counter = temp_counter + 1;
        %         RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_RDM_shortestpath'; temp_counter = temp_counter + 1;
        %         RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_RDM_communicability'; temp_counter = temp_counter + 1;
        %         RDMs{temp_counter} = temp2.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_RDM_directconnections'; temp_counter = temp_counter + 1;
        %         RDMs{temp_counter} = temp3.social_groups_RDM;              RDM_names{temp_counter} = 'social_groups_RDM'; temp_counter = temp_counter + 1;
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
        
        % Measuring the consistency of patterns between runs
        num_runs = length(data_files);
        curr_corr = nan(num_conditions, num_conditions, num_runs, num_runs);
        means_diag_minus_nondiag = [];
        means_diag_minus_nondiag_groups = [];
        for run1 = 1:num_runs - 1
            for run2 = run1+1:num_runs
                for i = 1:num_conditions
                    for j = 1:num_conditions
                        curr_corr(i,j,run1,run2) = corr(ds.samples(i + num_conditions*(run1-1),:)',ds.samples(j + num_conditions*(run2-1),:)');
                    end
                end
                
                % Calculating the mean of the diagonal minus nondiagonal elements - consistency across runs
                ccc = curr_corr(:, :, run1, run2); 
                means_diag_minus_nondiag = [means_diag_minus_nondiag, mean(diag(ccc)) - mean(ccc(find(~eye(num_conditions))))];
                % Testing the effect of social groups - the mean of cells corresponding to correlation of individuals to themselves, minus the mean of cells corresponding to the correlation of individuals to their social group
                groups_mat = zeros(24); groups_mat(1:6,1:6)=1; groups_mat(7:12,7:12)=1; groups_mat(13:18,13:18)=1; groups_mat(19:24,19:24)=1; groups_mat = groups_mat & ~eye(24);
                means_diag_minus_nondiag_groups = [means_diag_minus_nondiag_groups, mean(diag(ccc)) - mean(ccc(find(~groups_mat)))];
            end
        end
        all_means_diag_nondiag_betweenruns(r,s) = nanmean(means_diag_minus_nondiag);
        all_means_diag_nondiag_betweenruns_groups(r,s) = nanmean(means_diag_minus_nondiag_groups);        
    end
end



% Test significance
r = 1;      % The ROI we are looking at
curr_data = all_results_reg{r}'; curr_data = curr_data(~isnan(curr_data(:,1)),:);
% Calculate the overall repeated-measures anova
p_anova = anova_rm(curr_data);
% Conduct post-hoc tests - dependent-samples t-tests, corrected for multiple comparisons by Bonferroni-Hochberg FDR
num_columns = size(curr_data,2);
ps=nan(4); stats_all=nan(4);
for i=1:num_columns,for j=i+1:num_columns, [~,ps(i,j),~,STATS]=ttest(curr_data(:,i),curr_data(:,j)); stats_all(i,j)=STATS.tstat; end,end
ps=ps(~isnan(ps)); stats_all=stats_all(~isnan(stats_all));
corrected_ps_fdr = mafdr(ps,'BHFDR','true');
[corrected_ps_bonfholm,~] = bonf_holm(ps,0.05);

% Plotting the results
% Selecting the data
roi = 2;            % The ROI of interest
a=all_results_reg{roi}';
% Creating a bar graph
figure; hold on;
% bar(nanmean(a), 0.4, 'FaceColor',[.5 .5 .5]);
bar(nanmean(a), 0.4,'FaceColor',[68,114,196]/256);
% Removing the Y axis tickmarks
h = gca; h.YAxis.TickLength = [0,0];
% Creating a scatter plot
% colors = [repmat([1,0,0],18,1); repmat([0,0,1],18,1); repmat([0,1,0],18,1); repmat([1,1,0],18,1)];
% scatter([ones(1,18),ones(1,18)*2,ones(1,18)*3,ones(1,18)*4],a(:),200,colors,'.','jitter','on','jitterAmount',0.1,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% colors = repmat([0.3,0.3,0.3],18*4,1);
colors = repmat([0.25,0.25,0.75],18*4,1);
scatter([ones(1,18),ones(1,18)*2,ones(1,18)*3,ones(1,18)*4],a(:),150,colors,'.','jitter','on','jitterAmount',0.1)
% Adding error bars
h = errorbar(nanmean(a),nanstd(a)/sqrt(17),'+k');
h.CapSize = 12;  h.LineWidth = 1;   % Changing the bars' length and width
% Changing the axes' appearance
xlim([0.3,4.7]), ylim([-0.15,0.4])  % Changing the X and Y axis length
line(xlim, [0 0],'Color',[.5,.5,.5]);   % Adding an X axis line
set(gca,'box','off','xcolor','w','ycolor',[.5,.5,.5])   % Changing the axes' colors


% Test significance with FDR correction across ROIs
rdm = 1;
relevant_ROIs = 1:10;
curr_data = [];
for r=1:length(relevant_ROIs)
    curr_data(:,r) = all_results_reg{r}(rdm,:)';
end
[~,all_ps,~,STATS] = ttest(curr_data, []);
fdr_corrected_ps = mafdr(all_ps,'BHFDR','true');


% Test significance of reliability of patterns
relevant_ROIs = [1, 2, 9, 10];
[~,p1,~,STATS1] = ttest(all_means_diag_nondiag_betweenruns(relevant_ROIs, :)',[],'tail','right');
p1_corrected = mafdr(p1,'BHFDR','true');
[~,p2,~,STATS2] = ttest(all_means_diag_nondiag_betweenruns_groups(relevant_ROIs, :)',[],'tail','right');
p2_corrected = mafdr(p2,'BHFDR','true');



%% Check correlation in ROIs found in the searchlight analysis (original correlation values in regions identified by regressing out the variables)
num_iterations = 10000;     % Number of iterations for the permutation testing
searchlight_radius = 3;
RSA_filenames = {'RSA_searchlight_facebook_distance_percentcommonfriends_RDM_radius3.nii',...       % Social distance by proportion of common friends
    'RSA_searchlight_responses_self_proximity_RDM_radius3.nii',...                                  % Dissimilarity in responses to personal affiliation questions
    'RSA_searchlight_responses_personality_RDM_radius3.nii',...                                     % Dissimilarity in responses to personality traits questions
    'RSA_searchlight_responses_appearance_RDM_radius3.nii',...                                      % Dissimilarity in responses to appearance questions
    };
% Defining output filenames for each analysis
RSA_output_filenames = {['group_analysis_RSA_searchlight_facebook_percentcommonfriends_radius' num2str(searchlight_radius) '_othermats_regressed_percent_iters' num2str(num_iterations) '.nii'],...
    ['group_analysis_RSA_searchlight_self_proximity_radius' num2str(searchlight_radius) '_othermats_regressed_percent_iters' num2str(num_iterations) '.nii'],...
    ['group_analysis_RSA_searchlight_personality_radius' num2str(searchlight_radius) '_othermats_regressed_percent_iters' num2str(num_iterations) '.nii'],...
    ['group_analysis_RSA_searchlight_appearance_radius' num2str(searchlight_radius) '_othermats_regressed_percent_iters' num2str(num_iterations) '.nii'],...
    };
stat_threshold = 2.33;
for o = 1:length(RSA_output_filenames)
    RSA_output_filename = RSA_output_filenames{o};    
    current_RSA_mask = spm_read_vols(spm_vol(fullfile(all_subjects_analysis_results_dir, RSA_output_filename)));
    current_RSA_mask = current_RSA_mask >= stat_threshold;
    
    for i = 1:length(RSA_filenames)
        RSA_filename = RSA_filenames{i};
        
        % Reading the RSA analysis results files from all subjects
        all_RSA_results = [];
        for s = 1:length(subject_names)
            subj = subject_names(s).name;
            current_subject_results_dir = fullfile(all_subjects_analysis_results_dir, subj);
            RSA_result_file = dir(fullfile(current_subject_results_dir, RSA_filename));
            % Converting to CosmoMVPA dataset format
            if ~isempty(RSA_result_file)
                curr_subj_RSA_result = spm_read_vols(spm_vol(fullfile(current_subject_results_dir, RSA_result_file(1).name)));
                all_regions_subjects_RSA_corr(s,i,o) = nanmean(curr_subj_RSA_result(current_RSA_mask ~= 0));
            end
        end
    end
end

r = 1;      % The ROI we are looking at
subjects_to_remove = 15;
curr_data = all_regions_subjects_RSA_corr(:,:,r); curr_data(subjects_to_remove,:) = [];
% Calculate the overall repeated-measures anova
p_anova = anova_rm(curr_data);
% Conduct post-hoc tests - dependent-samples t-tests, corrected for multiple comparisons by Bonferroni-Hochberg FDR
num_columns = size(curr_data,2);
ps=nan(4); 
for i=1:num_columns,for j=i+1:num_columns, [~,ps(i,j)]=ttest(curr_data(:,i),curr_data(:,j)); end,end
ps=ps(~isnan(ps));
corrected_ps_fdr = mafdr(ps,'BHFDR','true');
[corrected_ps_bonfholm,~] = bonf_holm(ps,0.05);


% Plotting the results
% Reading the data
r = 1;            % The ROI of interest          
a = all_regions_subjects_RSA_corr(:,:,r);
% Creating a bar graph
figure; hold on;
bar(nanmean(a), 0.4, 'FaceColor', [68,114,196]/256);
% Removing the Y axis tickmarks
h = gca; h.YAxis.TickLength = [0,0];
% Creating a scatter plot
colors = repmat([0.25,0.25,0.75],18*4,1);
scatter([ones(1,18),ones(1,18)*2,ones(1,18)*3,ones(1,18)*4],a(:),150,colors,'.','jitter','on','jitterAmount',0.1)
% Adding error bars
h = errorbar(nanmean(a),nanstd(a)/sqrt(17),'+k');
h.CapSize = 12;  h.LineWidth = 1;   % Changing the bars' length and width
% Changing the axes' appearance
xlim([0.3,4.7]), ylim([-0.15,0.4])  % Changing the X and Y axis length
line(xlim, [0 0],'Color',[.5,.5,.5]);   % Adding an X axis line
set(gca,'box','off','xcolor','w','ycolor',[.5,.5,.5])   % Changing the axes' colors 
% Adding lines for partial correlation coefficients, on the same graph
b = all_results_reg{10 + r}';
h = errorbar(nanmean(b),[0,0,0,0],'.r');
h.CapSize = 20;  h.LineWidth = 1.5;


%% Find clusters and their peak coordinates, sizes and anatomical labels
r = {};
r{1} = spm_read_vols(spm_vol('C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_facebook_percentcommonfriends_radius3_othermats_regressed_percent_iters10000.nii'));
r{2} = spm_read_vols(spm_vol('C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_self_proximity_radius3_othermats_regressed_percent_iters10000.nii'));
r{3} = spm_read_vols(spm_vol('C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_personality_radius3_othermats_regressed_percent_iters10000.nii'));
aicha = reslice_data('C:\Users\michaelpeer1\Dropbox (Personal)\Michael_scripts\Templates\AICHA.nii','C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_personality_radius3_othermats_regressed_percent_iters10000.nii',0);

roi = 1;
curr_region = r{roi};
curr_region(curr_region<2.33) = 0;
% Finding the voxels of interest
v = find(curr_region>=2.33);
[coordsx,coordsy,coordsz] = ind2sub(size(curr_region), v);
% Separating to clusters using the spm_clusters function
A = spm_clusters([coordsx,coordsy,coordsz]');
num_clusters = length(unique(A));
sizes_clusters = []; peak_values_clusters = []; cluster_label = [];
for i = 1:num_clusters
    sizes_clusters(i) = length(find(A==i)); 
    peak_values_clusters(i) = max(curr_region(v(A==i)));
    curr_voxels_aicha_values = aicha(v(find(A==i)));
    cluster_label(i) = mode(curr_voxels_aicha_values(curr_voxels_aicha_values~=0));
end
% TO FIND CLUSTER PEAK MNI COORDINATES - LOOK MANUALLY FOR THE VOXELS USING MRICRON
