% Analysis of individual questions instead of whole runs

%% Starting brainvoyager from Matlab
bvqx = actxserver('BrainVoyager.BrainVoyagerScriptAccess.1');


%% Defining working directories and general parameters
subjects_analysis_dir = 'E:\Social_networks_analysis\Analysis_2mm_smoothing\';          % A directory where the preprocessed data will be put
subject_names = dir(fullfile(subjects_analysis_dir,'18*')); %subject_names = subject_names(3:end);       % Getting the names of the subjects directories

num_individuals = 24;

%% Creating new PRT files with separation to questions
for s = 1:length(subject_names)
    subj = subject_names(s).name;
    disp(subj);
    current_subj_analysis_dir = fullfile(subjects_analysis_dir, subj);
    
    prt_files = dir(fullfile(current_subj_analysis_dir, '*individuals*.prt'));
    
    for i = 1:length(prt_files)
        current_prt = xff(fullfile(prt_files(i).folder, prt_files(i).name));
        current_prt.NrOfConditions = current_prt.NrOfConditions + num_individuals;
        for cond = 1:num_individuals
            % Add a new condition with the timing of question 2
            current_prt.AddCond(['q2_person' num2str(cond)], current_prt.Cond(cond+1).OnOffsets(2,:));
            % Remove the timing of question 2 from the original condition
            current_prt.Cond(cond+1).OnOffsets = current_prt.Cond(cond+1).OnOffsets(1,:);
        end
        % Saving the new PRT file
        current_prt.SaveAs(fullfile(current_subj_analysis_dir, [subj '_questions' num2str(i) '.prt']));
    end
end


%% Creating design matrices (SDM files) for each subject, from the PRT files
disp('Creating single study design matrices (.sdm files)')
for s = 1:length(subject_names)
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);
    
    % Loading the normalized anatomical (VMR) file
    vmr_file = dir(fullfile(output_dir, '*MNI.vmr'));
    vmr_file = fullfile(output_dir, vmr_file(1).name);
    vmr = bvqx.OpenDocument(vmr_file);
    
    % Loading the normalized functional (VTC) files
    vtc_files = dir(fullfile(output_dir,'*.vtc'));
    vtc_files = {vtc_files.name}; for i=1:length(vtc_files), vtc_files{i} = fullfile(output_dir, vtc_files{i}); end
    
    % Loading the experimental protocol (PRT) files
    prt_files = dir(fullfile(output_dir,'*questions*.prt'));
    prt_files = {prt_files.name}; for i=1:length(prt_files), prt_files{i} = fullfile(output_dir, prt_files{i}); end
    
    % Loading the motion parameters files, to be used as additional GLM regressors
    motion_sdm_files = dir(fullfile(output_dir,'*3DMC.sdm'));
    motion_sdm_files = {motion_sdm_files.name}; for i=1:length(motion_sdm_files), motion_sdm_files{i} = fullfile(output_dir, motion_sdm_files{i}); end
    
    % Creating the design matrix (SDM) file for each run: each person (name) is a
    % different predictor (regressor)
    for i=1:length(vtc_files)
        vtc_filename = vtc_files{i}; prt_filename = prt_files{i};
        SDM_name = fullfile(output_dir, [subj '_questions' num2str(i) '.sdm']);
        vmr.LinkVTC(vtc_filename);
        vmr.LinkStimulationProtocol(prt_filename);
        
        vmr.ClearDesignMatrix;
        vmr.AddPredictor('person1'); vmr.AddPredictor('person2'); vmr.AddPredictor('person3'); vmr.AddPredictor('person4');
        vmr.AddPredictor('person5'); vmr.AddPredictor('person6'); vmr.AddPredictor('person7'); vmr.AddPredictor('person8');
        vmr.AddPredictor('person9'); vmr.AddPredictor('person10'); vmr.AddPredictor('person11'); vmr.AddPredictor('person12');
        vmr.AddPredictor('person13'); vmr.AddPredictor('person14'); vmr.AddPredictor('person15'); vmr.AddPredictor('person16');
        vmr.AddPredictor('person17'); vmr.AddPredictor('person18'); vmr.AddPredictor('person19'); vmr.AddPredictor('person20');
        vmr.AddPredictor('person21'); vmr.AddPredictor('person22'); vmr.AddPredictor('person23'); vmr.AddPredictor('person24');
        vmr.AddPredictor('q2_person1'); vmr.AddPredictor('q2_person2'); vmr.AddPredictor('q2_person3'); vmr.AddPredictor('q2_person4');
        vmr.AddPredictor('q2_person5'); vmr.AddPredictor('q2_person6'); vmr.AddPredictor('q2_person7'); vmr.AddPredictor('q2_person8');
        vmr.AddPredictor('q2_person9'); vmr.AddPredictor('q2_person10'); vmr.AddPredictor('q2_person11'); vmr.AddPredictor('q2_person12');
        vmr.AddPredictor('q2_person13'); vmr.AddPredictor('q2_person14'); vmr.AddPredictor('q2_person15'); vmr.AddPredictor('q2_person16');
        vmr.AddPredictor('q2_person17'); vmr.AddPredictor('q2_person18'); vmr.AddPredictor('q2_person19'); vmr.AddPredictor('q2_person20');
        vmr.AddPredictor('q2_person21'); vmr.AddPredictor('q2_person22'); vmr.AddPredictor('q2_person23'); vmr.AddPredictor('q2_person24');
        vmr.AddPredictor('Instructions');
        vmr.SetPredictorValuesFromCondition('person1', 'person1', 1.0); vmr.SetPredictorValuesFromCondition('person2', 'person2', 1.0); vmr.SetPredictorValuesFromCondition('person3', 'person3', 1.0); vmr.SetPredictorValuesFromCondition('person4', 'person4', 1.0);
        vmr.SetPredictorValuesFromCondition('person5', 'person5', 1.0); vmr.SetPredictorValuesFromCondition('person6', 'person6', 1.0); vmr.SetPredictorValuesFromCondition('person7', 'person7', 1.0); vmr.SetPredictorValuesFromCondition('person8', 'person8', 1.0);
        vmr.SetPredictorValuesFromCondition('person9', 'person9', 1.0); vmr.SetPredictorValuesFromCondition('person10', 'person10', 1.0); vmr.SetPredictorValuesFromCondition('person11', 'person11', 1.0); vmr.SetPredictorValuesFromCondition('person12', 'person12', 1.0);
        vmr.SetPredictorValuesFromCondition('person13', 'person13', 1.0); vmr.SetPredictorValuesFromCondition('person14', 'person14', 1.0); vmr.SetPredictorValuesFromCondition('person15', 'person15', 1.0); vmr.SetPredictorValuesFromCondition('person16', 'person16', 1.0);
        vmr.SetPredictorValuesFromCondition('person17', 'person17', 1.0); vmr.SetPredictorValuesFromCondition('person18', 'person18', 1.0); vmr.SetPredictorValuesFromCondition('person19', 'person19', 1.0); vmr.SetPredictorValuesFromCondition('person20', 'person20', 1.0);
        vmr.SetPredictorValuesFromCondition('person21', 'person21', 1.0); vmr.SetPredictorValuesFromCondition('person22', 'person22', 1.0); vmr.SetPredictorValuesFromCondition('person23', 'person23', 1.0); vmr.SetPredictorValuesFromCondition('person24', 'person24', 1.0);
        vmr.SetPredictorValuesFromCondition('q2_person1', 'q2_person1', 1.0); vmr.SetPredictorValuesFromCondition('q2_person2', 'q2_person2', 1.0); vmr.SetPredictorValuesFromCondition('q2_person3', 'q2_person3', 1.0); vmr.SetPredictorValuesFromCondition('q2_person4', 'q2_person4', 1.0);
        vmr.SetPredictorValuesFromCondition('q2_person5', 'q2_person5', 1.0); vmr.SetPredictorValuesFromCondition('q2_person6', 'q2_person6', 1.0); vmr.SetPredictorValuesFromCondition('q2_person7', 'q2_person7', 1.0); vmr.SetPredictorValuesFromCondition('q2_person8', 'q2_person8', 1.0);
        vmr.SetPredictorValuesFromCondition('q2_person9', 'q2_person9', 1.0); vmr.SetPredictorValuesFromCondition('q2_person10', 'q2_person10', 1.0); vmr.SetPredictorValuesFromCondition('q2_person11', 'q2_person11', 1.0); vmr.SetPredictorValuesFromCondition('q2_person12', 'q2_person12', 1.0);
        vmr.SetPredictorValuesFromCondition('q2_person13', 'q2_person13', 1.0); vmr.SetPredictorValuesFromCondition('q2_person14', 'q2_person14', 1.0); vmr.SetPredictorValuesFromCondition('q2_person15', 'q2_person15', 1.0); vmr.SetPredictorValuesFromCondition('q2_person16', 'q2_person16', 1.0);
        vmr.SetPredictorValuesFromCondition('q2_person17', 'q2_person17', 1.0); vmr.SetPredictorValuesFromCondition('q2_person18', 'q2_person18', 1.0); vmr.SetPredictorValuesFromCondition('q2_person19', 'q2_person19', 1.0); vmr.SetPredictorValuesFromCondition('q2_person20', 'q2_person20', 1.0);
        vmr.SetPredictorValuesFromCondition('q2_person21', 'q2_person21', 1.0); vmr.SetPredictorValuesFromCondition('q2_person22', 'q2_person22', 1.0); vmr.SetPredictorValuesFromCondition('q2_person23', 'q2_person23', 1.0); vmr.SetPredictorValuesFromCondition('q2_person24', 'q2_person24', 1.0);
        vmr.SetPredictorValuesFromCondition('Instructions', 'Instructions', 1.0);
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person1'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person2'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person3'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person4');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person5'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person6'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person7'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person8');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person9'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person10'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person11'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person12');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person13'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person14'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person15'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person16');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person17'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person18'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person19'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person20');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person21'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person22'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person23'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person24');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person1'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person2'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person3'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person4');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person5'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person6'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person7'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person8');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person9'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person10'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person11'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person12');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person13'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person14'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person15'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person16');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person17'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person18'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person19'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person20');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person21'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person22'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person23'); vmr.ApplyHemodynamicResponseFunctionToPredictor('q2_person24');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('Instructions');
        vmr.SDMContainsConstantPredictor = 0;
        vmr.SaveSingleStudyGLMDesignMatrix(SDM_name);
        
        % Adding motion parameters to the design matrix to clean the data of motion-related activations - 24 motion parameters (Friston 1996 model)
        sdm_inc_confounds = xff(SDM_name);
        sdm_motion = xff(motion_sdm_files{i});
        sdm_inc_confounds.NrOfPredictors = sdm_inc_confounds.NrOfPredictors + 24;
        sdm_inc_confounds.PredictorColors = [sdm_inc_confounds.PredictorColors; repmat(sdm_motion.PredictorColors, 4, 1)];
        sdm_inc_confounds.PredictorNames = [sdm_inc_confounds.PredictorNames, repmat(sdm_motion.PredictorNames, 1, 4)];
        % Calculating the 24 motion related parameters (Friston 24-parameter autoregressive model: 6 motion parameters, the motion parameters at the previous time point, and the squares of these 12 values - Friston et al. 1996, Movement-related effects in fMRI time-series, Magn Reson Med)
        motion_params = sdm_motion.SDMMatrix;
        motion_params = [motion_params, [zeros(1,6); motion_params(1:end-1,:)]];    % Adding the temporal derivatives - equal to the motion parameters at the previous time point
        motion_params = [motion_params, motion_params.^2];                          % Adding the square values of the motion parameters and their derivatives
        sdm_inc_confounds.SDMMatrix = [sdm_inc_confounds.SDMMatrix motion_params];
        % Saving the new SDM
        sdm_inc_confounds.SaveAs(SDM_name);
        
        % Checking if any run has excessive motion (>3mm)
        if ~isempty(find(sdm_motion.SDMMatrix > 3, 1))
            disp(['subject has excessive motion in run i, of ' num2str(max(sdm_motion.SDMMatrix(:))) ' mm!']);
        end
    end
end


%% Computing a GLM for each run separately

disp('Creating GLMs for each run separately (for MVPA)')
for s=1:length(subject_names)
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);
    
    % Loading normalized functional data (VTC) files
    vtc_files = dir(fullfile(output_dir,'*.vtc'));
    vtc_files = {vtc_files.name}; for i=1:length(vtc_files), vtc_files{i} = fullfile(output_dir, vtc_files{i}); end
    
    % Loading design matrix (SDM) files
    sdm_files = dir(fullfile(output_dir,'*questions*.sdm'));
    sdm_files = {sdm_files.name}; for i=1:length(sdm_files), sdm_files{i} = fullfile(output_dir, sdm_files{i}); end
    
    % Loading the normalized anatomical (VMR) file
    vmr_file = dir(fullfile(output_dir, '*MNI.vmr'));
    vmr_file = fullfile(output_dir, vmr_file(1).name);
    vmr = bvqx.OpenDocument(vmr_file);
    
    vmr.CorrectForSerialCorrelations = 1;
    
    for i = 1:length(vtc_files)
        % Compute GLM for each SDM file
        vtc_filename = vtc_files{i};
        vmr.LinkVTC(vtc_filename);
        
        sdm_filename = sdm_files{i};
        disp([vtc_filename ', ' sdm_filename]);
        
        vmr.LoadSingleStudyGLMDesignMatrix(sdm_filename);
        vmr.ComputeSingleStudyGLM;
        vmr.SaveGLM(fullfile(output_dir, [subj '_glm_questions' num2str(i) '.glm']));
    end
end



%% Computing from each of the GLMs the t values for the relevant regressors using NeuroElf, and saving to a vmp file
disp('Creating a t-map from each GLM, for each predictor (for MVPA)')
for s = 1:length(subject_names)
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);
    
    % Identifying the GLM files
    glm_files = dir(fullfile(output_dir,'*_glm_questions*.glm'));
    glm_files = {glm_files.name}; for i = 1:length(glm_files), glm_files{i} = fullfile(output_dir, glm_files{i}); end
    
    num_conditions = 48;    % 24 names of people provided by the subject, 2 questions
    
    for i = 1:length(glm_files)
        try
            current_glm = xff(glm_files{i});
            
            % Creating a contrast for each of the first 48 conditions, ignoring the other conditions (confounds predictors, e.g. motion)
            current_contrast = zeros(current_glm.NrOfPredictors, num_conditions);
            for j = 1:num_conditions, current_contrast(j, j) = 1; end   % Putting 1 for each predictor in a separate contrast, resulting in a t-value for that predictor
            
            current_t_map = current_glm.FFX_tMap(current_contrast);
            
            % Saving the t-map
            current_t_map.SaveAs(fullfile(output_dir, [subj '_t_map_questions' num2str(i) '.vmp']));
            
            current_glm.ClearObject();
            current_t_map.ClearObject();
        catch
            disp(["subject " subj " has problematic glm run " num2str(i)]);
            % IF THERE ARE PROBLEMATIC RUNS - NEED TO CALCULATE THE T-MAPS
            % MANUALLY IN BRAINVOYAGER (BrainVoyager bug)
        end
    end
end




%% Performing RSA in ROIs, and checking the patterns' consistency across questions (instead of across runs) in ROIs

% Loading the ROIs
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
AAL_atlas = reslice_data('C:\Users\michaelpeer1\Dropbox\Michael_scripts\Templates\AAL_61x73x61_YCG.nii', filename_for_resampling, 0);
ROIs{counter} = single(AAL_atlas == 37); ROI_names{counter} = 'lHC'; counter = counter + 1;
ROIs{counter} = single(AAL_atlas == 38); ROI_names{counter} = 'rHC'; counter = counter + 1;
% Orientation study ROI
ROIs{counter} = reslice_data('C:\Users\michaelpeer1\Dropbox\Michael_scripts\Templates\Orientation_paper_activations\social_vs_control_and_rest_ROI.nii', filename_for_resampling, 0); ROI_names{counter} = 'orientation_social'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Users\michaelpeer1\Dropbox\Michael_scripts\Templates\Orientation_paper_activations\space_vs_control_and_rest_ROI.nii', filename_for_resampling, 0); ROI_names{counter} = 'orientation_space'; counter = counter + 1;
% Adding the ROI of the region identified in the searchlight - only for the pattern reliability analysis
ROIs{counter} = reslice_data('C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_facebook_percentcommonfriends_radius3_othermats_regressed_percent_iters10000.nii', filename_for_resampling, 0); ROIs{end} = ROIs{end}>=2.33; ROI_names{counter} = 'searchlight_facebookdist_regressed'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_self_proximity_radius3_othermats_regressed_percent_iters10000.nii', filename_for_resampling, 0); ROIs{end} = ROIs{end}>=2.33; ROI_names{counter} = 'searchlight_selfproximity_regressed'; counter = counter + 1;
ROIs{counter} = reslice_data('C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results\group_analysis_RSA_searchlight_personality_radius3_othermats_regressed_percent_iters10000.nii', filename_for_resampling, 0); ROIs{end} = ROIs{end}>=2.33; ROI_names{counter} = 'searchlight_personality_regressed'; counter = counter + 1;
% Write all ROIs to files
ROIs_directory = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\ROIs';
for i = 1:length(ROIs)
    save_mat_to_nifti(filename_for_resampling, ROIs{i}, fullfile(ROIs_directory, [ROI_names{i} '.nii']));
end

% Loading the experimental data
all_subjects_GLM_files_dir = 'E:\Social_networks_analysis\Analysis_2mm_smoothing';          % The per-run GLM results
all_subjects_RDMs_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\Dissimilarity_mats';      % The dissimilarity matrices
all_subjects_analysis_results_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_results';          % Directory for outlut
all_subjects_stimuli_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\Experiment_stimuli';
% Finding the subject names (directories)
subject_names = dir(all_subjects_GLM_files_dir); subject_names = subject_names(3:end); subject_names = subject_names([subject_names.isdir]);
% Defining the number of stimuli (individuals)
groups_num = 4; people_in_group_num = 6;
num_conditions = groups_num * people_in_group_num;

% Computing pattern consistency within ROI
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
        current_subject_stimuli_dir = fullfile(all_subjects_stimuli_dir, subj);
        
        % Reading the data for analysis
        data_files = dir(fullfile(current_subject_GLM_dir, '*t_map_questions*.vmp'));       % GLM files - t values
        
        % Reading GLM data from the specific ROI and creating dataset
        all_ds = cell(1, length(data_files));
        for i = 1:length(data_files)
            ds_current = cosmo_fmri_dataset(fullfile(current_subject_GLM_dir, data_files(i).name), 'mask', current_ROI_mask_filename);  % Reading the current GLM file to CosmoMVPA dataset format
            ds_current = cosmo_slice(ds_current, 1:num_conditions*2, 1);          % Getting only the first 48 conditions (the run's two questions) - the rest are not interesting (motion parameters, etc.)
            ds_current.sa.targets = [1:num_conditions, 1:num_conditions]';                        % Setting the targets (condition numbers) to be 1-24
            ds_current.sa.chunks = [ones(num_conditions, 1) * (i*2-1); ones(num_conditions,1)*(i*2)];                 % Setting the chunks (independent data runs) to be the run number
            all_ds{i} = ds_current;
        end
        ds = cosmo_stack(all_ds);                   % Combine all runs to a single dataset
        ds = cosmo_remove_useless_data(ds);         % Remove voxels with no information from the data
        
        % Loading the dissimilarity matrices
        temp1 = load(fullfile(current_subject_RDMs_dir, 'facebook_distance_percentcommonfriends_RDM.mat'));
        temp3 = load(fullfile(current_subject_RDMs_dir, 'subject_responses_RDMs.mat'));
        temp_counter = 1;
        RDMs{temp_counter} = temp1.facebook_distance_RDM;          RDM_names{temp_counter} = 'facebook_distance_percentcommonfriends_RDM'; temp_counter = temp_counter + 1;
        RDMs{temp_counter} = temp3.responses_self_proximity_RDM;   RDM_names{temp_counter} = 'responses_self_proximity_RDM'; temp_counter = temp_counter + 1;
        RDMs{temp_counter} = temp3.responses_personality_RDM;      RDM_names{temp_counter} = 'responses_personality_RDM'; temp_counter = temp_counter + 1;
        RDMs{temp_counter} = temp3.responses_appearance_RDM;       RDM_names{temp_counter} = 'responses_appearance_RDM'; temp_counter = temp_counter + 1;
        % Removing places that are NaN values in one of the matrices
        sum_RDMs = zeros(size(RDMs{1})); for i = 1:length(RDMs), sum_RDMs = sum_RDMs + RDMs{i}; end   % Creating a matrix that is a sum of the others, to identify locations with nan values in one of the matrices
        for i = 1:length(RDMs), RDMs{i}(isnan(sum_RDMs)) = nan; end
        
        % Figuring out the questions' order
        current_subject_stimuli_files = dir(fullfile(current_subject_stimuli_dir, 'PEOPLE*.sce'));
        questions_order = [];
        for i = 1:length(current_subject_stimuli_files)
            f = strfind(current_subject_stimuli_files(i).name,'_');
            q1 = str2double(current_subject_stimuli_files(i).name(f(end-1)+1:f(end)-1));
            q2 = str2double(current_subject_stimuli_files(i).name(f(end)+1:end-4));
            questions_order = [questions_order, q1, q2];
        end
                
        % Running the RSA for each dissimilarity matrix, with or without regression of other matrices from it
        if sum(RDMs{2}(:))~=0       % Making sure there's no matrix that is all zeros - subject's responses were not recorded
            for rdm = 1:length(RDMs)
                % Defining the current dissimilarity matrix
                current_RDM = RDMs{rdm};
                if nansum(current_RDM(:)) ~= 0     % Ignoring empty matrices - subjects for which no data was collected, or another problem
                    % Averaging patterns across runs to get one pattern for
                    % each individual, EXCLUDING THE QUESTIONS ON THIS
                    % SPECIFIC FACTOR
                    num_conditions = length(unique(ds.sa.targets));
                    ds_average_runs = ds;
                    
                    % Removing the irrelevant questions
                    if rdm == 1      % Facebook distance - keeping all questions
                        questions_to_retain = 1:12;
                    elseif rdm == 2      % Self-proximity - questions 1,2,3,4
                        questions_to_retain = find(questions_order>4);
                    elseif rdm == 3      % Personality - questions 5,6,7,8
                        questions_to_retain = find(questions_order<=4 | questions_order>=9);
                    elseif rdm == 4      % Appearance - questions 9,10,11,12
                        questions_to_retain = find(questions_order<9);                       
                    end
                    samples_to_retain = [];
                    for i = 1:length(questions_to_retain)
                        samples_to_retain = [samples_to_retain, questions_to_retain(i)*24-23:questions_to_retain(i)*24];
                    end
                    ds_average_runs = cosmo_slice(ds_average_runs, samples_to_retain);
                    
                    % Averaging across the remaining patterns for each individual
                    for i = 1:num_conditions
                        ds_average_runs.samples(i,:) = mean(ds.samples(i:num_conditions:size(ds.samples,1),:));
                    end
                    ds_average_runs = cosmo_slice(ds_average_runs, 1:num_conditions);
                    
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
        num_runs = length(data_files) * 2;
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



