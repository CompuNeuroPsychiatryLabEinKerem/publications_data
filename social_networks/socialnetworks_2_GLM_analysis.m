% Running a GLM analysis in Brainvoyager
% Each person has one regressor, including all of this person's name
% appearances during the task


%% Starting brainvoyager from Matlab
bvqx = actxserver('BrainVoyager.BrainVoyagerScriptAccess.1');


%% Defining working directories
subjects_analysis_dir = 'F:\Social_networks_analysis\Analysis_2mm_smoothing\';          % A directory where the preprocessed data will be put
logfiles_parent_dir = 'F:\Social_networks_analysis\Experiments_stimuli\';               % A directory containing the logfiles from each subject
subject_names = dir(subjects_analysis_dir); subject_names = subject_names(3:end);       % Getting the names of the subjects directories



%% Creating BrainVoyager protocol (PRT) files for each subject, from the experiment log files

template_prt_file = 'F:\Social_networks_analysis\Experiments_stimuli\General_stimuli\empty_prt_individuals.prt';   % an empty PRT files to be used as a template - containing 4 social groups, 6 people in each

for s = 1:length(subject_names)
    subj = subject_names(s).name;
    disp(subj);
    current_logfile_dir = fullfile(logfiles_parent_dir, subj);
    current_subj_analysis_dir = fullfile(subjects_analysis_dir, subj);
    
    current_logfiles = dir(fullfile(current_logfile_dir, '*.log'));     % All of the subjects' logfiles
    for i = 1:length(current_logfiles)
        current_logfile = fullfile(current_logfiles(i).folder, current_logfiles(i).name);
        new_prt_name_individuals = fullfile(current_subj_analysis_dir, [subj '_individuals' num2str(i) '.prt']);
        % Converting the log files to PRT files
        socialnetworks_convert_logfile_to_prt(current_logfile, new_prt_name_individuals, template_prt_file);
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
    prt_files = dir(fullfile(output_dir,'*individuals*.prt'));
    prt_files = {prt_files.name}; for i=1:length(prt_files), prt_files{i} = fullfile(output_dir, prt_files{i}); end
    
    % Loading the motion parameters files, to be used as additional GLM regressors
    motion_sdm_files = dir(fullfile(output_dir,'*3DMC.sdm'));
    motion_sdm_files = {motion_sdm_files.name}; for i=1:length(motion_sdm_files), motion_sdm_files{i} = fullfile(output_dir, motion_sdm_files{i}); end
    
    % Creating the design matrix (SDM) file for each run: each person (name) is a
    % different predictor (regressor)
    for i=1:length(vtc_files)
        vtc_filename = vtc_files{i}; prt_filename = prt_files{i};
        SDM_name = fullfile(output_dir, [subj '_individuals' num2str(i) '.sdm']);
        vmr.LinkVTC(vtc_filename);
        vmr.LinkStimulationProtocol(prt_filename);
        
        vmr.ClearDesignMatrix;
        vmr.AddPredictor('person1'); vmr.AddPredictor('person2'); vmr.AddPredictor('person3'); vmr.AddPredictor('person4'); 
        vmr.AddPredictor('person5'); vmr.AddPredictor('person6'); vmr.AddPredictor('person7'); vmr.AddPredictor('person8');
        vmr.AddPredictor('person9'); vmr.AddPredictor('person10'); vmr.AddPredictor('person11'); vmr.AddPredictor('person12');
        vmr.AddPredictor('person13'); vmr.AddPredictor('person14'); vmr.AddPredictor('person15'); vmr.AddPredictor('person16');
        vmr.AddPredictor('person17'); vmr.AddPredictor('person18'); vmr.AddPredictor('person19'); vmr.AddPredictor('person20');
        vmr.AddPredictor('person21'); vmr.AddPredictor('person22'); vmr.AddPredictor('person23'); vmr.AddPredictor('person24');
        vmr.AddPredictor('Instructions');
        vmr.SetPredictorValuesFromCondition('person1', 'person1', 1.0); vmr.SetPredictorValuesFromCondition('person2', 'person2', 1.0); vmr.SetPredictorValuesFromCondition('person3', 'person3', 1.0); vmr.SetPredictorValuesFromCondition('person4', 'person4', 1.0); 
        vmr.SetPredictorValuesFromCondition('person5', 'person5', 1.0); vmr.SetPredictorValuesFromCondition('person6', 'person6', 1.0); vmr.SetPredictorValuesFromCondition('person7', 'person7', 1.0); vmr.SetPredictorValuesFromCondition('person8', 'person8', 1.0);
        vmr.SetPredictorValuesFromCondition('person9', 'person9', 1.0); vmr.SetPredictorValuesFromCondition('person10', 'person10', 1.0); vmr.SetPredictorValuesFromCondition('person11', 'person11', 1.0); vmr.SetPredictorValuesFromCondition('person12', 'person12', 1.0);
        vmr.SetPredictorValuesFromCondition('person13', 'person13', 1.0); vmr.SetPredictorValuesFromCondition('person14', 'person14', 1.0); vmr.SetPredictorValuesFromCondition('person15', 'person15', 1.0); vmr.SetPredictorValuesFromCondition('person16', 'person16', 1.0);
        vmr.SetPredictorValuesFromCondition('person17', 'person17', 1.0); vmr.SetPredictorValuesFromCondition('person18', 'person18', 1.0); vmr.SetPredictorValuesFromCondition('person19', 'person19', 1.0); vmr.SetPredictorValuesFromCondition('person20', 'person20', 1.0);
        vmr.SetPredictorValuesFromCondition('person21', 'person21', 1.0); vmr.SetPredictorValuesFromCondition('person22', 'person22', 1.0); vmr.SetPredictorValuesFromCondition('person23', 'person23', 1.0); vmr.SetPredictorValuesFromCondition('person24', 'person24', 1.0);
        vmr.SetPredictorValuesFromCondition('Instructions', 'Instructions', 1.0);
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person1'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person2'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person3'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person4'); 
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person5'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person6'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person7'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person8');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person9'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person10'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person11'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person12');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person13'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person14'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person15'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person16');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person17'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person18'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person19'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person20');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('person21'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person22'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person23'); vmr.ApplyHemodynamicResponseFunctionToPredictor('person24');
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
    sdm_files = dir(fullfile(output_dir,'*individuals*.sdm'));
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
        vmr.SaveGLM(fullfile(output_dir, [subj '_glm_individuals' num2str(i) '.glm']));
    end
end



%% Computing from each of the GLMs the t values for the relevant regressors using NeuroElf, and saving to a vmp file
disp('Creating a t-map from each GLM, for each predictor (for MVPA)')
for s = 1:length(subject_names)
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);
    
    % Identifying the GLM files
    glm_files = dir(fullfile(output_dir,'*_glm_individuals*.glm'));
    glm_files = {glm_files.name}; for i = 1:length(glm_files), glm_files{i} = fullfile(output_dir, glm_files{i}); end
        
    num_conditions = 24;    % 24 names of people provided by the subject
    
    for i = 1:length(glm_files)
        try
            current_glm = xff(glm_files{i});
            
            % Creating a contrast for each of the first 24 conditions, ignoring the other conditions (confounds predictors, e.g. motion)
            current_contrast = zeros(current_glm.NrOfPredictors, num_conditions);
            for j = 1:num_conditions, current_contrast(j, j) = 1; end   % Putting 1 for each predictor in a separate contrast, resulting in a t-value for that predictor
            
            current_t_map = current_glm.FFX_tMap(current_contrast);
            
            % Saving the t-map
            current_t_map.SaveAs(fullfile(output_dir, [subj '_t_map_individuals' num2str(i) '.vmp']));
            
            current_glm.ClearObject();
            current_t_map.ClearObject();
        catch
            disp(["subject " subj " has problematic glm run " num2str(i)]);
            % IF THERE ARE PROBLEMATIC RUNS - NEED TO CALCULATE THE T-MAPS
            % MANUALLY IN BRAINVOYAGER (BrainVoyager bug)
        end
    end
end
