% This is performed through the BrainVoyager 20.6 scripting platform, along with NeuroElf v1.1
% (after running the preprocessing pipeline script)
%
% Each subject directory should contain VMR (anatomical) and FMR (functional) files
% In addition, you need to copy the relevant PRT (experiment protocol) files to the subject directory, which should be named in the order of the functional files


%% Enabling brainvoyager scripting from path
% bvqx = actxserver('BrainVoyagerQX.BrainVoyagerQXScriptAccess.1');     % old brainvoyager command (before BV20.4)
bvqx = actxserver('BrainVoyager.BrainVoyagerScriptAccess.1');



%% Defining working directories and getting subjects' directories names
subjects_output_data_dir = 'F:\Scales_of_space_project\Preprocessed_data\';     % the directory where subjects' preprocessed data will be saved
% getting subjects' directories names
subject_names = dir(subjects_output_data_dir); subject_names = subject_names(3:end);


%% MANUAL STAGE: EXTRACT CSF AND WM masks
% 1. White matter mask: in BrainVoyager, under segmentation tools (under 3D 
% volume tools), choose first intensity values of min=150, max=255, grow 
% region, options, define VOI, name: WM_mask.
% 2. CSF mask: choose intensity values of min=0, max=10, bounding box: x
% between 90-160, y between 100-200, z between 100-150. Grow region,
% options, define VOI, name: CSF_mask.
% 3. Save both VOIs as "WM_CSF_vois.voi".


%% Creating single study design matrices (SDM files) for each functional run
disp('Creating single study design matrices (.sdm files)')
for s = 1:length(subject_names)
    disp(s);
    current_subj = subject_names(s).name;
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    
    % loading the normalized anatomical (VMR) file
    vmr_file = dir(fullfile(output_dir, '*MNI.vmr'));
    vmr_file = fullfile(output_dir, vmr_file(1).name);
    vmr = bvqx.OpenDocument(vmr_file);
    
    % loading the normalized functional (VTC) files
    vtc_files = dir(fullfile(output_dir,'*.vtc'));
    vtc_files = {vtc_files.name};
    for i=1:length(vtc_files)
        vtc_files{i} = fullfile(output_dir, vtc_files{i});
    end
    
    % loading PRT files (need to create them separately using another script, and manually copy them to the subject directory)
    prt_files = dir(fullfile(output_dir, '*.prt'));
    prt_files = {prt_files.name}; 
    for i=1:length(prt_files)
        prt_files{i} = fullfile(output_dir, prt_files{i});
    end
    prt_files(cellfun(@(x) ~isempty(strfind(x,'full')), prt_files))=[];   % remove PRTs with all conditions, if existing
    prt_files(cellfun(@(x) ~isempty(strfind(x,'w_confounds')), prt_files))=[];   % remove PRTs with confounds, if existing
    
    % loading motion correction parameters files (created in the motion correction stage)
    motion_sdm_files = dir(fullfile(output_dir,'*3DMC.sdm'));
    motion_sdm_files = {motion_sdm_files.name}; 
    for i=1:length(motion_sdm_files)
        motion_sdm_files{i} = fullfile(output_dir, motion_sdm_files{i});
    end
    
    % Loading the WM and CSF VOI masks
    WM_CSF_mask = xff(fullfile(output_dir,'WM_CSF_vois.voi'));    
    
    % create SDM for each functional run
    for i=1:length(vtc_files)
        vtc_filename = vtc_files{i}; 
        prt_filename = prt_files{i};
        output_SDM_name = [vtc_filename(1:strfind(vtc_filename, '_SCCTBL_3DMCTS')-1) '.sdm'];
        
        vmr.LinkVTC(vtc_filename);  % linking current VTC file to the VMR
        vmr.LinkStimulationProtocol(prt_filename);  % linking current PRT file to the VMR
        disp(vtc_filename), disp(prt_filename)
        
        % Defining the design matrix parameters with the six predictors of scales
        vmr.ClearDesignMatrix;
        vmr.AddPredictor('room'); vmr.AddPredictor('building'); vmr.AddPredictor('neighborhood');
        vmr.AddPredictor('city'); vmr.AddPredictor('country'); vmr.AddPredictor('continent');
        vmr.SetPredictorValuesFromCondition('room', 'room', 1.0); vmr.SetPredictorValuesFromCondition('building', 'building', 1.0); vmr.SetPredictorValuesFromCondition('neighborhood', 'neighborhood', 1.0);
        vmr.SetPredictorValuesFromCondition('city', 'city', 1.0); vmr.SetPredictorValuesFromCondition('country', 'country', 1.0); vmr.SetPredictorValuesFromCondition('continent', 'continent', 1.0);
        vmr.ApplyHemodynamicResponseFunctionToPredictor('room'); vmr.ApplyHemodynamicResponseFunctionToPredictor('building'); vmr.ApplyHemodynamicResponseFunctionToPredictor('neighborhood');
        vmr.ApplyHemodynamicResponseFunctionToPredictor('city'); vmr.ApplyHemodynamicResponseFunctionToPredictor('country'); vmr.ApplyHemodynamicResponseFunctionToPredictor('continent');
        vmr.SDMContainsConstantPredictor = 0;
        vmr.SaveSingleStudyGLMDesignMatrix(output_SDM_name);    % Saving the resulting SDM files
        
        % adding motion correction parameters to the SDM as confound predictors
        sdm_inc_motion = xff(output_SDM_name);   % loading the SDM using neuroelf
        sdm_motion = xff(motion_sdm_files{i});      % loading the motion parameters using neuroelf
        sdm_inc_motion.NrOfPredictors = sdm_inc_motion.NrOfPredictors + sdm_motion.NrOfPredictors;
        sdm_inc_motion.PredictorColors = [sdm_inc_motion.PredictorColors; sdm_motion.PredictorColors];
        sdm_inc_motion.PredictorNames = [sdm_inc_motion.PredictorNames sdm_motion.PredictorNames];
        sdm_inc_motion.SDMMatrix = [sdm_inc_motion.SDMMatrix sdm_motion.SDMMatrix];
        
        % Adding the WM and CSF time courses to the SDM as confound predictors
        vtc = xff(vtc_filename);
        WM_CSF_timecourses = vtc.VOITimeCourse(WM_CSF_mask);
        sdm_inc_motion.NrOfPredictors = sdm_inc_motion.NrOfPredictors + 2;
        sdm_inc_motion.PredictorColors = [sdm_inc_motion.PredictorColors; sdm_inc_motion.PredictorColors(end-1:end, :)];
        sdm_inc_motion.PredictorNames = [sdm_inc_motion.PredictorNames, 'WM_mean', 'CSF_mean'];
        sdm_inc_motion.SDMMatrix = [sdm_inc_motion.SDMMatrix zscore(WM_CSF_timecourses)];

        % saving the SDM file again, this time with the motion parameters and WM and CSF predictors included
        sdm_inc_motion.SaveAs(output_SDM_name);  
    end
end



%% Creating a multi-study design matrix (MDM) from all runs of the subject, and computing the corresponding individual-subject GLM
disp('Creating multi-study design matrices (.mdm files) and GLM files')
for s=1:length(subject_names)
    disp(s);
    current_subj = subject_names(s).name;
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    
    % loading the normalized functional (VTC) files
    vtc_files = dir(fullfile(output_dir,'*.vtc'));
    vtc_files = {vtc_files.name}; for i=1:length(vtc_files), vtc_files{i} = fullfile(output_dir, vtc_files{i}); end
    
    % loading the normalized anatomical (VMR) file
    vmr_file = dir(fullfile(output_dir, '*MNI.vmr'));
    vmr_file = fullfile(output_dir, vmr_file(1).name);
    vmr = bvqx.OpenDocument(vmr_file);
    
    % Creating MDM file for each subject - all runs combined, but without the control run
    vmr.ClearMultiStudyGLMDefinition;
    for i=1:length(vtc_files)
        vtc_filename = vtc_files{i};
        current_SDM_name = [vtc_filename(1:strfind(vtc_filename, '_SCCTBL_3DMCTS')-1) '.sdm'];
        if isempty(strfind(vtc_filename, 'control'))     % checking that it is not a control run
            vmr.AddStudyAndDesignMatrix(vtc_filename, current_SDM_name);
        end
    end
    % saving MDM
    MDM_name = fullfile(output_dir, [current_subj '_MNI_combined_no_control.mdm']);
    vmr.SaveMultiStudyGLMDefinitionFile(MDM_name);
    
    % computing and saving the GLM from all the subject's combined runs (without control run)
    GLM_name = fullfile(output_dir, [current_subj '_MNI_combined_no_control.glm']);
    % defining GLM estimation parameters
    vmr.SeparationOfStudyPredictors = 0;
    vmr.CorrectForSerialCorrelations = 1;
    vmr.PSCTransformStudies = 1;
    % computing and saving the GLM
    vmr.ComputeMultiStudyGLM;
    vmr.SaveGLM(GLM_name);
    
    
    % Creating another MDM file including the control run, with separate study predictors
    for i=1:length(vtc_files)
        if ~isempty(strfind(vtc_files{i},'control'))
            % adding the control run to the existing MDM
            vtc_filename = vtc_files{i};
            current_SDM_name = [vtc_filename(1:strfind(vtc_filename, '_SCCTBL_3DMCTS')-1) '.sdm'];
            vmr.AddStudyAndDesignMatrix(vtc_filename, current_SDM_name);
        end
    end
    % saving MDM with control run
    MDM_name = fullfile(output_dir, [current_subj '_MNI_separate_study_predictors.mdm']);
    vmr.SaveMultiStudyGLMDefinitionFile(MDM_name);
    
end



%% Performing a group analysis

% defining a directory for the group analysis results
subjects_group_analysis_dir = 'F:\Scales_of_space_project\Analysis_group\';

% creating group MDM file with all runs of all subjects, using NeuroElf
MDM_group = xff('new:mdm');
MDM_group.SeparatePredictors = 2;
% Adding data from each subjects' MDM
for s=1:length(subject_names)
    current_subj = subject_names(s).name;
    current_data_dir = fullfile(subjects_output_data_dir, current_subj);
    % loading subject-specific MDM file (including all runs of this subject)
    MDM_current = xff(fullfile(current_data_dir, [current_subj '_MNI_combined_no_control.mdm']));
    % adding to group MDM file
    MDM_group.XTC_RTC = [MDM_group.XTC_RTC; MDM_current.XTC_RTC];
end
% saving the new group MDM
MDM_name = fullfile(subjects_group_analysis_dir, 'all_subjs_MNI_combined_no_control.mdm');
MDM_group.SaveAs(MDM_name);


% COMPUTE THE GLM MANUALLY USING BRAINVOYAGER
% (name - 'all_subjs_MNI_combined_no_control.glm', CorrectForSerialCorrelations, PSCTransformStudies, RFX)




%% %% Linking VTC and PRT files, for event-related averaging

disp('Linking VTC and PRT files')
for s=1:length(subject_names)
    subj = subject_names(s).name;
    disp(subj);
    output_dir = fullfile(subjects_output_data_dir, subj);
    
    vtc_files = dir(fullfile(output_dir,'*.vtc'));
    vtc_files = {vtc_files.name}; for i=1:length(vtc_files), vtc_files{i} = fullfile(output_dir, vtc_files{i}); end

    prt_files = dir(fullfile(output_dir, '*.prt'));
    prt_files = {prt_files.name}; for i=1:length(prt_files), prt_files{i} = fullfile(output_dir, prt_files{i}); end
    prt_files(cellfun(@(x) ~isempty(strfind(x,'full')), prt_files)) = [];   % remove PRTs with all conditions, if existing
    prt_files(cellfun(@(x) ~isempty(strfind(x,'w_confounds')), prt_files))=[];   % remove PRTs with confounds, if existing
    
    disp([length(vtc_files) length(prt_files)])
    
    for i=1:length(vtc_files)
        curr_vtc = xff(vtc_files{i});
        curr_vtc.NrOfLinkedPRTs = 1;
        curr_vtc.NameOfLinkedPRT = prt_files{i};
        curr_vtc.NameOfLinkedPRT(strfind(curr_vtc.NameOfLinkedPRT,'\')) = '/';
        curr_vtc.Save;
        curr_vtc.ClearObject;
    end
end

% Manual - create event-related averaging using the MDM of all subjects

