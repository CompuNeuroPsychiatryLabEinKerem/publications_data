%% Defining working directories and getting subjects' directories names
subjects_output_data_dir = 'F:\Scales_of_space_project\Preprocessed_data\';     % the directory where subjects' preprocessed data will be saved
% getting subjects' directories names
subject_names = dir(subjects_output_data_dir); subject_names = subject_names(3:end);

% Directory with group analysis results
subjects_group_analysis_dir = 'F:\Scales_of_space_project\Analysis_group\';

% Directory with all subjects' log files
all_subjects_experiments_stimuli_dir = 'F:\Scales_of_space_project\Experiments_stimuli_and_log_files\';

% Excel file with all of the subjects' ratings of difficulty, emotional valence and familiarity for each location
excel_file_questionnaire_allsubjs = 'F:\Scales_of_space_project\behavioral_questionnaires_rescaling.xlsx';



%% Analyzing the effect of potential confounds, by creating SDM files for all subjects with potential confounds as parametrically modulated regressors
% This part uses the "create_sdm_with_other_factors.m" script

for s=1:length(subject_names)
    disp(subject_names(s).name)
    current_subject_analysis_dir = fullfile(subjects_output_data_dir, subject_names(s).name);
    current_subject_logfiles_dir = fullfile(all_subjects_experiments_stimuli_dir, subject_names(s).name);
    
    % finding the PRT files with all locations, with 'full' in their filename (12 locations, not grouped according to scales into 6 scales)
    all_locations_prt_files = dir(fullfile(current_subject_analysis_dir, '*full*.prt'));
    all_locations_prt_files = {all_locations_prt_files.name};
    for i=1:length(all_locations_prt_files)
        all_locations_prt_files{i} = fullfile(current_subject_analysis_dir, all_locations_prt_files{i});
    end
    
    % finding the SDM files created earlier, with predictors only for the six scales
    original_sdm_files = dir(fullfile(current_subject_analysis_dir,'*.sdm'));
    original_sdm_files = {original_sdm_files.name};
    for i=1:length(original_sdm_files)
        original_sdm_files{i} = fullfile(current_subject_analysis_dir, original_sdm_files{i});
    end
    original_sdm_files(cellfun(@(x) ~isempty(strfind(x,'3DMC')), original_sdm_files))=[];
    original_sdm_files(cellfun(@(x) ~isempty(strfind(x,'w_confounds')), original_sdm_files))=[];
    original_sdm_files(cellfun(@(x) ~isempty(strfind(x,'full')), original_sdm_files))=[];
    original_sdm_files(cellfun(@(x) ~isempty(strfind(x,'w_controlcond')), original_sdm_files))=[];
    original_sdm_files(cellfun(@(x) ~isempty(strfind(x,'oneconfound')), original_sdm_files))=[];
    
    % finding the appropriate log file in the experiments stimuli directory
    logfiles = dir(fullfile(current_subject_logfiles_dir, '*.log'));
    logfiles = {logfiles.name};
    for i=1:length(logfiles)
        logfiles{i} = fullfile(current_subject_logfiles_dir, logfiles{i});
    end
    logfiles(cellfun(@(x) ~isempty(strfind(x,'training')), logfiles))=[];
    logfiles(cellfun(@(x) ~isempty(strfind(x,'localizer')), logfiles))=[];
    logfiles(cellfun(@(x) ~isempty(strfind(x,'old')), logfiles))=[];
    
    % creating the SDMs including the other factors' predictors
    for i = 1:length(all_locations_prt_files)
        [~, output_file] = fileparts(original_sdm_files{i});
        output_file = fullfile(current_subject_analysis_dir, [output_file(1:11) '_w_confounds' output_file(12:end) '.sdm']);
        create_sdm_with_other_factors(all_locations_prt_files{i}, original_sdm_files{i}, logfiles{i}, excel_file_questionnaire_allsubjs, subject_names(s).name, output_file);
    end
end

% add all runs to group analysis MDM
m = xff('new:mdm');
for s=1:length(subject_names)
    current_subject_analysis_dir = fullfile(subjects_output_data_dir, subject_names(s).name);
    sdm_files = dir(fullfile(current_subject_analysis_dir,'*w_confounds*.sdm'));
    sdm_files = {sdm_files.name}; for i=1:length(sdm_files), sdm_files{i} = fullfile(current_subject_analysis_dir, sdm_files{i}); end
    sdm_files(cellfun(@(x) ~isempty(strfind(x,'control')), sdm_files))=[];   % remove control run files, if existing
    
    vtc_files = dir(fullfile(current_subject_analysis_dir,'*.vtc'));
    vtc_files = {vtc_files.name}; for i=1:length(vtc_files), vtc_files{i} = fullfile(current_subject_analysis_dir, vtc_files{i}); end
    vtc_files(cellfun(@(x) ~isempty(strfind(x,'control')), vtc_files))=[];   % remove control run files, if existing
    
    disp([length(sdm_files) length(vtc_files)])
    
    for i=1:length(sdm_files)
        m.XTC_RTC = [m.XTC_RTC; {vtc_files{i} sdm_files{i}}];
    end
end
m.NrOfStudies = size(m.XTC_RTC, 1);
m.SaveAs(fullfile(subjects_group_analysis_dir, 'all_subjs_combined_no_control_w_confounds.mdm'));

% MANUAL - run VOI GLM with the MDM, after exclusion of subjects not having some of the questionnaire responses


% Create group analysis MDM for each of the confounds separately, and for the confounds with the scales
confound_names = {'difficulty', 'emotion', 'familiarity', 'strategy_looking', 'strategy_map', 'strategy_moving', 'strategy_lines', 'strategy_feeling', 'strategy_1pp', 'strategy_3pp', 'RT'};
for c = 1:length(confound_names)
    m = xff('new:mdm');
    m_w_scales = xff('new:mdm');
    for s = 1:length(subject_names)
        current_subject_analysis_dir = fullfile(subjects_output_data_dir, subject_names(s).name);
        sdm_files = dir(fullfile(current_subject_analysis_dir, ['*_oneconfound_*', confound_names{c},'*.sdm']));        
        sdm_files_w_scales = dir(fullfile(current_subject_analysis_dir, ['*_scalesandoneconfound_*', confound_names{c},'*.sdm']));        
        if ~isempty(sdm_files)              % Checking that the data exists for this subject
            sdm_files = {sdm_files.name}; for i=1:length(sdm_files), sdm_files{i} = fullfile(current_subject_analysis_dir, sdm_files{i}); end
            sdm_files_w_scales = {sdm_files_w_scales.name}; for i=1:length(sdm_files_w_scales), sdm_files_w_scales{i} = fullfile(current_subject_analysis_dir, sdm_files_w_scales{i}); end
            sdm_files(cellfun(@(x) ~isempty(strfind(x,'control')), sdm_files))=[];   % remove control run files, if existing
            sdm_files_w_scales(cellfun(@(x) ~isempty(strfind(x,'control')), sdm_files_w_scales))=[];   % remove control run files, if existing
            
            vtc_files = dir(fullfile(current_subject_analysis_dir,'*.vtc'));
            vtc_files = {vtc_files.name}; for i=1:length(vtc_files), vtc_files{i} = fullfile(current_subject_analysis_dir, vtc_files{i}); end
            vtc_files(cellfun(@(x) ~isempty(strfind(x,'control')), vtc_files))=[];   % remove control run files, if existing
            
            disp([length(sdm_files) length(vtc_files)])
            
            for i=1:length(sdm_files)
                m.XTC_RTC = [m.XTC_RTC; {vtc_files{i} sdm_files{i}}];
                m_w_scales.XTC_RTC = [m_w_scales.XTC_RTC; {vtc_files{i} sdm_files_w_scales{i}}];
            end
        end
    end
    m.NrOfStudies = size(m.XTC_RTC, 1); m_w_scales.NrOfStudies = size(m_w_scales.XTC_RTC, 1);
    m.SaveAs(fullfile(subjects_group_analysis_dir, ['all_subjs_oneconfound_',confound_names{c},'.mdm']));
    m_w_scales.SaveAs(fullfile(subjects_group_analysis_dir, ['all_subjs_scalesandoneconfound_',confound_names{c},'.mdm']));
end

% MANUAL STAGE - COMPUTE GROUP GLM FOR EACH CONFOUND SEPARATELY



%% Analysis of activations vs. the control run activations

% Creating SDMs with control condition included
disp('Creating single study design matrices (.sdm files) with control condition')
for s=1:length(subject_names)
    disp(s);
    subj = subject_names(s).name;
    current_subject_dir = fullfile(subjects_output_data_dir, subj);
    
    % read existing SDM files
    sdm_files = dir(fullfile(current_subject_dir, '*.sdm')); sdm_files = fullfile(current_subject_dir, {sdm_files.name});
    sdm_files(cellfun(@(x) ~isempty(strfind(x,'full')), sdm_files)) = [];   % remove full sdm files, if existing
    sdm_files(cellfun(@(x) ~isempty(strfind(x,'3DMC')), sdm_files)) = [];   % remove motion sdm files, if existing
    sdm_files(cellfun(@(x) ~isempty(strfind(x,'w_controlcond')), sdm_files)) = [];   % remove existing SDM files with control condition, if existing
    sdm_files(cellfun(@(x) ~isempty(strfind(x,'w_confounds')), sdm_files)) = [];   % remove existing SDM files with control condition, if existing
    sdm_files(cellfun(@(x) ~isempty(strfind(x,'oneconfound')), sdm_files)) = [];   % remove existing SDM files with control condition, if existing
    
    % changing design matrices to include the control conditions and saving them to new SDM files
    for i=1:length(sdm_files)
        curr_sdm = xff(sdm_files{i});   % reading the current SDM file
        curr_sdm.NrOfPredictors = 18;    % 6 scales, 6 control, 6 motion
        curr_sdm.FirstConfoundPredictor = 13;
        curr_sdm.PredictorNames = [curr_sdm.PredictorNames(1:6) ...
            {'room_control' 'buil_control' 'neig_control' 'city_control' 'coun_control' 'cont_control'} ...
            curr_sdm.PredictorNames(7:12)];
        curr_sdm.PredictorColors = [curr_sdm.PredictorColors(1:6,:); curr_sdm.PredictorColors(1:6,:); curr_sdm.PredictorColors(7:12,:)];
        % if the current run is not a control run
        if ~strcmp(sdm_files{i}(end-10:end-4), 'control')
            curr_sdm.SDMMatrix = [curr_sdm.SDMMatrix(:,1:6) zeros(size(curr_sdm.SDMMatrix,1),6) curr_sdm.SDMMatrix(:,7:12)];
            curr_sdm.RTCMatrix = [curr_sdm.RTCMatrix(:,1:6) zeros(size(curr_sdm.RTCMatrix,1),6)];
            curr_sdm.SaveAs([sdm_files{i}(1:end-6) '_w_controlcond_' sdm_files{i}(end-4:end)]);
            % else if the current run is the control run
        else
            curr_sdm.SDMMatrix = [zeros(size(curr_sdm.SDMMatrix,1),6) curr_sdm.SDMMatrix(:,1:6) curr_sdm.SDMMatrix(:,7:12)];
            curr_sdm.RTCMatrix = [zeros(size(curr_sdm.RTCMatrix,1),6) curr_sdm.RTCMatrix(:,1:6)];
            curr_sdm.SaveAs([sdm_files{i}(1:end-12) '_w_controlcond_' sdm_files{i}(end-10:end)]);
        end
        curr_sdm.ClearObject;
    end
end

% Creating MDMs and GLMs, including the control run
disp('Creating multi-study design matrices (.mdm files) and GLM files')
for s=1:length(subject_names)
    disp(s);
    subj = subject_names(s).name;
    current_subject_dir = fullfile(subjects_output_data_dir, subj);
    
    % reading existing MDM file
    mdm_file = dir(fullfile(current_subject_dir, '*separate_study_predictors.mdm')); mdm_file = fullfile(current_subject_dir, mdm_file.name);
    mdm_file_w_control = [mdm_file(1:end-30) '_w_controlcond.mdm'];     % the name of the new MDM to create
    
    % reading existing SDM files with control condition
    sdm_files = dir(fullfile(current_subject_dir, '*_w_controlcond_*.sdm')); sdm_files = fullfile(current_subject_dir, {sdm_files.name});
    
    % changing MDM file contents to the SDMs with control condition, and saving
    curr_mdm = xff(mdm_file);
    curr_mdm.XTC_RTC(:,2) = sdm_files';
    curr_mdm.SaveAs(mdm_file_w_control);
    curr_mdm.ClearObject;
end


% creating group MDM, including the control condition
group_MDM_path = fullfile(subjects_group_analysis_dir, 'group_MDM_w_controlcond.mdm');
subj = subject_names(1).name; current_subject_dir = fullfile(subjects_output_data_dir, subj);
mdm_file = dir(fullfile(current_subject_dir, '*_w_controlcond.mdm')); mdm_file = fullfile(current_subject_dir, mdm_file.name);
mdm_group = xff(mdm_file);
for s=2:length(subject_names)
    disp(s);
    subj = subject_names(s).name;
    current_subject_dir = fullfile(subjects_output_data_dir, subj);
    mdm_file = dir(fullfile(current_subject_dir, '*_w_controlcond.mdm')); mdm_file = fullfile(current_subject_dir, mdm_file.name);
    curr_mdm = xff(mdm_file);
    mdm_group.XTC_RTC = [mdm_group.XTC_RTC; curr_mdm.XTC_RTC];
    curr_mdm.ClearObject;
end
mdm_group.NrOfStudies = size(mdm_group.XTC_RTC, 1);
mdm_group.SaveAs(group_MDM_path);


% CHANGE MANUALLY TO REMOVE SDM FILES OF SUBJECTS WITHOUT CONTROL RUN



%% Correlation between parameters and scales
[~, ~, xls_quest] = xlsread(excel_file_questionnaire_allsubjs);
LINES_DIFFICULTY = 3:20;
LINES_EMOTION = 25:39;
LINES_FAMILIARITY = 45:59;
LINES_STRATEGY_LOOKING = 64:81;
LINES_STRATEGY_MAP = 86:103;
LINES_STRATEGY_MOVING = 109:126;
LINES_STRATEGY_LINES = 131:148;
LINES_STRATEGY_FEELING = 153:170;
LINES_STRATEGY_1PP = 175:190;
LINES_STRATEGY_3PP = 195:210;
num_conditions = 12;


difficulty = cell2mat(xls_quest(LINES_DIFFICULTY,2:1+num_conditions));
emotion = cell2mat(xls_quest(LINES_EMOTION,2:1+num_conditions));
familiarity = cell2mat(xls_quest(LINES_FAMILIARITY,2:1+num_conditions));
strategy_1pp = cell2mat(xls_quest(LINES_STRATEGY_1PP,2:1+num_conditions));
strategy_3pp = cell2mat(xls_quest(LINES_STRATEGY_3PP,2:1+num_conditions));
strategy_feeling = cell2mat(xls_quest(LINES_STRATEGY_FEELING,2:1+num_conditions));
strategy_lines = cell2mat(xls_quest(LINES_STRATEGY_LINES,2:1+num_conditions));
strategy_looking = cell2mat(xls_quest(LINES_STRATEGY_LOOKING,2:1+num_conditions));
strategy_map = cell2mat(xls_quest(LINES_STRATEGY_MAP,2:1+num_conditions));
strategy_moving = cell2mat(xls_quest(LINES_STRATEGY_MOVING,2:1+num_conditions));

corr_difficulty = corr(difficulty', [1;1;2;2;3;3;4;4;5;5;6;6]);
corr_emotion = corr(emotion', [1;1;2;2;3;3;4;4;5;5;6;6]);
corr_familiarity = corr(familiarity', [1;1;2;2;3;3;4;4;5;5;6;6]);
corr_strategy_1pp = corr(strategy_1pp(:, 1:2:end)', [1;2;3;4;5;6]);
corr_strategy_3pp = corr(strategy_3pp(:, 1:2:end)', [1;2;3;4;5;6]);
corr_strategy_feeling = corr(strategy_feeling(:, 1:2:end)', [1;2;3;4;5;6]);
corr_strategy_lines = corr(strategy_lines(:, 1:2:end)', [1;2;3;4;5;6]);
corr_strategy_looking = corr(strategy_looking(:, 1:2:end)', [1;2;3;4;5;6]);
corr_strategy_map = corr(strategy_map(:, 1:2:end)', [1;2;3;4;5;6]);
corr_strategy_moving = corr(strategy_moving(:, 1:2:end)', [1;2;3;4;5;6]);
