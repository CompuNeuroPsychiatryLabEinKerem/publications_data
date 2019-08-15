function create_sdm_with_other_factors(original_allconds_prt, original_scales_sdm, logfile_for_RT, excel_file_questionnaire_allsubjs, subject_name, output_file_scales_confounds)
% create_sdm_with_other_factors(original_allconds_prt, original_scales_sdm, logfile_for_RT, excel_file_questionnaire_allsubjs, subject_name, output_file)
%
% Create an SDM file with all confounds as parametric modulators - response time, reported difficulty, reported emotional valence, reported familiarity
% Inputs:
% - original_allconds_prt - a "full" prt file (with all 12 locations)
% - original_scales_sdm - an existing sdm file with only 6 scales (no separation between the two locations inside scale)
% - logfile_for_RT - presentation log file of the corresponding run, to extract the response times from
% - excel_file_questionnaire_allsubjs - an excel file summarizing the results of all subjects' questionnaires (difficulty, emotional valence and familiarity)
% - subject_name, as in the sdm/prt/directory name
% - output_file - a name for the resulting SDM file with all scales and confounds (another file will be saved with the confounds only)


% read the excel file to extract current subject's rankings, and z-scoring them
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

difficulty = zeros(1, num_conditions); emotion = zeros(1, num_conditions); familiarity = zeros(1, num_conditions);
strategy_looking = zeros(1, num_conditions); strategy_map = zeros(1, num_conditions); strategy_moving = zeros(1, num_conditions);
strategy_lines = zeros(1, num_conditions); strategy_feeling = zeros(1, num_conditions); 
strategy_1pp = zeros(1, num_conditions); strategy_3pp = zeros(1, num_conditions);
for i=1:length(LINES_DIFFICULTY)
    if strcmp(subject_name, xls_quest{LINES_DIFFICULTY(i),1}(1:11))
        difficulty = zscore(cell2mat(xls_quest(LINES_DIFFICULTY(i),2:1+num_conditions)));
    end
end
for i=1:length(LINES_EMOTION)
    if strcmp(subject_name, xls_quest{LINES_EMOTION(i),1}(1:11))
        emotion = zscore(cell2mat(xls_quest(LINES_EMOTION(i),2:1+num_conditions)));
    end
end
for i=1:length(LINES_FAMILIARITY)
    if strcmp(subject_name, xls_quest{LINES_FAMILIARITY(i),1}(1:11))
        familiarity = zscore(cell2mat(xls_quest(LINES_FAMILIARITY(i),2:1+num_conditions)));
    end
end
for i=1:length(LINES_STRATEGY_LOOKING)
    if strcmp(subject_name, xls_quest{LINES_STRATEGY_LOOKING(i),1}(1:11))
        strategy_looking = zscore(cell2mat(xls_quest(LINES_STRATEGY_LOOKING(i),2:1+num_conditions)));
    end
end
for i=1:length(LINES_STRATEGY_MAP)
    if strcmp(subject_name, xls_quest{LINES_STRATEGY_MAP(i),1}(1:11))
        strategy_map = zscore(cell2mat(xls_quest(LINES_STRATEGY_MAP(i),2:1+num_conditions)));
    end
end
for i=1:length(LINES_STRATEGY_MOVING)
    if strcmp(subject_name, xls_quest{LINES_STRATEGY_MOVING(i),1}(1:11))
        strategy_moving = zscore(cell2mat(xls_quest(LINES_STRATEGY_MOVING(i),2:1+num_conditions)));
    end
end
for i=1:length(LINES_STRATEGY_LINES)
    if strcmp(subject_name, xls_quest{LINES_STRATEGY_LINES(i),1}(1:11))
        strategy_lines = zscore(cell2mat(xls_quest(LINES_STRATEGY_LINES(i),2:1+num_conditions)));
    end
end
for i=1:length(LINES_STRATEGY_FEELING)
    if strcmp(subject_name, xls_quest{LINES_STRATEGY_FEELING(i),1}(1:11))
        strategy_feeling = zscore(cell2mat(xls_quest(LINES_STRATEGY_FEELING(i),2:1+num_conditions)));
    end
end
for i=1:length(LINES_STRATEGY_1PP)
    if strcmp(subject_name, xls_quest{LINES_STRATEGY_1PP(i),1}(1:11))
        strategy_1pp = zscore(cell2mat(xls_quest(LINES_STRATEGY_1PP(i),2:1+num_conditions)));
    end
end
for i=1:length(LINES_STRATEGY_3PP)
    if strcmp(subject_name, xls_quest{LINES_STRATEGY_3PP(i),1}(1:11))
        strategy_3pp = zscore(cell2mat(xls_quest(LINES_STRATEGY_3PP(i),2:1+num_conditions)));
    end
end


% read RT from presentation log file
fid = fopen(logfile_for_RT);
text = textscan(fid,'%s');
text = text{1};
fclose(fid);
% find pulse indices
f = strfind(text,'Pulse'); indices_pulse = []; for i=1:length(f), if f{i}==1, indices_pulse = [indices_pulse i]; end,end
% find responses
RTs = nan(length(indices_pulse),1);
for i=1:length(indices_pulse)-1
    r = strfind(text(indices_pulse(i):indices_pulse(i+1)), 'Response');
    r = find(cellfun(@(x) ~isempty(x), r));
    if ~isempty(r)
        RTs(i) = str2num(text{indices_pulse(i) + r(1) + 2});    % if there are two responses, taking the first RT of the two
    end
end
RTs(~isnan(RTs)) = zscore(RTs(~isnan(RTs)));
RTs(isnan(RTs)) = 0;


% read the original PRT (with no confounds, only times for each condition)
prt = xff(original_allconds_prt);
new_prt = xff('new:prt');
new_prt.FileVersion = 3;
new_prt.ResolutionOfTime = 'Volumes';
new_prt.NrOfConditions = 4;

all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*difficulty(i-1)] ]; end
new_prt.AddCond('difficulty', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*emotion(i-1)] ]; end
new_prt.AddCond('emotion', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*familiarity(i-1)] ]; end
new_prt.AddCond('familiarity', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_looking(i-1)] ]; end
new_prt.AddCond('strategy_looking', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_map(i-1)] ]; end
new_prt.AddCond('strategy_map', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_moving(i-1)] ]; end
new_prt.AddCond('strategy_moving', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_lines(i-1)] ]; end
new_prt.AddCond('strategy_lines', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_feeling(i-1)] ]; end
new_prt.AddCond('strategy_feeling', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_1pp(i-1)] ]; end
new_prt.AddCond('strategy_1pp', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_3pp(i-1)] ]; end
new_prt.AddCond('strategy_3pp', all_onoffsets);
% adding RT
all_onoffsets = []; for i=1:length(RTs), if RTs(i)~=0, all_onoffsets = [all_onoffsets; [i i RTs(i)]]; end, end
new_prt.AddCond('RT', all_onoffsets);
% new_prt.SaveAs([output_file(1:end-4) '.prt']);


% Create a new SDM with only the confound predictors
sdm_confounds = new_prt.CreateSDM(struct('nvol', 200, 'prtr', 2500, 'rcond', 0));
n = 2;
if ~strcmp(sdm_confounds.PredictorNames{n},'difficulty x p1')    % no difficulty ranking for this subject
    disp('no difficulty ranking for this subject!');
    sdm_confounds.PredictorNames = [sdm_confounds.PredictorNames(1:n-1) 'difficulty x p1' sdm_confounds.PredictorNames(n:end)];
    sdm_confounds.PredictorColors = [sdm_confounds.PredictorColors(1:n-1,:); sdm_confounds.PredictorColors(1,:); sdm_confounds.PredictorColors(n:end,:)];
    sdm_confounds.SDMMatrix = [sdm_confounds.SDMMatrix(:,1:n-1) zeros(size(sdm_confounds.SDMMatrix,1),1) sdm_confounds.SDMMatrix(:,n:end)];
    sdm_confounds.RTCMatrix = [sdm_confounds.RTCMatrix(:,1:n-1) zeros(size(sdm_confounds.RTCMatrix,1),1) sdm_confounds.RTCMatrix(:,n:end)];
end
n = n + 2;
if ~strcmp(sdm_confounds.PredictorNames{n},'emotion x p1')    % no difficulty ranking for this subject
    disp('no emotion ranking for this subject!');
    sdm_confounds.PredictorNames = [sdm_confounds.PredictorNames(1:n-1) 'emotion x p1' sdm_confounds.PredictorNames(n:end)];
    sdm_confounds.PredictorColors = [sdm_confounds.PredictorColors(1:n-1,:); sdm_confounds.PredictorColors(1,:); sdm_confounds.PredictorColors(n:end,:)];
    sdm_confounds.SDMMatrix = [sdm_confounds.SDMMatrix(:,1:n-1) zeros(size(sdm_confounds.SDMMatrix,1),1) sdm_confounds.SDMMatrix(:,n:end)];
    sdm_confounds.RTCMatrix = [sdm_confounds.RTCMatrix(:,1:n-1) zeros(size(sdm_confounds.RTCMatrix,1),1) sdm_confounds.RTCMatrix(:,n:end)];
end
n = n + 2;
if ~strcmp(sdm_confounds.PredictorNames{n},'familiarity x p1')    % no difficulty ranking for this subject
    disp('no familiarity ranking for this subject!');
    sdm_confounds.PredictorNames = [sdm_confounds.PredictorNames(1:n-1) 'familiarity x p1' sdm_confounds.PredictorNames(n:end)];
    sdm_confounds.PredictorColors = [sdm_confounds.PredictorColors(1:n-1,:); sdm_confounds.PredictorColors(1,:); sdm_confounds.PredictorColors(n:end,:)];
    sdm_confounds.SDMMatrix = [sdm_confounds.SDMMatrix(:,1:n-1) zeros(size(sdm_confounds.SDMMatrix,1),1) sdm_confounds.SDMMatrix(:,n:end)];
    sdm_confounds.RTCMatrix = [sdm_confounds.RTCMatrix(:,1:n-1) zeros(size(sdm_confounds.RTCMatrix,1),1) sdm_confounds.RTCMatrix(:,n:end)];
end
n = n + 2;
if ~strcmp(sdm_confounds.PredictorNames{n},'strategy_looking x p1')    % no difficulty ranking for this subject
    disp('no strategy_looking ranking for this subject!');
    sdm_confounds.PredictorNames = [sdm_confounds.PredictorNames(1:n-1) 'strategy_looking x p1' sdm_confounds.PredictorNames(n:end)];
    sdm_confounds.PredictorColors = [sdm_confounds.PredictorColors(1:n-1,:); sdm_confounds.PredictorColors(1,:); sdm_confounds.PredictorColors(n:end,:)];
    sdm_confounds.SDMMatrix = [sdm_confounds.SDMMatrix(:,1:n-1) zeros(size(sdm_confounds.SDMMatrix,1),1) sdm_confounds.SDMMatrix(:,n:end)];
    sdm_confounds.RTCMatrix = [sdm_confounds.RTCMatrix(:,1:n-1) zeros(size(sdm_confounds.RTCMatrix,1),1) sdm_confounds.RTCMatrix(:,n:end)];
end
n = n + 2;
if ~strcmp(sdm_confounds.PredictorNames{n},'strategy_map x p1')    % no difficulty ranking for this subject
    disp('no strategy_map ranking for this subject!');
    sdm_confounds.PredictorNames = [sdm_confounds.PredictorNames(1:n-1) 'strategy_map x p1' sdm_confounds.PredictorNames(n:end)];
    sdm_confounds.PredictorColors = [sdm_confounds.PredictorColors(1:n-1,:); sdm_confounds.PredictorColors(1,:); sdm_confounds.PredictorColors(n:end,:)];
    sdm_confounds.SDMMatrix = [sdm_confounds.SDMMatrix(:,1:n-1) zeros(size(sdm_confounds.SDMMatrix,1),1) sdm_confounds.SDMMatrix(:,n:end)];
    sdm_confounds.RTCMatrix = [sdm_confounds.RTCMatrix(:,1:n-1) zeros(size(sdm_confounds.RTCMatrix,1),1) sdm_confounds.RTCMatrix(:,n:end)];
end
n = n + 2;
if ~strcmp(sdm_confounds.PredictorNames{n},'strategy_moving x p1')    % no difficulty ranking for this subject
    disp('no strategy_moving ranking for this subject!');
    sdm_confounds.PredictorNames = [sdm_confounds.PredictorNames(1:n-1) 'strategy_moving x p1' sdm_confounds.PredictorNames(n:end)];
    sdm_confounds.PredictorColors = [sdm_confounds.PredictorColors(1:n-1,:); sdm_confounds.PredictorColors(1,:); sdm_confounds.PredictorColors(n:end,:)];
    sdm_confounds.SDMMatrix = [sdm_confounds.SDMMatrix(:,1:n-1) zeros(size(sdm_confounds.SDMMatrix,1),1) sdm_confounds.SDMMatrix(:,n:end)];
    sdm_confounds.RTCMatrix = [sdm_confounds.RTCMatrix(:,1:n-1) zeros(size(sdm_confounds.RTCMatrix,1),1) sdm_confounds.RTCMatrix(:,n:end)];
end
n = n + 2;
if ~strcmp(sdm_confounds.PredictorNames{n},'strategy_lines x p1')    % no difficulty ranking for this subject
    disp('no strategy_lines ranking for this subject!');
    sdm_confounds.PredictorNames = [sdm_confounds.PredictorNames(1:n-1) 'strategy_lines x p1' sdm_confounds.PredictorNames(n:end)];
    sdm_confounds.PredictorColors = [sdm_confounds.PredictorColors(1:n-1,:); sdm_confounds.PredictorColors(1,:); sdm_confounds.PredictorColors(n:end,:)];
    sdm_confounds.SDMMatrix = [sdm_confounds.SDMMatrix(:,1:n-1) zeros(size(sdm_confounds.SDMMatrix,1),1) sdm_confounds.SDMMatrix(:,n:end)];
    sdm_confounds.RTCMatrix = [sdm_confounds.RTCMatrix(:,1:n-1) zeros(size(sdm_confounds.RTCMatrix,1),1) sdm_confounds.RTCMatrix(:,n:end)];
end
n = n + 2;
if ~strcmp(sdm_confounds.PredictorNames{n},'strategy_feeling x p1')    % no difficulty ranking for this subject
    disp('no strategy_feeling ranking for this subject!');
    sdm_confounds.PredictorNames = [sdm_confounds.PredictorNames(1:n-1) 'strategy_feeling x p1' sdm_confounds.PredictorNames(n:end)];
    sdm_confounds.PredictorColors = [sdm_confounds.PredictorColors(1:n-1,:); sdm_confounds.PredictorColors(1,:); sdm_confounds.PredictorColors(n:end,:)];
    sdm_confounds.SDMMatrix = [sdm_confounds.SDMMatrix(:,1:n-1) zeros(size(sdm_confounds.SDMMatrix,1),1) sdm_confounds.SDMMatrix(:,n:end)];
    sdm_confounds.RTCMatrix = [sdm_confounds.RTCMatrix(:,1:n-1) zeros(size(sdm_confounds.RTCMatrix,1),1) sdm_confounds.RTCMatrix(:,n:end)];
end
n = n + 2;
if ~strcmp(sdm_confounds.PredictorNames{n},'strategy_1pp x p1')    % no difficulty ranking for this subject
    disp('no strategy_1pp ranking for this subject!');
    sdm_confounds.PredictorNames = [sdm_confounds.PredictorNames(1:n-1) 'strategy_1pp x p1' sdm_confounds.PredictorNames(n:end)];
    sdm_confounds.PredictorColors = [sdm_confounds.PredictorColors(1:n-1,:); sdm_confounds.PredictorColors(1,:); sdm_confounds.PredictorColors(n:end,:)];
    sdm_confounds.SDMMatrix = [sdm_confounds.SDMMatrix(:,1:n-1) zeros(size(sdm_confounds.SDMMatrix,1),1) sdm_confounds.SDMMatrix(:,n:end)];
    sdm_confounds.RTCMatrix = [sdm_confounds.RTCMatrix(:,1:n-1) zeros(size(sdm_confounds.RTCMatrix,1),1) sdm_confounds.RTCMatrix(:,n:end)];
end
n = n + 2;
if ~strcmp(sdm_confounds.PredictorNames{n},'strategy_3pp x p1')    % no difficulty ranking for this subject
    disp('no strategy_3pp ranking for this subject!');
    sdm_confounds.PredictorNames = [sdm_confounds.PredictorNames(1:n-1) 'strategy_3pp x p1' sdm_confounds.PredictorNames(n:end)];
    sdm_confounds.PredictorColors = [sdm_confounds.PredictorColors(1:n-1,:); sdm_confounds.PredictorColors(1,:); sdm_confounds.PredictorColors(n:end,:)];
    sdm_confounds.SDMMatrix = [sdm_confounds.SDMMatrix(:,1:n-1) zeros(size(sdm_confounds.SDMMatrix,1),1) sdm_confounds.SDMMatrix(:,n:end)];
    sdm_confounds.RTCMatrix = [sdm_confounds.RTCMatrix(:,1:n-1) zeros(size(sdm_confounds.RTCMatrix,1),1) sdm_confounds.RTCMatrix(:,n:end)];
end


% using only the needed predictors - removing the ones without parametric modulation
sdm_confounds.NrOfPredictors = 11;
sdm_confounds.FirstConfoundPredictor = 5;
predictor_general = sdm_confounds.SDMMatrix(:, 1);      % Saving the predictor without parametric modulation
sdm_confounds.PredictorColors = sdm_confounds.PredictorColors([2 4 6 8 10 12 14 16 18 20 22],:);
sdm_confounds.PredictorNames = sdm_confounds.PredictorNames([2 4 6 8 10 12 14 16 18 20 22]);
sdm_confounds.SDMMatrix = sdm_confounds.SDMMatrix(:,[2 4 6 8 10 12 14 16 18 20 22]);
sdm_confounds.RTCMatrix = sdm_confounds.RTCMatrix(:,[2 4 6 8 10 12 14 16 18 20 22]);

% read the original SDM and add the confounds
sdm_inc_confounds = xff(original_scales_sdm);
sdm_inc_confounds.NrOfPredictors = sdm_inc_confounds.NrOfPredictors + sdm_confounds.NrOfPredictors;
sdm_inc_confounds.PredictorColors = [sdm_inc_confounds.PredictorColors(1:6,:); sdm_confounds.PredictorColors; sdm_inc_confounds.PredictorColors(7:end,:)];
sdm_inc_confounds.PredictorNames = [sdm_inc_confounds.PredictorNames(1:6) sdm_confounds.PredictorNames sdm_inc_confounds.PredictorNames(7:end)];
sdm_inc_confounds.SDMMatrix = [sdm_inc_confounds.SDMMatrix(:,1:6) sdm_confounds.SDMMatrix sdm_inc_confounds.SDMMatrix(:,7:end)];
sdm_inc_confounds.FirstConfoundPredictor = 11;
sdm_inc_confounds.SaveAs(output_file_scales_confounds);

% Save each confound separately as an sdm, with the original motion predictors from the data
sdm_orig = xff(original_scales_sdm);
sdm_one_confound = xff(original_scales_sdm);        % An sdm with one parametrically modulated predictor
sdm_one_confound.NrOfPredictors = sdm_one_confound.NrOfPredictors - 6 + 2;
sdm_one_confound.FirstConfoundPredictor = 3;
sdm_scales_and_oneconfound = xff(original_scales_sdm);      % An sdm with the 6 scales and a parametrically modulated predictor
sdm_scales_and_oneconfound.NrOfPredictors = sdm_scales_and_oneconfound.NrOfPredictors + 1;
sdm_scales_and_oneconfound.FirstConfoundPredictor = 8;
for i = 1:size(sdm_confounds.SDMMatrix, 2)
    if sum(sdm_confounds.SDMMatrix(:,i)) ~= 0           % Verifying that the confound predictor exists for this subject
        % Creating an sdm with only this variable and motion predictors
        sdm_one_confound.PredictorColors = [sdm_confounds.PredictorColors(i,:); sdm_confounds.PredictorColors(i,:); sdm_orig.PredictorColors(7:end,:)];
        sdm_one_confound.PredictorNames = ['predictor_general' sdm_confounds.PredictorNames(i) sdm_orig.PredictorNames(7:end)];
        sdm_one_confound.SDMMatrix = [predictor_general sdm_confounds.SDMMatrix(:,i) sdm_orig.SDMMatrix(:,7:end)];
        sdm_one_confound.SaveAs([original_scales_sdm(1:end-4), '_oneconfound_', sdm_confounds.PredictorNames{i}(1:end-5), '.sdm']);
        % Creating an sdm with the 6 scales and this variable and motion predictors
        sdm_scales_and_oneconfound.PredictorColors = [sdm_orig.PredictorColors(1:6,:); sdm_confounds.PredictorColors(i,:); sdm_orig.PredictorColors(7:end,:)];
        sdm_scales_and_oneconfound.PredictorNames = [sdm_orig.PredictorNames(1:6) sdm_confounds.PredictorNames(i) sdm_orig.PredictorNames(7:end)];
        sdm_scales_and_oneconfound.SDMMatrix = [sdm_orig.SDMMatrix(:,1:6) sdm_confounds.SDMMatrix(:,i) sdm_orig.SDMMatrix(:,7:end)];
        sdm_scales_and_oneconfound.SaveAs([original_scales_sdm(1:end-4), '_scalesandoneconfound_', sdm_confounds.PredictorNames{i}(1:end-5), '.sdm']);        
    end
end

% saving the PRT file with the confounds but without z-scoring of RT, to be able to compute average RT per condition

% read RT from presentation log file
fid = fopen(logfile_for_RT);
text = textscan(fid,'%s');
text = text{1};
fclose(fid);
% find pulse indices
f = strfind(text,'Pulse'); indices_pulse = []; for i=1:length(f), if f{i}==1, indices_pulse = [indices_pulse i]; end,end
% find responses
RTs = nan(length(indices_pulse),1);
for i=1:length(indices_pulse)-1
    r = strfind(text(indices_pulse(i):indices_pulse(i+1)), 'Response');
    r = find(cellfun(@(x) ~isempty(x), r));
    if ~isempty(r)
        RTs(i) = str2num(text{indices_pulse(i) + r(1) + 2});    % if there are two responses, taking the first RT of the two
    end
end
% RTs(~isnan(RTs)) = zscore(RTs(~isnan(RTs)));
RTs(isnan(RTs)) = 0;
% read the original PRT (with no confounds, only times for each condition)
prt = xff(original_allconds_prt);
new_prt = xff('new:prt');
new_prt.FileVersion = 3;
new_prt.ResolutionOfTime = 'Volumes';
new_prt.NrOfConditions = 4;
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*difficulty(i-1)] ]; end
new_prt.AddCond('difficulty', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*emotion(i-1)] ]; end
new_prt.AddCond('emotion', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*familiarity(i-1)] ]; end
new_prt.AddCond('familiarity', all_onoffsets);
% adding RT
all_onoffsets = []; for i=1:length(RTs), if RTs(i)~=0, all_onoffsets = [all_onoffsets; [i i RTs(i)]]; end, end
new_prt.AddCond('RT', all_onoffsets);

all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_looking(i-1)] ]; end
new_prt.AddCond('strategy_looking', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_map(i-1)] ]; end
new_prt.AddCond('strategy_map', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_moving(i-1)] ]; end
new_prt.AddCond('strategy_moving', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_lines(i-1)] ]; end
new_prt.AddCond('strategy_lines', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_feeling(i-1)] ]; end
new_prt.AddCond('strategy_feeling', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_1pp(i-1)] ]; end
new_prt.AddCond('strategy_1pp', all_onoffsets);
all_onoffsets = []; for i=2:1+num_conditions, all_onoffsets = [all_onoffsets; [prt.Cond(i).OnOffsets ones(size(prt.Cond(i).OnOffsets,1),1)*strategy_3pp(i-1)] ]; end
new_prt.AddCond('strategy_3pp', all_onoffsets);

new_prt.SaveAs([output_file_scales_confounds(1:end-4) '.prt']);


