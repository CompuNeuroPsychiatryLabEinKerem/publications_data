function socialnetworks_convert_logfile_to_prt(logfile_to_convert, output_prt_filename, template_prt_file)
% A function for conversion of presentation logfiles to BrainVoyager protocol (PRT) files


% Open and read the log file
fid = fopen(logfile_to_convert);
text = textscan(fid,'%s');
text = text{1};
fclose(fid);


% Find the MRI pulses in the log file
f = strfind(text,'Pulse');
indices_pulse = []; 
for i = 1:length(f)
    if f{i} == 1
        indices_pulse = [indices_pulse i]; 
    end
end


% Fill in the conditions
conditions = {};
for i = 1:length(indices_pulse)
    if indices_pulse(i) + 7 < length(text)            % Excluding the last pulse
        if length(text{indices_pulse(i)+7}) >= 2    % Ignoring user responses (user responses are single digits)
            conditions{i} = text{indices_pulse(i) + 7};
        else                                        % If it is a user response, skip to next thing happening
            conditions{i} = text{indices_pulse(i) + 13};
        end
    end
end
conditions{length(indices_pulse)} = '30';   % fill in the last condition


% Change responses to rest
locations_responses_in_conditions = find(cellfun(@length, conditions) == 1);
if length(locations_responses_in_conditions) == 1
    conditions{locations_responses_in_conditions} = '30';
elseif length(locations_responses_in_conditions) > 1
    for i=1:length(locations_responses_in_conditions), conditions{locations_responses_in_conditions(i)} = '30'; end
end


% Find when a 'people' condition occured
conds_all_times = find(cellfun(@(x) strcmp(x(1:2),'pe'),conditions));


% Get the stimulus number (from 1 to 24) (e.g. from 'people_22')
conds_stimuli_nums = cellfun(@(x) x(8:end), conditions(conds_all_times),'UniformOutput',0);
conds = cellfun(@str2num, conds_stimuli_nums);      % Convert to numbers 


% Record the condition number instead of each stimulus
conditions_new = conditions;
for i=1:length(conds_all_times)
    conditions_new{conds_all_times(i)} = num2str(conds(i));
end


% Convert all empty stimuli into the previous ones
for i=1:length(conditions_new)
    if strcmp(conditions_new{i}, '30')
        conditions_new{i} = conditions_new{i-1};
    end
end


% Convert all rest types to 'rest' label
for i=1:length(conditions_new)
    if length(conditions_new{i}) > 2
        if strcmp(conditions_new{i}(end-3:end), 'rest')
            conditions_new{i} = 'rest';
        end
    end
end


% Create PRT file according to template
current_prt = xff(template_prt_file);

conditions_onoffsets = cell(current_prt.NrOfConditions,1);
i=1;
conditions_new{end+1} = 0;  % add to avoid running over edge of array in the while loop
while i<length(conditions_new)-1
    if strcmp(conditions_new{i},'rest')
        condition_onset_current = i;
        while strcmp(conditions_new{i+1},'rest')
            i = i+1;
        end
        condition_offset_current = i;
        conditions_onoffsets{1} = [conditions_onoffsets{1} ; condition_onset_current condition_offset_current];
        
    elseif length(conditions_new{i})>10
        if strcmp(conditions_new{i}(1:12),'instructions')
            condition_onset_current = i;
            while length(conditions_new{i+1})>10
                if strcmp(conditions_new{i+1}(1:12),'instructions')
                    i = i+1;
                end
            end
            condition_offset_current = i;
            conditions_onoffsets{current_prt.NrOfConditions} = [conditions_onoffsets{current_prt.NrOfConditions} ; condition_onset_current condition_offset_current];
        end
        
    else
        current_condition = str2num(conditions_new{i});
        condition_onset_current = i;
        condition_offset_current = i+1;
        conditions_onoffsets{current_condition+1} = [conditions_onoffsets{current_condition+1}; condition_onset_current condition_offset_current];
        i = i+1;
    end
    
    i = i+1;
end


for c=1:current_prt.NrOfConditions
    current_prt.Cond(c).OnOffsets = conditions_onoffsets{c};
    current_prt.Cond(c).NrOfOnOffsets = size(conditions_onoffsets{c},1);
end

% Saving the resulting PRT file
current_prt.SaveAs(output_prt_filename);

