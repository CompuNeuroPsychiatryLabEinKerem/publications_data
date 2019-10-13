%% Creating dissimilarity matrices for each subject

% Defining working directories
all_subjects_files_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\Graphml_files';          % Contains GraphML format files of subjects' social network, extracted using the LostCircles tool (https://lostcircles.com), and the subjects' provided excel file with the stimuli names (organized in a separate folder for each subject)
RDMs_output_files_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\Dissimilarity_mats';      % A directory to store all of the dissimilarity matrices created here (organized in a separate folder for each subject)
experiment_logfiles_dir = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\Experiment_stimuli';    % A directory with the log files of the experiment (organized in a separate folder for each subject)
% Finding the subjects' names (directories)
subject_names = dir(all_subjects_files_dir); subject_names = subject_names(3:end); subject_names = subject_names([subject_names.isdir]);

% Defining general analysis parameters
groups_num = 4;             % Numbber of social groups provided by the subjects
people_in_group_num = 6;    % Number of individuals (names) in each group provided by the subjects
num_people = groups_num * people_in_group_num;      % Total number of individuals names (stimuli)


% Creating the dissimilarity matrices - looping over subjects
for s = 1:length(subject_names)
    subj = subject_names(s).name;
    disp(s)
    
    % Defining subject-specific directories
    current_subject_dir = fullfile(all_subjects_files_dir, subj);
    current_subject_RDM_dir = fullfile(RDMs_output_files_dir, subj);
    
    
    %% Getting the subject's Facebook social network data from her/his downloaded graphml file, which was created using the Lost Circles chrome plugin
    % Getting the 24 individuals that the subject saw during the experiment:
    % Reading the subjects' excel file with the individuals and their organization into social groups (24 individuals, 6 in each group)
    excel_names_path = dir(fullfile(current_subject_dir, '*friends_groups*.xlsx'));
    
    % Getting the names of the individuals
    [~, ~, people_names] = xlsread(fullfile(current_subject_dir, excel_names_path(1).name),['B3:E', num2str(2 + people_in_group_num)]);
    people_names = reshape(people_names, 1, groups_num * people_in_group_num);
    
    % Read the grapml file that contains subjects' social network connections
    graphml_filename = dir(fullfile(current_subject_dir, '*.graphml'));
    graph = xml2struct(fullfile(current_subject_dir, graphml_filename(1).name));
    graph = graph.Children(7);
    
    % Get all of the nodes (individuals) and their connections (edges) from the graph
    nodes = {}; edges = [];
    for i = 1:length(graph.Children)
        if strcmp(graph.Children(i).Name,'node')
            nodes{end+1} = graph.Children(i).Children(2).Children.Data;
        elseif strcmp(graph.Children(i).Name,'edge')
            node1_current = str2num(graph.Children(i).Attributes(2).Value(2:end));
            node2_current = str2num(graph.Children(i).Attributes(3).Value(2:end));
            edges(node1_current+1, node2_current+1) = 1;    % adding one since the nodes are numbered from zero
            edges(node2_current+1, node1_current+1) = 1;    % adding one since the nodes are numbered from zero
        end
    end
    
    % get the node number for each of the 24 individuals provided by the subject
    people_names = lower(people_names); nodes = lower(nodes);    % convert all names to lower case to be able to compare them
    people_numbers = zeros(1, num_people);
    for i = 1:num_people
        num_current_person = find(ismember(nodes, people_names(i)));
        if length(num_current_person) == 1
            people_numbers(i) = num_current_person;
        else    % Ignoring individuals whose name has no Facebook correlate, or more than one Facebook correlates
            people_numbers(i) = nan;
        end
    end
    
    
    %% Creating a social network distance dissimilarity matrix - distances calculated according to the proportion of mutual friends out of all friends (embeddedness measure, Granovetter 1973, Bapna & Sundararajan 2017)
    facebook_distance_RDM_percentcommonfriends = zeros(num_people);
    for i = 1:num_people
        for j = 1:num_people
            if ~isnan(people_numbers(i)) && ~isnan(people_numbers(j))
                % Finding the number of mutual friends
                facebook_distance_RDM_percentcommonfriends(i,j) = length(intersect(find(edges(people_numbers(i),:)), find(edges(people_numbers(j),:))));
                % Adding the direct connection between the two individuals, if it exists - as only second-level connections (friends of friends) are counted in the previous stage
                if sum(find(edges(people_numbers(i),:)) == people_numbers(j)) ~= 0
                    facebook_distance_RDM_percentcommonfriends(i,j) = facebook_distance_RDM_percentcommonfriends(i,j) + 1;
                end
                % Calculating the proportion of mutual friends out of all friends of both individuals
                facebook_distance_RDM_percentcommonfriends(i,j) = facebook_distance_RDM_percentcommonfriends(i,j) * 2 / (length(find(edges(people_numbers(i),:))) + length(find(edges(people_numbers(j),:))));
            else    % If one of the names does not exist in the network, put NaN in the matrix at this location
                facebook_distance_RDM_percentcommonfriends(i,j) = nan;
            end
        end
    end
    % Removing values along the diagonal
    for i = 1:num_people, facebook_distance_RDM_percentcommonfriends(i,i) = max(facebook_distance_RDM_percentcommonfriends(:)); end
    % Normalizing the matrix values to range 0-1
    facebook_distance_RDM_percentcommonfriends = facebook_distance_RDM_percentcommonfriends / max(facebook_distance_RDM_percentcommonfriends(:));
    % Converting to distance (dissimilarity) matrix
    facebook_distance_RDM_percentcommonfriends = 1 - facebook_distance_RDM_percentcommonfriends;
    % Saving the RDM
    facebook_distance_RDM = facebook_distance_RDM_percentcommonfriends;
    save(fullfile(current_subject_RDM_dir, 'facebook_distance_percentcommonfriends_RDM.mat'), 'facebook_distance_RDM');
    
    
    
    %% Creating social network distance dissimilarity matrices according to other distance measures
    % Distance measure - total number of mutual friends (1st+2nd level connections between each pair of individuals)
    facebook_distance_RDM_numcommonfriends = zeros(num_people);
    for i = 1:num_people
        for j = 1:num_people
            if ~isnan(people_numbers(i)) && ~isnan(people_numbers(j))
                % Finding the number of mutual friends
                facebook_distance_RDM_numcommonfriends(i,j) = length(intersect(find(edges(people_numbers(i),:)), find(edges(people_numbers(j),:))));
                % Adding the direct connection between the two individuals, if it exists - as only second-level connections (friends of friends) are counted in the previous stage
                if sum(find(edges(people_numbers(i),:)) == people_numbers(j)) ~= 0
                    facebook_distance_RDM_numcommonfriends(i,j) = facebook_distance_RDM_numcommonfriends(i,j) + 1;
                end
            else    % If one of the names does not exist in the network, put NaN in the matrix at this location
                facebook_distance_RDM_numcommonfriends(i,j) = nan;
            end
        end
    end
    % Removing values along the diagonal
    for i = 1:num_people, facebook_distance_RDM_numcommonfriends(i,i) = max(facebook_distance_RDM_numcommonfriends(:)); end
    % Normalizing the matrix values to range 0-1
    facebook_distance_RDM_numcommonfriends = facebook_distance_RDM_numcommonfriends / max(facebook_distance_RDM_numcommonfriends(:));
    % Converting to distance (dissimilarity) matrix
    facebook_distance_RDM_numcommonfriends = 1 - facebook_distance_RDM_numcommonfriends;
    % Saving the RDM
    facebook_distance_RDM = facebook_distance_RDM_numcommonfriends;
    save(fullfile(current_subject_RDM_dir, 'facebook_distance_numcommonfriends_RDM.mat'), 'facebook_distance_RDM');
    
    
    % Distance measure - existence of a direct connection between individuals (binary measure)
    facebook_distance_RDM_directconnections = zeros(num_people);
    for i = 1:num_people
        for j = 1:num_people
            if ~isnan(people_numbers(i)) && ~isnan(people_numbers(j))
                % Looking if there is an edge between the individuals
                facebook_distance_RDM_directconnections(i,j) = edges(people_numbers(i), people_numbers(j));
            else    % If one of the names does not exist in the network, put NaN in the matrix at this location
                facebook_distance_RDM_directconnections(i,j) = nan;
            end
        end
    end
    % Normalizing the matrix values to range 0-1
    for i = 1:num_people, facebook_distance_RDM_directconnections(i,i) = max(facebook_distance_RDM_directconnections(:)); end
    % Converting to distance (dissimilarity) matrix
    facebook_distance_RDM_directconnections = 1 - facebook_distance_RDM_directconnections;
    % Saving the RDM
    facebook_distance_RDM = facebook_distance_RDM_directconnections;
    save(fullfile(current_subject_RDM_dir, 'facebook_distance_RDM_directconnections.mat'), 'facebook_distance_RDM');

    
    % Distance measure - shortest path between individuals along the network graph
    facebook_distance_RDM_shortestpath = zeros(num_people);
    % Calculating the shortest path between all individuals, using an algorithm from the Brain Connectivity Toolbox (Rubinov & Sporns, 2010)
    shortest_path_full_mat = distance_bin(edges);
    for i = 1:num_people
        for j = 1:num_people
            if ~isnan(people_numbers(i)) && ~isnan(people_numbers(j))
                % Getting the shortest path between the individuals from the matrix created above
                facebook_distance_RDM_shortestpath(i,j) = shortest_path_full_mat(people_numbers(i), people_numbers(j));
            else    % If one of the names does not exist in the network, put NaN in the matrix at this location
                facebook_distance_RDM_shortestpath(i,j) = nan;
            end
        end
    end
    % Normalizing the matrix values to range 0-1
    facebook_distance_RDM_shortestpath = facebook_distance_RDM_shortestpath / max(facebook_distance_RDM_shortestpath(:));
    % Saving the RDM
    facebook_distance_RDM = facebook_distance_RDM_shortestpath;
    save(fullfile(current_subject_RDM_dir, 'facebook_distance_RDM_shortestpath.mat'), 'facebook_distance_RDM');
    
    
    % Distance measure -  communicability (graph theory measure, Estrada et al. 2008)
    facebook_distance_RDM_communicability = zeros(num_people);
    % Calculating the communicability between all graph nodes (individuals in the network), according to Estrada et al. 2008
    communicability_full_mat = expm(edges);
    for i = 1:num_people
        for j = 1:num_people
            if ~isnan(people_numbers(i)) && ~isnan(people_numbers(j))
                % Getting the communicability between the individuals from the matrix created above
                facebook_distance_RDM_communicability(i,j) = communicability_full_mat(people_numbers(i), people_numbers(j));
            else    % If one of the names does not exist in the network, put NaN in the matrix at this location
                facebook_distance_RDM_communicability(i,j) = nan;
            end
        end
    end
    % Removing values along the diagonal
    for i = 1:num_people, facebook_distance_RDM_communicability(i,i) = max(facebook_distance_RDM_communicability(:)); end
    % Normalizing the matrix values to range 0-1
    facebook_distance_RDM_communicability = facebook_distance_RDM_communicability / max(facebook_distance_RDM_communicability(:));
    % Converting to distance (dissimilarity) matrix
    facebook_distance_RDM_communicability = 1 - facebook_distance_RDM_communicability;
    % Saving the RDM
    facebook_distance_RDM = facebook_distance_RDM_communicability;
    save(fullfile(current_subject_RDM_dir, 'facebook_distance_RDM_communicability.mat'), 'facebook_distance_RDM');
    
    
    
    %% Creating dissimilarity matices by subjects' responses to the questions on personal affiliation, personality and appearance
    % Getting the subjects' responses to the questions from the experiment log files
    current_subject_logfiles_dir = fullfile(experiment_logfiles_dir, subj);
    current_logfiles = dir(fullfile(current_subject_logfiles_dir, '*.log'));
    for f = 1:length(current_logfiles)      % Going over each logfile and reading the subject's responses
        % Read the current logfile
        logfile = fullfile(current_logfiles(f).folder, current_logfiles(f).name);
        fid = fopen(logfile);
        text = textscan(fid,'%s');
        text = text{1};
        fclose(fid);
        
        % Identifying the questions asked in the current logfile - two questions per logfile (run)
        insructions_idices = find(strncmp(text, 'instructions_', 13) == 1);     % Locations in the logfile where a question was displayed
        quest1 = text(insructions_idices(1)); quest1 = str2double(quest1{:}(14:end));
        quest2 = text(insructions_idices(3)); quest2 = str2double(quest2{:}(14:end));
        
        % Identifying the logfile locations where individuals names and subject responses are recorded
        people_idices = find(strncmp(text, 'people_', 7) == 1);
        Response_idices = find(strncmp(text, 'Response', 8) == 1);
        Response_vec = str2double(text(Response_idices+1));
        
        % if a response is missing - replace it with '0' (to be later replaced to NaN)
        if length(people_idices) > length(Response_idices)
            fixed_Response_vec = zeros(length(people_idices),1);
            indices = [people_idices;numel(text)];
            for i3 = 1:length(indices)-1
                for i5 = indices(i3):indices(i3+1)-1
                    if strncmp(text{i5}, 'Response', 8)
                        fixed_Response_vec(i3) = str2double(text{i5+1});
                        break
                    end
                end
            end
            Response_vec = fixed_Response_vec;
        end
        quest_vec = [ones(length(people_idices)/2, 1) * quest1; ones(length(people_idices)/2, 1) * quest2];
        % Build the matrix
        for i2 = 1:length(people_idices)
            people = text(people_idices(i2));
            people = str2double(people{:}(8:end));
            all_responses_mat(people,quest_vec(i2)) = Response_vec(i2);
        end
    end
    % Replace all places with 0 in the response with NaN
    all_responses_mat(all_responses_mat == 0) = nan;
    
    % Computing dissimilarity according to different types of questions, using Euclidean distance between answers, and ignoring NaN values (unanswered questions)
    responses_self_proximity_RDM = zeros(24); responses_personality_RDM = zeros(24); responses_appearance_RDM = zeros(24); responses_all_RDM = zeros(24); 
    for i = 1:num_people
        for j = 1:num_people
            responses_self_proximity_RDM(i,j) = sqrt(nansum((all_responses_mat(i,1:4) - all_responses_mat(j,1:4)).^2));     % Personal affiliation questions
            responses_personality_RDM(i,j) = sqrt(nansum((all_responses_mat(i,5:8) - all_responses_mat(j,5:8)).^2));        % Personality traits questions
            responses_appearance_RDM(i,j) = sqrt(nansum((all_responses_mat(i,9:12) - all_responses_mat(j,9:12)).^2));       % Appearance questions
            responses_all_RDM(i,j) = sqrt(nansum((all_responses_mat(i,:) - all_responses_mat(j,:)).^2));                    % Overall response similarity
        end
    end
    
    % Saving the resulting RDMs
    save(fullfile(current_subject_RDM_dir, 'subject_responses_RDMs.mat'), 'responses_self_proximity_RDM', 'responses_personality_RDM', 'responses_appearance_RDM', 'responses_all_RDM');
end

