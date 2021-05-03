%% Test correlation between popularity in network and strength of relations
s = 1;
subj = subject_names(s).name;
disp(s)

% Defining subject-specific directories
current_subject_dir = fullfile(all_subjects_files_dir, subj);
current_subject_RDM_dir = fullfile(RDMs_output_files_dir, subj);

% Reading the subjects' excel file with the individuals and their organization into social groups (24 individuals, 6 in each group)
excel_names_path = dir(fullfile(current_subject_dir, '*friends_groups*.xlsx'));

% Getting the names of the individuals
[~, ~, people_names] = xlsread(fullfile(current_subject_dir, excel_names_path(1).name),['B3:E', num2str(2 + people_in_group_num)]);
people_names = reshape(people_names, 1, groups_num * people_in_group_num);
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
facebook_distance_RDM_percentcommonfriends = zeros(length(nodes));

for i=1:length(edges)
    for j = 1:length(edges)
        facebook_distance_RDM_percentcommonfriends(i,j) = length(intersect(find(edges(i,:)), find(edges(j,:))));
        facebook_distance_RDM_percentcommonfriends(i,j) = facebook_distance_RDM_percentcommonfriends(i,j) * 2 / (length(find(edges(i,:))) + length(find(edges(j,:))));
        %         facebook_distance_RDM_percentcommonfriends(i,j) = (facebook_distance_RDM_percentcommonfriends(i,j) / length(find(edges(i,:))) + facebook_distance_RDM_percentcommonfriends(i,j) / length(find(edges(j,:)))) / 2;
    end
end
for i = 1:num_people, facebook_distance_RDM_percentcommonfriends(i,i) = max(facebook_distance_RDM_percentcommonfriends(:));
end
facebook_distance_RDM_percentcommonfriends = facebook_distance_RDM_percentcommonfriends / max(facebook_distance_RDM_percentcommonfriends(:));
% Converting to distance (dissimilarity) matrix
facebook_distance_RDM_percentcommonfriends = 1 - facebook_distance_RDM_percentcommonfriends;
figure;imagesc(facebook_distance_RDM_percentcommonfriends)

% Calculating the relation between popularity and connectivity strength
m = mean(facebook_distance_RDM_percentcommonfriends);
e = sum(edges);
figure;scatter(m,e)


%% Compute relation between ratings of friendship strength and Facebook connectivity
ratings_folder = 'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Old\Test_pilot_corr_dist_ratings';
ratings_files = {'Michael_ratings.xlsx','Bar_ratings.xlsx','Mordechai_ratings.xlsx','Neta_ratings.xlsx'};
graphml_files = {'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Old\Data_pilot\Michael\Facebook_data\michael.peer_new.graphml',...
    'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Old\Data_pilot\Bar\Facebook_data\bar_tamir_nopics.graphml',...
    'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Old\Data_pilot\Mordechai\Facebook_data\mordechai.hayman.graphml',...
    'C:\Users\michaelpeer1\Google Drive\bar_mord_shared_files\Experiment_data\Graphml_files\180328HA_NE\neta.hayman.graphml'};

all_corr = zeros(1, length(ratings_files));
all_corr_groups = zeros(1, length(ratings_files));
for s=1:length(ratings_files)
    disp(s)
    
    % Reading graphml data
    graph = xml2struct(graphml_files{s});
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
    [~, ~, curr_xls_data] = xlsread(fullfile(ratings_folder, ratings_files{s}),'B2:Z25');
    people_names = curr_xls_data(:,1);
    curr_ratings = curr_xls_data(:,2:end); curr_ratings = cell2mat(curr_ratings);
    num_people = length(people_names);
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
    facebook_distance_RDM_percentcommonfriends = zeros(length(nodes));
    
    for i=1:length(edges)
        for j = 1:length(edges)
            facebook_distance_RDM_percentcommonfriends(i,j) = length(intersect(find(edges(i,:)), find(edges(j,:))));
            facebook_distance_RDM_percentcommonfriends(i,j) = facebook_distance_RDM_percentcommonfriends(i,j) * 2 / (length(find(edges(i,:))) + length(find(edges(j,:))));
            %         facebook_distance_RDM_percentcommonfriends(i,j) = (facebook_distance_RDM_percentcommonfriends(i,j) / length(find(edges(i,:))) + facebook_distance_RDM_percentcommonfriends(i,j) / length(find(edges(j,:)))) / 2;
        end
    end
    for i = 1:num_people, facebook_distance_RDM_percentcommonfriends(i,i) = max(facebook_distance_RDM_percentcommonfriends(:));
    end
    facebook_distance_RDM_percentcommonfriends = facebook_distance_RDM_percentcommonfriends / max(facebook_distance_RDM_percentcommonfriends(:));
    % Converting to distance (dissimilarity) matrix
    facebook_distance_RDM_percentcommonfriends = 1 - facebook_distance_RDM_percentcommonfriends;
    facebook_distance_RDM_percentcommonfriends = facebook_distance_RDM_percentcommonfriends(people_numbers, people_numbers);
    
    % Comparing ratings to Facebook data
    nondiag_locs = find(tril(ones(num_people),-1));
    all_corr(s) = corr(curr_ratings(nondiag_locs), facebook_distance_RDM_percentcommonfriends(nondiag_locs));
    % Computing within social groups only
    g = zeros(24); g(1:6,1:6)=1;g(7:12,7:12)=1;g(13:18,13:18)=1;g(19:24,19:24)=1;
    gg = find(tril(g,-1));
    all_corr_groups(s) = corr(curr_ratings(gg), facebook_distance_RDM_percentcommonfriends(gg));
end
