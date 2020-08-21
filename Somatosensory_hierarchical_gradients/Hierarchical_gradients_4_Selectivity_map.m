% Noam Saadon-Grosman
% August 2020
% This code creates selectivity surface map from the output of "Selectivity_computation"
%%% input: some files are in BrainVoyager format, for matrices extracted from
%%% other formats exclude the xff command. Fill in names and directories
% Template BrainVoyager map
% Somatosensory_response_tMap
% file with POIs (patch of intrest)
% Selectivity data from "Selectivity_computation"- .mat file
%%% output: BrainVoyager map= selectivity averaged across subjects and stimulus directions (lip and toe);
clear all;
hemi='RH';
threshold=0.01/8;
poidir='C:\Users\OWNER';% POIs directory
mapdir='C:\Users\OWNER';% template map directory
output='C:\Users\OWNER';% output directory
% load data and a template map
if isempty(strfind(hemi,'LH'))
    poiName='POIs_RH.poi';
    mapName='Template_Map_RH_fsaverage.smp';
    MaskDir='C:\Users\OWNER'; % somatosensory_response_tMap directory RH;
    %PA=22; %31pv area that has only one significant vertex and will be removed
else
    poiName='POIs_LH.poi';
    mapName='Template_Map_LH_fsaverage.smp';
    MaskDir='C:\Users\OWNER'; % somatosensory_response_tMap directory LH;
    %PA=5; % PEF area that has only one significant vertex and will be removed
end
Group_Mask=['Somatosensory_response_tMap_' hemi '.smp']; % somatosensory_response_tMap name RH;
cd('C:\Users\OWNER') % directory of selectivity data
load(['selectivity_lip_and_toe_' hemi '.mat'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(poidir)
POIs=xff(poiName);
POIs_data=POIs.POI;

% load somatosensory_response_tMap
cd(MaskDir);
mask=xff(Group_Mask);
maskData=mask.Map.SMPData;
for p=1:length(POIs_data)
    to_keep=[];
    for i=1:length(POIs_data(p).Vertices)
        if maskData(POIs_data(p).Vertices(i))~=0
            to_keep=[to_keep,i];
            MaskedV(p).masked=to_keep;
        end
    end
end
% loop on all POIs and average across subjects
All=dataval;
for p=1:length(All)
    if p==PA;
        All(p).widthAsubject=0;
    else
        for i=1:size(All(p).width,1) % loop on all vertices in POI
            temp=All(p).width(i,:);
            tempR=All(p).Pvalue(MaskedV(p).masked(i),:); % p vales
            above=tempR<threshold;
            tempthresh=temp(above); %apply significance threshold
            tempNZ=tempthresh(tempthresh~=0);
            tempNZnan=tempNZ(~isnan(tempNZ));
            temp_selectivity=sqrt(2)./tempNZnan; % selectivity as defined
            All(p).widthAsubject(i)=mean(temp_selectivity); % average across subjects and protocols
        end
    end
end


% save map
cd(mapdir)
smp=xff(mapName);
smp_data=smp.Map.SMPData;
new_smp=smp;
new_smp.Map.SMPData=zeros(size(smp_data));
% assign values within the map by POIs
for p=1:length(POIs_data)
    new_smp.Map.SMPData(POIs_data(p).Vertices(MaskedV(p).masked))=All(p).widthAsubject;
end

cd(output)
saveName=[hemi '_selectivity_map'];
new_smp.Map.Name=saveName;
new_smp.saveas([saveName '.smp']);
