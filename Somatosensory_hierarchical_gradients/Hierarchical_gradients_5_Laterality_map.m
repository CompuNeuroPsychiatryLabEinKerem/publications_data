% Noam Saadon-Grosman
% August 2020
% This code creates laterality surface map from the output of "Laterality_computation"
%%% input: some files are in BrainVoyager format, for matrices extracted from
%%% other formats exclude the xff command. Fill in names and directories
% Template BrainVoyager map
% Somatosensory_response_tMap
% file with POIs (patch of intrest)
% Laterality data from "Laterality_computation"- .mat file
%%% output: BrainVoyager map= laterality averaged across subjects and stimulus directions (lip and toe);
clear all;
close all;
poidir='C:\Users\OWNER';% POIs directory
mapdir='C:\Users\OWNER';% template map directory
output='C:\Users\OWNER';% output directory
hemi='LH';
threshold=0.01/16;
cd('C:\Users\OWNER') % directory of laterality data
load(['Laterality_lip_and_toe_' hemi '.mat'])
% loop on all POIs and average across subjects
dataval=datavalMs;
for p=1:length(dataval)
    for i=1:length(dataval(p).Laterality) % loop on all vertices in POI
        temp=dataval(p).Laterality(i,:);
        tempR=min(min(dataval(p).pValuesContra(i,:,:)),min(dataval(p).pValuesIpsi(i,:,:))); 
        tempR=squeeze(tempR(1,1,:));
        tempR=tempR';
        tempthresh=temp(tempR<threshold); %apply significance threshold
        tempNZ=tempthresh(tempthresh~=0); 
        tempNZnan=tempNZ(~isnan(tempNZ));
        dataval(p).LateralityAsubject(i)=mean(tempNZnan); % average across subject and stimulus direction (lip and toe)
        
    end
end

% load POI data and template map
if isempty(strfind(hemi,'LH'))
    poiName='POIs_RH.poi';
    mapName='Template_Map_RH_fsaverage.smp';
    MaskDir='C:\Users\OWNER'; % somatosensory_response_tMap directory RH;
else
    poiName='POIs_LH.poi';
    mapName='Template_Map_LH_fsaverage.smp';
    MaskDir='C:\Users\OWNER'; % somatosensory_response_tMap directory LH;
end
Group_Mask=['Somatosensory_response_tMap_' hemi '.smp']; % somatosensory_response_tMap name RH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% load map
cd(mapdir)
smp=xff(mapName);
smp_data=smp.Map.SMPData;
new_smp=smp;
new_smp.Map.SMPData=zeros(size(smp_data));
% assign values within the map by POIs
for p=1:length(POIs_data)
    new_smp.Map.SMPData(POIs_data(p).Vertices(MaskedV(p).masked))=dataval(p).LateralityAsubject;
end
