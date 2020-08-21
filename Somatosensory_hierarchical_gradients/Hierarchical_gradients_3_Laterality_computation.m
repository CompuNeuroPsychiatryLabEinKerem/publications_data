% Noam Saadon-Grosman
% August 2020
% This code computes laterality in each vertex of each subject in each
% stimulus direction that passed the threshold for significant
% somatosensory response. Laterality is defined as the normalized difference
%between contralateral and ipsilateral response
% The output stracture is arranged by parcellation areas.
%%% input: Files are in BrainVoyager format, for matrices extracted from
%%% other formats exclude the xff command. Fill in names and directories
%MTCs- (mesh time course) a nxm matrix witn n-number of volumes (TR), m-mesh
%vertices. Both ipsi and contralateral
%RTC- response prediction (boxar function convolved with the HRF)
% Somatosensory_response_tMap
% file with POIs (patch of intrest)
%%% output: A matlab .mat file 
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed parameters %%%%%%%%%%%%%%%%%%%%
lagmax=8; % number of TRs within stimulus block 
threshold=0.01/16; % significance threshold
hemi='LH'; % hemisphere
protocol={'lip','toe'};
%%%%%%%%%%%%%%%%%% names and directories %%%%%%%%%%%%%%%%%%%%%%%%%
mtcdirMain='C:\Users\OWNER';% MTCs directory
poidir='C:\Users\OWNER';% POIs directory
rtcdir='C:\Users\OWNER';% RTC directory
MAINoutput='C:\Users\OWNER';% output directory
if isempty(strfind(hemi,'LH'))
    poiName='POIs_RH.poi'; 
    MaskDir='C:\Users\OWNER'; % somatosensory_response_tMap directory RH;
    output=[MAINoutput '\RH'];
    contra_BS='Left_BS';
    ipsi_BS='Right_BS';
else
    poiName='POIs_LH.poi'; 
    MaskDir='C:\Users\OWNER'; % somatosensory_response_tMap directory LH;
    output=[MAINoutput '\LH'];
    contra_BS='Right_BS';
    ipsi_BS='Left_BS';
end
Group_Mask=['Somatosensory_response_tMap_' hemi '.smp']; % somatosensory_response_tMap name RH;
rtcName ='First_predictor_one_BS.rtc';
%%%%%%%%%%%%%%%%%%%%%%%% POI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(poidir)
POIs=xff(poiName);
POIs_data=POIs.POI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load template map
cd(mapdir)
smp=xff(mapName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RTC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(rtcdir);
rtc= xff(rtcName);
rtcData= rtc.RTCMatrix;
 %%%%%%%%%%%%%%%%%%%%%% MTC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MTCs_lip_contra= dir(fullfile([mtcdirMain 'Start_' protocol{1} '\' contra_BS],'*.mtc'));
    MTCs_toe_contra = dir(fullfile([mtcdirMain 'Start_' protocol{2} '\' contra_BS],'*.mtc'));
    MTCs_lip_ipsi= dir(fullfile([mtcdirMain 'Start_' protocol{1} '\' ipsi_BS '\Ipsilateral'],'*.mtc'));
    MTCs_toe_ipsi = dir(fullfile([mtcdirMain 'Start_' protocol{2} '\' ipsi_BS '\Ipsilateral'],'*.mtc'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open the group map to Mask with for significance
cd(MaskDir);
mask=xff(Group_Mask);
maskData=mask.Map.SMPData;
for p=1:length(POIs_data)
    to_keep=[];
    for i=1:length(POIs_data(p).Vertices)
        if maskData(POIs_data(p).Vertices(i))~=0
            to_keep=[to_keep,i];
            MaskedV(p).masked=to_keep;
            lkeep(p)=length(to_keep);
            loriginal(p)=length(POIs_data(p).Vertices);
        end
    end
end

% Loop on all MTCs in the folder
for m=1:(length(MTCs_lip_contra)+length(MTCs_toe_contra))
    if m<=length(MTCs_lip_contra)  % check if lip or toe folder
        f=1;
        MTCs_contra=MTCs_lip_contra;
        MTCs_ipsi=MTCs_lip_ipsi;
        n=m;
    else
        f=2;
        MTCs_contra=MTCs_toe_contra;
        MTCs_ipsi=MTCs_toe_ipsi;
        n=m-length(MTCs_lip_contra);
    end
    % load data
    cd([mtcdirMain 'Start_' protocol{f} '\' contra_BS])
    mtcNameC=MTCs_contra(n).name;
    mtcC=xff(mtcNameC);
    mtcdataC=mtcC.MTCData;
    cd([mtcdirMain 'Start_' protocol{f} '\' ipsi_BS '\Ipsilateral'])
    mtcNameI=MTCs_ipsi(n).name;
    mtcI=xff(mtcNameI);
    mtcdataI=mtcI.MTCData;
    % Loop on all POIs
    for p=1:length(POIs_data)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% MTC POI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        poiV=POIs_data(p).Vertices;
        %%%%%%%%%%%%%%%%%%%%%%%%%%% calculate selectivity %%%%%%%%%%%%%%%%
        POIdataC=mtcdataC(:,poiV);
        POIdataI=mtcdataI(:,poiV);
        rValuesContra=[]; rValuesIpsi=[]; pValuesContra=[]; pValuesIpsi=[];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correlation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for lag=1:lagmax 
            shiftedRtc= circshift(rtcData,lag);
            [rContra pContra]= corr(POIdataC, shiftedRtc);
            [rIpsi pIpsi]= corr(POIdataI, shiftedRtc);
            rValuesContra = [rValuesContra, rContra];
            rValuesIpsi = [rValuesIpsi,rIpsi];
            pValuesContra = [pValuesContra, pContra];
            pValuesIpsi = [pValuesIpsi,pIpsi];
        end
        rValuesContra(rValuesContra<0)=0;
        rValuesIpsi(rValuesIpsi<0)=0;
        % find max values and p min
        minpContra=min(pValuesContra');
        minpIpsi=min(pValuesIpsi');
        [maxcorContra lagContra]=max(rValuesContra');
        [maxcorIpsi lagIpsi]=max(rValuesIpsi');
        POI_Laterality=(maxcorContra-maxcorIpsi)./(maxcorContra+maxcorIpsi);
        toAverage=[];
        for i=1:size(POIdataC,2)
            if min(minpContra(i),minpIpsi(i))<threshold&&maskData(POIs_data(p).Vertices(i))~=0
                toAverage=[toAverage,POI_Laterality(i)];
            end
        end
    
        if isempty(toAverage)
            toAverage=0;
            sub_weight_lat(p,m)=0;
        else
          toAverage(isnan(toAverage))=0;
          sub_weight_lat(p,m)=length(toAverage(toAverage~=0))./(length(MaskedV(p).masked));
          MEAN_Laterality=mean(toAverage(toAverage~=0)); % for each subject averaged leterality within POI across significant vertices
        end
        
        % save all POI data with Group significance filter
        datavalMs(p).name=POIs_data(p).Name;
        datavalMs(p).rValuesContra(:,:,m)=rValuesContra(MaskedV(p).masked,:);
        datavalMs(p).rValuesIpsi(:,:,m)=rValuesIpsi(MaskedV(p).masked,:);
        datavalMs(p).pValuesContra(:,:,m)=pValuesContra(MaskedV(p).masked,:);
        datavalMs(p).pValuesIpsi(:,:,m)=pValuesIpsi(MaskedV(p).masked,:);
        datavalMs(p).Laterality(:,m)=POI_Laterality(MaskedV(p).masked);
        datavalMs(p).MEAN_Laterality(m)=MEAN_Laterality; 
    end
end
cd(MAINoutput);
save(['Laterality_' protocol{1} '_and_' protocol{2} '_' hemi],'datavalMs');

