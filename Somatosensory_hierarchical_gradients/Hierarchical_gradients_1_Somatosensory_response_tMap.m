% Noam Saadon-Grosman
% August 2020
clear all;
close all;
% This code takes mesh (surface) time courses of all participants and the two stimulus
% directions (star_lip and start_toe) and computes a random effect group t-values map
% corrected for multiple comparisons and masked for significance.
%%% input: Files are in BrainVoyager format, for matrices extracted from
%%% other formats exclude the xff command. Fill in the right names and directories
%MTCs- (mesh time course) a nxm matrix witn n-number of volumes (TR), m-mesh
%vertices
%RTC- response prediction (boxar function convolved with the HRF)
% Template BrainVoyager map
%%% output: Somatosensory_response_tMap
%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed parameters %%%%%%%%%%%%%%%%%%%%
lagmax=8; % number of TRs within stimulus block 
body_side='RIGHT_BS';
hemi='LH';
DF=137-2; % degrees of freedom
Sig=0.05/16; % significance threshold
%%%%%%%%%%%%%%%%%% names and directories %%%%%%%%%%%%%%%%%%%%%%%%%
rtcdir='C:\Users\OWNER';% RTC directory
rtcName='First_predictor_one_BS.rtc';
mtcdir_lip='C:\Users\OWNER';% MTCs start_lip directory
mtcdir_toe='C:\Users\OWNER';% MTCs start_toe directory
smpdir='C:\Users\OWNER';% template map directory
smpPname='Template_GLM_LH_fsaverage.smp';
output='C:\Users\OWNER';% output directory
saveName=['Somatosensory_response_tMap_' hemi];% name to save output map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RTC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(rtcdir);
rtc = xff(rtcName);
rtcData= rtc.RTCMatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start lip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(mtcdir_lip)
MTCs_lip = dir(fullfile(mtcdir_lip ,'*.mtc'));
    tempMTC=xff(MTCs_lip(1).name);
    tempdata=tempMTC.MTCData;
    nVoxels=size(tempdata,2);
% define storage matrices
Rvsub_lip=zeros(nVoxels,lagmax,length(MTCs_lip)); % 3D matrix with: vertices x all correlation values x all subjects
Pvsub_lip=zeros(nVoxels,lagmax,length(MTCs_lip));
% loop on all vtcs\subjects
for v=1:length(MTCs_lip)
    tempMTC=xff(MTCs_lip(v).name);
    tempdata=tempMTC.MTCData;
     rValues = [];
     pValues = [];
        for lag=1:lagmax
            shiftedRtc = circshift(rtcData,lag);
            [r p] = corr(tempdata, shiftedRtc);
            rValues = [rValues, r];
            pValues= [pValues, p];
        end
        Rvsub_lip(:,:,v)=rValues;
        Pvsub_lip(:,:,v)=pValues;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start toe %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the exact procedure as in start lip, needs to be done separately due to
% the opposite direction
cd(mtcdir_toe)
MTCs_toe = dir(fullfile(mtcdir_toe ,'*.mtc'));
% define storage matrices
Rvsub_toe=zeros(nVoxels,lagmax,length(MTCs_toe));
Pvsub_toe=zeros(nVoxels,lagmax,length(MTCs_lip));
% loop on all mtcs\subjects
for v=1:length(MTCs_toe)
    tempMTC=xff(MTCs_toe(v).name);
    tempdata=tempMTC.MTCData;
     rValuest = [];
     pValuest = [];
        for lag=1:lagmax
            shiftedRtc = circshift(rtcData,lag);
            [rt pt] = corr(tempdata, shiftedRtc);
            rValuest = [rValuest, rt];
            pValuest = [pValuest, pt];
        end
        Rvsub_toe(:,:,v)=rValuest;
        Pvsub_toe(:,:,v)=pValuest;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save pvalues %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a matrix of all p values (the minimum out of the 8) for each TC-
% first all subjects in start_lip and then all subjects in start_toe
PPvsub_lip=permute(Pvsub_lip,[2,1,3]);
AllpValue_lip=squeeze(min(PPvsub_lip));
PPvsub_toe=permute(Pvsub_toe,[2,1,3]);
AllpValue_toe=squeeze(min(PPvsub_toe));
AllpValue=[AllpValue_lip,AllpValue_toe];
cd(output)
save(['AllpValue_' hemi],'AllpValue');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% average lip and toe %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in order to generate an averaged map we need to flip start toe
% correlation values
Rvsub_toe_flip=Rvsub_toe(:,(lagmax:-1:1),:);
% average lip and toe correlation values
Rvsub_toeP=permute(Rvsub_toe_flip,[2,1,3]);
Rvsub_lipP=permute(Rvsub_lip,[2,1,3]);
R=[];
R(1,:,:,:)=Rvsub_toeP;
R(2,:,:,:)=Rvsub_lipP;
R_liptoe=mean(R);
R_liptoe=squeeze(R_liptoe(1,:,:,:));
%%%%%%%%%%%%%%%%%%%%%%% statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transform r values to t values
t_liptoe=(R_liptoe).*sqrt(DF./(1-(R_liptoe.^2))); 
% for each vertex of each lag apply t-test on all subject's t-values 
pgroup=zeros(size(t_liptoe,1),size(t_liptoe,2));
tgroup=zeros(size(t_liptoe,1),size(t_liptoe,2));
for i=1:size(t_liptoe,1)
    for j=1:size(t_liptoe,2)
        temp=t_liptoe(i,j,:);
        temp=squeeze(temp(1,:,:));
        [h,p,ci,stats]=ttest(temp,0,'Tail','right');
        pgroup(i,j)=p;
        tgroup(i,j)=stats.tstat;
    end
end
%%%%%%%%%%%%%%%%%%%%%%% correct for multiple comparisons %%%%%%%%%%%%%%%%%
pgroupFDR=zeros(size(t_liptoe,1),size(t_liptoe,2));
% apply FDR on each lag seperatly
for l=1:size(pgroup,1)
    templag=pgroup(l,:);
    pgroupFDR(l,:)=mafdr(templag,'BHFDR',true); 
end
% choose the lag with the smallest probability
ChosenTvalue=zeros(1,size(pgroup,2));
for i=1:size(pgroup,2)
    tempL=pgroupFDR(:,i);
    tempt=tgroup(:,i);
    [minp indp]=min(tempL);
    ChosenTvalue(i)=tempt(indp);
end
% mask T values with p significance
ChosenTvalue(Chosen>Sig)=0;
% create map
 toMap=ChosenTvalue; % set back map dimentions
 cd(smpdir);
        smp=xff(smpPname);
        new_smp=smp;
        new_smp.Map.SMPData=toMap;
        new_smp.Map.Name=saveName;
        snameO=[saveName '.smp'];
        cd(output);
        new_smp.SaveAs(snameO); 
               
