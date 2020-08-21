% Noam Saadon-Grosman
% August 2020
% This code computes selectivity in each vertex of each subject in each
% stimulus direction that passed the threshold for significant
% somatosensory response. Selectivity is defined as 1/standard deviation of
% a Gaussian fitted to the response Event Related Averaging. The output
% stracture is arranged by parcellation areas.
%%% input: Files are in BrainVoyager format, for matrices extracted from
%%% other formats exclude the xff command. Fill in names and directories
%MTCs- (mesh time course) a nxm matrix witn n-number of volumes (TR), m-mesh
%vertices
%RTC- response prediction (boxar function convolved with the HRF)
% Somatosensory_response_tMap
% file with POIs (patch of intrest)
%%% output: A matlab .mat file 
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed parameters %%%%%%%%%%%%%%%%%%%%
rep=7; % number of stimulus-rest block repetitions in a run
intervals=9:18:137; % windows of stimulus-rest block
time_window=18; % number of TRs in each stimulus-rest block repetition
lagmax=8; % number of TRs within stimulus block 
hemi='RH'; 
protocol={'lip','toe'};
body_side='RBS';
dirBS='Right';
%%%%%%%%%%%%%%%%%% names and directories %%%%%%%%%%%%%%%%%%%%%%%%%
MTCs_lipdir='C:\Users\OWNER';% MTCs start_lip directory
MTCs_toedir='C:\Users\OWNER';% MTCs start_toe directory
poidir='C:\Users\OWNER';% POIs directory
rtcdir='C:\Users\OWNER';% RTC directory
output='C:\Users\OWNER';% output directory
if isempty(strfind(dirBS,'Right')) % left body stimulation 
    MaskDir='C:\Users\OWNER'; % somatosensory_response_tMap directory RH;
    poiName='POIs_RH.poi'; 
else
    MaskDir='C:\Users\OWNER'; % somatosensory_response_tMap directory LH;
    poiName='POIs_LH.poi'; 
end
    Group_Mask=['Somatosensory_response_tMap_' hemi '.smp']; % somatosensory_response_tMap name RH;

rtcName='First_predictor_one_BS.rtc';
%%%%%%%%%%%%%%%%%%%%%% MTC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(MTCs_lipdir)
MTCs_lip= dir(fullfile(MTCs_lipdir,'*.mtc'));
cd(MTCs_toedir)
MTCs_toe = dir(fullfile(MTCs_toedir,'*.mtc'));
%%%%%%%%%%%%%%%%%%%%%%%% POI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(poidir)
POIs=xff(poiName);
POIs_data=POIs.POI; 

%%%%%%%%%%%%%%%%%% somatosensory response map %%%%%%%%%%%%%%%%%%%
% open somatosensory response map
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RTC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(rtcdir);
rtc = xff(rtcName);
rtcData= rtc.RTCMatrix;
% Loop on all MTCs in the folder
for m=1:(length(MTCs_lip)+length(MTCs_toe))
    tic
    if m<=length(MTCs_lip)  % check if lip or toe
        f=1;
        cd(MTCs_lipdir)
        MTCs=MTCs_lip;
        n=m;
    else
        f=2;
        cd(MTCs_toedir)
        MTCs=MTCs_toe;
        n=m-length(MTCs_lip);
    end
    cd(MTCs_lip)
    mtcName=MTCs(n).name;
    mtc=xff(mtcName);
    mtcdata=mtc.MTCData;
    
    % Loop on all POIs
    for p=1:length(POIs_data)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% MTC POI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        poiV=POIs_data(p).Vertices(MaskedV(p).masked);
        %%%%%%%%%%%%%%%%%%%%%%%%%%% calculate selectivity %%%%%%%%%%%%%%%%
        POIdata=mtcdata(:,poiV);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correlation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rValues = []; % cross correlation values
        for lag=1:lagmax 
            shiftedRtc = circshift(rtcData,lag);
            r = corr(POIdata, shiftedRtc);
            rValues = [rValues, r];
        end
        % storage correlation values
        % if the MTC is start toe, switch the direction of correlation
        % values that it would fit start lip
        if f==1
            POIrValues=rValues';
        else
            POIrValues=rValues(:,lagmax:(-1):1)';
        end
        [maxcor maxlag]=max(POIrValues);
                
        % Event Related Averaging
        repTC=zeros(rep,time_window,length(POIdata)); % cut repetitions
        nPOIdata=zscore(POIdata); % normalize all vertices
        repTC=[];
        if size(nPOIdata,2)==1
            for j=1:rep % loop on repetitions
                repTC(j,:)=nPOIdata(intervals(j):intervals(j)+(time_window-1))';
            end
            mrepTC=mean(repTC); % mean time windows across repetitions
            mrepTC=squeeze(mrepTC(1,:));
        else
            for j=1:rep % loop on repetitions
                repTC(j,:,:)=nPOIdata(intervals(j):intervals(j)+(time_window-1),:);
            end
            mrepTC=mean(repTC); % mean time windows across repetitions
            mrepTC=squeeze(mrepTC(1,:,:));
        end
        % selectivity
        FailFit=0;
        amplitude=zeros(1,size(mrepTC,2));
        selectivity=zeros(1,size(mrepTC,2));
        preference=zeros(1,size(mrepTC,2));
        r_square=zeros(1,size(mrepTC,2));
        parfor i=1:size(mrepTC,2)
            temp=mrepTC(:,i)+abs(min(mrepTC(:,i)));
            if ~isempty(temp)
                [xData, yData] = prepareCurveData((1:length(temp))',temp)
                try
                    ft = fittype( 'gauss1' );
                    opts = fitoptions( ft );
                    opts.Display = 'Off';
                    opts.Lower = [-Inf -Inf 0];
                    opts.StartPoint = [0.604662089339839 6 1.78290936063073];
                    opts.Upper = [Inf Inf Inf];
                    [fitresult, gof] = fit( xData, yData, ft, opts )
                    amplitude(i)=fitresult.a1;
                    selectivity(i)=fitresult.c1;
                    preference(i)=fitresult.b1;
                    r_square(i)=gof.rsquare;
                catch
                    FailFit=FailFit+1
                end
            end
        end
        
               
        % save all POI data
        dataval(p).name=POIs_data(p).Name;
        dataval(p).rValues(:,:,m)=POIrValues;
        dataval(p).lag(:,m)=maxlag;
        dataval(p).cor(:,m)=maxcor;
        dataval(p).width_ERA(:,m)=selectivity;
        dataval(p).Prefered_ERA(:,m)=preference;
        dataval(p).FailFit(:,m)=FailFit;
        dataval(p).amplitude(:,m)=amplitude;
        dataval(p).r_square(:,m)=r_square;
    end
    toc
end
    
cd(output);    
save(['selectivity_' protocol{1} '_and_' protocol{2} '_' hemi],'dataval');
