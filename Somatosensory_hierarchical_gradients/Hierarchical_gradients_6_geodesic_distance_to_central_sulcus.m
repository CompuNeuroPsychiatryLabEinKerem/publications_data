% Noam Saadon-Grosman
% August 2020
% This code computes the geodesic distance between each vertex to the closest vertex in area
% 3a. We use an online code (https://code.google.com/archive/p/geodesic/) to compute the distance, separately for each
% gross anatomical region. The computation is for the right and then for
% the left hemisphere
% It also computes selectivity and laterality values within distance
% quantiles
%%% input: some files are in BrainVoyager format, for matrices extracted from
%%% other formats exclude the xff command. Fill in names and directories
% Selectivity ans Laterality maps
% surface mesh
% file with POIs (patch of intrest)- area 3a and gross anatomical regions
clear all;
close all;
% Define directories and upload data
output='C:\Users\OWNER'; % output directory
% Load POIs of all 4 gross anatomical regions and area 3a
cd('C:\Users\OWNER\') % POIs directory
POI3a_RH=xff('3a_boundaries_of_CS_RH.poi'); % area 3a
POI3a_LH=xff('3a_boundaries_of_CS_LH.poi'); 
POIs_RH=xff('HCP_RH_somatosensory_All_4_gross.poi'); 
POIs_LH=xff('HCP_LH_somatosensory_All_4_gross.poi'); 
% area 3a
A3aPOI_RH=POI3a_RH.POI(1);
A3avertices_RH=A3aPOI_RH.Vertices;
A3aPOI_LH=POI3a_LH.POI(1);
A3avertices_LH=A3aPOI_LH.Vertices;
% load surface information- fsaverage half distance between pial and white
% matter surfaces
cd('C:\Users\OWNER\') % surface directory
srf_RH=load('fsaverage_mid_rh_bv.mat');
srf_LH=load('fsaverage_mid_lh_bv.mat');
% for geodesic algorithm:
% RH
surface_RH.X=srf_RH.coords(:,1);
surface_RH.Y=srf_RH.coords(:,2);
surface_RH.Z=srf_RH.coords(:,3);
surface_RH.TRIV=srf_RH.triangles+1;
shapeRH.surface=surface_RH;
% LH
surface_LH.X=srf_LH.coords(:,1);
surface_LH.Y=srf_LH.coords(:,2);
surface_LH.Z=srf_LH.coords(:,3);
surface_LH.TRIV=srf_LH.triangles+1;
shapeLH.surface=surface_LH;

% load selectivity and laterality maps 
cd('C:\Users\OWNER\') % selectivity maps directory
SelMap_RH=xff('RH_selectivity_map.smp');
SelMap_LH=xff('LH_selectivity_map.smp');
cd('C:\Users\OWNER\') % laterality maps directory
LatMap_RH=xff('RH_laterality_map.smp');
LatMap_LH=xff('LH_laterality_map.smp');
SelData_RH=SelMap_RH.Map.SMPData;
LatData_RH=LatMap_RH.Map.SMPData;
SelData_LH=SelMap_LH.Map.SMPData;
LatData_LH=LatMap_LH.Map.SMPData;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% group mask %%%%%%%%%%%%%%%%%%%%%%%%%

for p=1:length(POIs_LH.POI)
    to_keepRH=[];
    to_keepLH=[];
    for i=1:length(POIs_RH.POI(p).Vertices)
        if SelData_RH(POIs_RH.POI(p).Vertices(i))~=0
            to_keepRH=[to_keepRH,i];
            MaskedV_RH(p).masked=to_keepRH;
        end
    end
    for i=1:length(POIs_LH.POI(p).Vertices)
        if SelData_LH(POIs_LH.POI(p).Vertices(i))~=0
            to_keepLH=[to_keepLH,i];
            MaskedV_LH(p).masked=to_keepLH;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop on all 4 gross anatomical regions, in each one find the geodesic distance
% of each vertex to all vertices in area 3a and choose the closest path,
% save all in a stracture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geodesic code directory
cd('C:\Users\OWNER\') % code files directory
% compute distance to all vertices in the mesh to area 3a
distancesAllmesh_RH= geodesics_to_all(shapeRH, A3avertices_RH');
figure(1)
for p=1:length(POIs_RH.POI)
    tempPOI=POIs_RH.POI(p);
    POIvertices=tempPOI.Vertices;
    distances=distancesAllmesh_RH(POIvertices(MaskedV_RH(p).masked));
    POIGeoDist_RH(p).distance=distances;
    POIGeoDist_RH(p).name=tempPOI.Name;
    % selectivity and laterality
    Seltemp=SelData_RH(POIvertices(MaskedV_RH(p).masked));
    Lattemp=LatData_RH(POIvertices(MaskedV_RH(p).masked));
    POIGeoDist_RH(p).Sel=Seltemp;
    POIGeoDist_RH(p).Lat=Lattemp;
    % visualization
    SelnonZ=find(Seltemp~=0);
    subplot(2,4,p)
    hold on
    plot(distances(SelnonZ),Seltemp(SelnonZ),'.','Color',[0 0 0])
    POIGeoDist_RH(p).CorSel=corr(distances(SelnonZ),Seltemp(SelnonZ));
    subplot(2,4,p+4)
    hold on
    LatnonZ=find(Lattemp~=0);
    plot(distances(LatnonZ),Lattemp(LatnonZ),'.','Color',[0 0 0])
    POIGeoDist_RH(p).CorLat=corr(distances(LatnonZ),Lattemp(LatnonZ));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distancesAllmesh_LH= geodesics_to_all(shapeLH, A3avertices_LH');
figure(2)
for p=1:length(POIs_LH.POI)
    tempPOI=POIs_LH.POI(p);
    POIvertices=tempPOI.Vertices;
    distances=distancesAllmesh_LH(POIvertices(MaskedV_LH(p).masked));
    POIGeoDist_LH(p).distance=distances;
    POIGeoDist_LH(p).name=tempPOI.Name;
    Seltemp=SelData_LH(POIvertices(MaskedV_LH(p).masked));
    Lattemp=LatData_LH(POIvertices(MaskedV_LH(p).masked));
    POIGeoDist_LH(p).Sel=Seltemp;
    POIGeoDist_LH(p).Lat=Lattemp;
    SelnonZ=find(Seltemp~=0);
    subplot(2,4,p)
    hold on
    plot(distances(SelnonZ),Seltemp(SelnonZ),'.','Color',[0 0 1])
    POIGeoDist_LH(p).CorSel=corr(distances(SelnonZ),Seltemp(SelnonZ));
    subplot(2,4,p+4)
    hold on
    LatnonZ=find(Lattemp~=0);
    plot(distances(LatnonZ),Lattemp(LatnonZ),'.','Color',[0 0 1])
    POIGeoDist_LH(p).CorLat=corr(distances(LatnonZ),Lattemp(LatnonZ));
end
close all;
% quantiles
for p=1:length(POIGeoDist_RH)
    test_data_RH=POIGeoDist_RH(p).distance;
    test_data_LH=POIGeoDist_LH(p).distance;
    YQ_RH = [quantile(test_data_RH,10) max(test_data_RH)];
    YQ_LH = [quantile(test_data_LH,10) max(test_data_LH)];
    for i=1:length(YQ_RH)-1
        disind_RH=find(YQ_RH(i)<=test_data_RH&test_data_RH<YQ_RH(i+1));
        disind_LH=find(YQ_LH(i)<=test_data_LH&test_data_LH<YQ_LH(i+1));
        allsel_RH=double(POIGeoDist_RH(p).Sel(disind_RH));
        alllat_RH=double(POIGeoDist_RH(p).Lat(disind_RH));
        allsel_LH=double(POIGeoDist_LH(p).Sel(disind_LH));
        alllat_LH=double(POIGeoDist_LH(p).Lat(disind_LH));
        seltemp_RH=allsel_RH(allsel_RH~=0);
        lattemp_RH=alllat_RH(alllat_RH~=0);
        seltemp_LH=allsel_LH(allsel_LH~=0);
        lattemp_LH=alllat_LH(alllat_LH~=0);
        sel_RH(i)=mean(seltemp_RH);
        lat_RH(i)=mean(lattemp_RH);
        sel_LH(i)=mean(seltemp_LH);
        lat_LH(i)=mean(lattemp_LH);
    end
    %%%%%%%%%%%%%%%%%%% linear regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mSel_RH=((mean(YQ_RH(1:end-1))*mean(sel_RH))-mean(YQ_RH(1:end-1).*sel_RH))/((mean(YQ_RH(1:end-1)))^2-mean(YQ_RH(1:end-1).^2));
    bSel_RH=mean(sel_RH)-mSel_RH*mean(YQ_RH(1:end-1));
    mLat_RH=((mean(YQ_RH(1:end-1))*mean(lat_RH))-mean(YQ_RH(1:end-1).*lat_RH))/((mean(YQ_RH(1:end-1)))^2-mean(YQ_RH(1:end-1).^2));
    bLat_RH=mean(lat_RH)-mLat_RH*mean(YQ_RH(1:end-1));
    mSel_LH=((mean(YQ_LH(1:end-1))*mean(sel_LH))-mean(YQ_LH(1:end-1).*sel_LH))/((mean(YQ_LH(1:end-1)))^2-mean(YQ_LH(1:end-1).^2));
    bSel_LH=mean(sel_LH)-mSel_LH*mean(YQ_LH(1:end-1));
    mLat_LH=((mean(YQ_LH(1:end-1))*mean(lat_LH))-mean(YQ_LH(1:end-1).*lat_LH))/((mean(YQ_LH(1:end-1)))^2-mean(YQ_LH(1:end-1).^2));
    bLat_LH=mean(lat_LH)-mLat_LH*mean(YQ_LH(1:end-1));
    POIGeoDist_LH(p).mSel=mSel_LH;
    POIGeoDist_LH(p).mLat=mLat_LH;
    POIGeoDist_RH(p).mSel=mSel_RH;
    POIGeoDist_RH(p).mLat=mLat_RH;

    figure(3)
    subplot(2,4,p)
    plot(YQ_RH(1:end-1),sel_RH,'.','MarkerSize',30)
    hold on
    plot((0:YQ_RH(end)+5),mSel_RH.*(0:YQ_RH(end)+5)+bSel_RH,'LineWidth',1.5,'Color',[0 0 0]);
    xlabel('Distance [mm]');
    subplot(2,4,p+4)
    plot(YQ_RH(1:end-1),lat_RH,'.','MarkerSize',30)
    hold on
    plot((0:YQ_RH(end)+5),mLat_RH.*(0:YQ_RH(end)+5)+bLat_RH,'LineWidth',1.5,'Color',[0 0 0]);
    figure(4)
    subplot(2,4,p)
    plot(YQ_LH(1:end-1),sel_LH,'.','MarkerSize',30)
    hold on
    plot((0:YQ_LH(end)+5),mSel_LH.*(0:YQ_LH(end)+5)+bSel_LH,'LineWidth',1.5,'Color',[0 0 0]);
    subplot(2,4,p+4)
    plot(YQ_LH(1:end-1),lat_LH,'.','MarkerSize',30)
    hold on
    plot((0:YQ_LH(end)+5),mLat_LH.*(0:YQ_LH(end)+5)+bLat_LH,'LineWidth',1.5,'Color',[0 0 0]);

end

cd(output)
save('Vertices_Sel_Lat_to_Geodesic' ,'distancesAllmesh_RH','distancesAllmesh_LH' ,'POIGeoDist_RH','POIGeoDist_LH');

