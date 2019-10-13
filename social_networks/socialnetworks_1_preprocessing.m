% Pre-processing of subjects' data in Brainvoyager


%% Starting brainvoyager from Matlab
bvqx = actxserver('BrainVoyager.BrainVoyagerScriptAccess.1');


%% Defining working directories
subjects_dicom_dir = 'F:\Social_networks_analysis\DICOMS\';     % A directory containing subdirectories for each subject, each containing one subject's dicom files - should be organized in advance
subjects_analysis_dir = 'F:\Social_networks_analysis\Analysis_2mm_smoothing\';      % A directory where the preprocessed data will be put
subject_names = dir(subjects_dicom_dir); subject_names = subject_names(3:end);      % Getting the names of the subjects directories
% Creating the analysis directory if it doesn't exist
if ~exist(subjects_analysis_dir, 'dir'), mkdir(subjects_analysis_dir); end


%% Defining MRI parameters
% Anatomical run
xres_T1=256; yres_T1=256;
nrslices_T1=160;
swap_T1=0;

% Functional runs
nrVolsInImg_func = 1;
byteswap_func = false;
sizeX_func = 64; sizeY_func = 64;
nrSlices_func = 37;
mosaicSizeX_func = 448; mosaicSizeY_func = 448;

% Smoothing level in mm
FWHM_smoothing = 2;


%% Renaming all DICOM files to BrainVoyager format
disp('Renaming DICOM files')
for s = 1:length(subject_names)     % Looping over all subjects
    subj = subject_names(s).name;
    curr_folder = fullfile(subjects_dicom_dir, subj);
    subfolders = dir(curr_folder); subfolders = subfolders(3:end);
    for d = 1:length(subfolders)
        dicom_dir = fullfile(curr_folder, subfolders(d).name);
        bvqx.RenameDicomFilesInDirectory(dicom_dir);
    end
end



%% creating BrainVoyager FMR/VMR files from functional and anatomical runs dicoms
disp('Creating FMR and VMR files')
for s = 1:length(subject_names)     % Looping over all subjects
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);
    % Creating the subject's directory if it doesn't exist yet
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    
    curr_folder = fullfile(subjects_dicom_dir, subj);
    subfolders = dir(curr_folder); subfolders = subfolders(3:end);
    counter_runs = 1; 
    for d = 1:length(subfolders)    % Looping over runs
        dicom_dir = fullfile(curr_folder, subfolders(d).name);
        current_files=dir(fullfile(dicom_dir,'*.dcm'));
        
        fileType = 'DICOM';
        firstFile = fullfile(dicom_dir, current_files(1).name);
        bytesperpixel = 2;
        
        if ~isempty(strfind(subfolders(d).name,'T1'))   % if it is an anatomical T1 file
            % creating VMR
            vmrproject = bvqx.CreateProjectVMR(fileType, firstFile, nrslices_T1, swap_T1, xres_T1, yres_T1, bytesperpixel);
            vmrproject.SaveAs(fullfile(output_dir,[subj '.vmr']));
            
        elseif ~isempty(strfind(subfolders(d).name,'PEOPLE'))    % if it is a functional file
            nrOfVols = length(current_files);
            skipVols = 0;
            createAMR = true;
            
            stcprefix = [subj '_' num2str(counter_runs)];
            fmr_output_filename = fullfile(output_dir,[subj '_' num2str(counter_runs) '.fmr']);
            counter_runs = counter_runs + 1;
            
            % create FMR
            fmr = bvqx.CreateProjectMosaicFMR(fileType, firstFile,...
                nrOfVols, skipVols, createAMR, nrSlices_func,...
                stcprefix, byteswap_func, mosaicSizeX_func, mosaicSizeY_func, bytesperpixel,...
                output_dir, nrVolsInImg_func, sizeX_func, sizeY_func);
            fmr.SaveAs(fmr_output_filename);
            
        end
    end
end



%% Functional data (FMR) preprocessing

% Slice timing correction
disp('Slice timing correction')
for s = 1:length(subject_names)     % Looping over all subjects
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);
    
    % finding relevant fmr files
    fmr_files = dir(fullfile(output_dir,'*.fmr')); 
    fmr_files = {fmr_files.name}; for i=1:length(fmr_files), fmr_files{i} = fullfile(output_dir, fmr_files{i}); end
    ix = 1:length(fmr_files);
    for i=1:length(ix)
        if strcmp(fmr_files{i}(end-11:end-4),'firstvol')
            ix(ix==i) = [];
        end
    end
    fmr_files = fmr_files(ix);    
        
    for i = 1:length(fmr_files)
        fmr = bvqx.OpenDocument(fmr_files{i});
        interpolation_type = 1;  % 1 - cubic spline interpolation
        % Performing the slice timing correction
        fmr.CorrectSliceTimingUsingTimeTable(interpolation_type);    
         % delete original files without deleting firstvol_as_anat
         delete(fmr_files{i});
        delete([fmr_files{i}(1:end-3) 'stc']);
    end
end


% Motion correction
disp('FMR motion correction')
for s = 1:length(subject_names)     % Looping over all subjects
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);
    
    % finding relevant fmr files
    fmr_files = dir(fullfile(output_dir,'*_SCCTBL.fmr'));
    fmr_files = {fmr_files.name}; for i=1:length(fmr_files), fmr_files{i} = fullfile(output_dir, fmr_files{i}); end

    for i=1:length(fmr_files)
        fmr = bvqx.OpenDocument(fmr_files{i});
        TargetVolume = 1;           % the number of volume in the series to align to
        Interpolation_type = 2;     % 0 and 1 - trilinear detection and trilinear interpolation, 2: trilinear detection and sinc interpolation or 3: sinc detection of motion and sinc interpolation
        UseFullDataSet = false;     % false is the default in the GUI
        MaxNumIterations = 100;     % 100 is the default in the GUI
        GenerateMovie = true;
        GenerateLogFile = true;     % creates a log file with the movement parameters
        % Performing the motion correction
        fmr.CorrectMotionEx(TargetVolume, Interpolation_type, UseFullDataSet,...
            MaxNumIterations, GenerateMovie, GenerateLogFile);
        fmr.Remove; % close or remove input FMR
    end

end



% High-pass filtering
disp('High pass filtering')
for s = 1:length(subject_names)     % Looping over all subjects
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);
    
    % finding relevant fmr files
    fmr_files = dir(fullfile(output_dir,'*_3DMCTS.fmr'));
    fmr_files = {fmr_files.name}; for i=1:length(fmr_files), fmr_files{i} = fullfile(output_dir, fmr_files{i}); end

    for i=1:length(fmr_files)
        fmr = bvqx.OpenDocument(fmr_files{i});
        NumCycles = 2;
        % Performing the high pass filtering
        fmr.TemporalHighPassFilterGLMFourier(NumCycles);
        fmr.Remove; % close or remove input FMR
    end
end


% Smoothing
disp('Smoothing')
for s = 1:length(subject_names)     % Looping over all subjects
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);
    
    % finding relevant fmr files
    fmr_files = dir(fullfile(output_dir,'*_THPGLMF2c.fmr'));
    fmr_files = {fmr_files.name}; for i=1:length(fmr_files), fmr_files{i} = fullfile(output_dir, fmr_files{i}); end

    for i = 1:length(fmr_files)
        fmr = bvqx.OpenDocument(fmr_files{i});
        % Performing the smoothing
        fmr.SpatialGaussianSmoothing(FWHM_smoothing, 'mm');
    end
end



%% Anatomical data (VMR) preprocessing

% Inhomogeneity correction and skull stripping
disp('Inhomogeneity correction and skull stripping')
for s = 1:length(subject_names)     % Looping over all subjects
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);

    vmr_file = dir(fullfile(output_dir,'*.vmr'));
    vmr_file = fullfile(output_dir, vmr_file(1).name); 
    vmr = bvqx.OpenDocument(vmr_file);
    % Performing the inhomogeneity correction
    vmr.CorrectIntensityInhomogeneities();
end


% Normalization to MNI space
disp('Normalizing')
for s = 1:length(subject_names)     % Looping over all subjects
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);

    vmr_file = dir(fullfile(output_dir,'*IIHC.vmr'));
    vmr_file = fullfile(output_dir, vmr_file(1).name);
    vmr = bvqx.OpenDocument(vmr_file);
    % Performing the normalization
    vmr.NormalizeToMNISpace();
end



%% MANUALLY CHECK NORMALIZATION RESULTS BEFORE CLOSING THE WINDOWS!!!
% (USE F8 TO SEE OVERLAY)


%% Co-registeration of functional and anatomical data

disp('Co-registering...')
for s = 1:length(subject_names)     % Looping over all subjects
    disp(s);
    subj = subject_names(s).name;
    output_dir = fullfile(subjects_analysis_dir, subj);
    
    % Finding relevant fmr files
    fmr_files = dir(fullfile(output_dir, '*SD3DSS*mm.fmr'));    % finding smoothed files
    vmr_file = dir(fullfile(output_dir,'*IIHC.vmr'));
    vmr_file = fullfile(output_dir, vmr_file(1).name);

    % Coregistration using intensity gradient-based matching
    useAttachedAMR = 1;
    for i=1:length(fmr_files)
        vmr = bvqx.OpenDocument(vmr_file);
        fmr_filename = fullfile(output_dir, fmr_files(i).name);
        % Performing the co-registration
        vmr.CoregisterFMRToVMR(fmr_filename, useAttachedAMR);
    end
end


%% MANUALLY CHECK COREGISTRATION RESULTS BEFORE CLOSING THE WINDOWS!!!
% IF PROBLEMATIC, REDO MANUALLY...


%% Creating BrainVoyager VTC files from each functional run (functional data co-registered and normalized to MNI space)

disp('Creating VTCs...')
for s = 1:length(subject_names)     % Looping over all subjects
    subj = subject_names(s).name;
    disp(subj);
    output_dir = fullfile(subjects_analysis_dir, subj);
        
    vmr_file = dir(fullfile(output_dir, '*MNI.vmr'));       % MNI-transformed anatomical file
    vmr_file = fullfile(output_dir, vmr_file(1).name);
    vmr = bvqx.OpenDocument(vmr_file);
    
    MNI_transform_file = dir(fullfile(output_dir,'*MNI_a12.trf')); MNI_transform_file = fullfile(output_dir, MNI_transform_file(1).name);
    
    fmr_files = dir(fullfile(output_dir,'*SD3DSS*mm.fmr'));     % Smoothed FMR files
    IA_files = dir(fullfile(output_dir,'*_IIHC_IA.trf'));       % Co-registration transformation parameters files
    FA_files = dir(fullfile(output_dir,'*_IIHC*_FA.trf'));      % Co-registration transformation parameters files
    
    for i=1:length(fmr_files)
        fmr_filename = fullfile(output_dir, fmr_files(i).name); 
        IA_filename = fullfile(output_dir, IA_files(i).name); 
        FA_filename = fullfile(output_dir, FA_files(i).name);        
        Datatype = 2;       % 1 - integer, 2 - float (GUI default)
        Resolution = 3;     % resolution relative to VMR - 3x3x3
        Interpolation = 1;    % 0 for nearest neighbor , 1 for trilinear (GUI default) , 2 for sinc interpolation
        Intensity_threshold = 100;    % intensity threshold for voxels
        VTC_name = [fmr_filename(1:end-4) '_MNI.vtc'];
        
        % Creating the VTC file
        vmr.CreateVTCInMNISpace(fmr_filename, IA_filename, FA_filename, MNI_transform_file, ...
            VTC_name, Datatype, Resolution, Interpolation, Intensity_threshold);
    end
end


