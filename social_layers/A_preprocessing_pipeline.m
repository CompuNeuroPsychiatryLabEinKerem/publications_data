% This pre-processing pipeline is performed through the BrainVoyager 20.6
% scripting platform, along with NeuroElf v1.1.


%% Enabling brainvoyager scripting from path
bvqx = actxserver('BrainVoyager.BrainVoyagerScriptAccess.1');


%% Defining working directories and getting subjects' directories names
% Preprocessing starts from DICOM files; each subject has a separate
% directory with their files inside the main DICOM directory. In each dicom
% directory there are subfolders for each experimental run and a subfolder
% for anatomical data.
subjects_dicom_dir = 'social_main\1_orig_dicom\';    % the directory with the DICOM files for each subject, arranged in separate subject directories
subjects_output_data_dir = 'social_main\2_BV_prepro\';     % the directory where subjects' preprocessed data will be saved

% getting subjects' directories names
subject_names = dir(subjects_dicom_dir); subject_names = subject_names(3:end);



%% renaming DICOM files to BrainVoyager format
disp('Renaming DICOM files')
for s=1:length(subject_names)
    current_subj = subject_names(s).name;
    current_subject_dir = fullfile(subjects_dicom_dir, current_subj);
    
%     going over subfolders (each subject dicom directory has separate
%     subfolders for each run, including anatomy)
    subfolders = dir(current_subject_dir); subfolders = subfolders(3:end);
    for d=1:length(subfolders)
        dicom_subfolder = fullfile(current_subject_dir, subfolders(d).name);
        bvqx.RenameDicomFilesInDirectory(dicom_subfolder);
    end
end



%% creating FMR and VMR files from each fMRI run
% (FMR is brainvoyager's format for functional data, VMR for anatomical data)

anatomical_run_dirname_part = 'T1';
functional_run_dirname_part = '_EGO';
control_run_dirname_part = '_LEX';


disp('Creating FMR and VMR files')
for s=1:length(subject_names)
    disp(s);
    current_subj = subject_names(s).name;
    
    % creating a directory for each subject in the preprocessed data directory
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    mkdir(output_dir);
    
    % going over runs and converting them to FMR or VMR according to run type (functional / anatomical)
    current_subject_dicom_dir = fullfile(subjects_dicom_dir, current_subj);
    subfolders = dir(current_subject_dicom_dir); subfolders = subfolders(3:end);
    run_number_counter = 1;
    
    for d=1:length(subfolders)
        % identifying all DICOM files
        current_dicom_subfolder = fullfile(current_subject_dicom_dir, subfolders(d).name);
        current_dicom_files = dir(fullfile(current_dicom_subfolder, '*.dcm'));
        firstFile = fullfile(current_dicom_subfolder, current_dicom_files(1).name);
        
        if ~isempty(strfind(subfolders(d).name, anatomical_run_dirname_part))   % This is a directory with anatomical (T1) data, by the name of the subdirectory
            % defining parameters for conversion
            fileType = 'DICOM';
            bytesperpixel = 2;
            swap=0;
            xres=256; yres=256;
            nrslices=160;
            
            % creating VMR file from dicoms
            vmrproject = bvqx.CreateProjectVMR(fileType, firstFile, nrslices, swap, xres, yres, bytesperpixel);
            vmrproject.SaveAs(fullfile(output_dir,[current_subj '.vmr']));
            
        elseif (~isempty(strfind(subfolders(d).name, functional_run_dirname_part)) || ~isempty(strfind(subfolders(d).name, control_run_dirname_part)))
            % This is a directory with functional (BOLD) data (including control) by the name of the subdirectory
            % defining parameters for conversion
            fileType = 'DICOM';
            bytesperpixel = 2;
            skipVols = 0;
            createAMR = true;
            nrVolsInImg = 1;
            byteswap = false;
            sizeX = 96; sizeY = 96;
            nrSlices = 60;
            mosaicSizeX = 768; mosaicSizeY = 768;
            
            % finding the number of volumes in this functional run
            nrOfVols = length(current_dicom_files);
            
            % identifying whether this is a regular experiment run or a run
            % with the control task
            if ~isempty(strfind(subfolders(d).name, control_run_dirname_part))  % This is a control task run
                stcprefix = [current_subj '_control'];
                fmr_output_filename = fullfile(output_dir,[current_subj '_control.fmr']);
            else    % This is a regular experiment run
                stcprefix = [current_subj '_' num2str(run_number_counter)];
                fmr_output_filename = fullfile(output_dir,[current_subj '_' num2str(run_number_counter) '.fmr']);
                run_number_counter = run_number_counter + 1;
            end
            
            % create the FMR file
            fmr = bvqx.CreateProjectMosaicFMR(fileType, firstFile,...
                nrOfVols, skipVols, createAMR, nrSlices,...
                stcprefix, byteswap, mosaicSizeX, mosaicSizeY, bytesperpixel,...
                output_dir, nrVolsInImg, sizeX, sizeY);
            fmr.SaveAs(fmr_output_filename);
            
        end
    end
end




%% PRE-PROCESSING OF FUNCTIONAL DATA

% Slice timing correction
disp('Slice timing correction')
for s=1:length(subject_names)
    disp(s);
    current_subj = subject_names(s).name;
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    
    % finding relevant fmr files
    fmr_files = dir(fullfile(output_dir, '*.fmr'));
    fmr_files = {fmr_files.name};
    for i=1:length(fmr_files)
        fmr_files{i} = fullfile(output_dir, fmr_files{i});
    end
    % removing FMR files which were created for the first volume only
    ix = 1:length(fmr_files);
    for i=1:length(ix)
        if strcmp(fmr_files{i}(end-11:end-4),'firstvol')
            ix(ix==i) = [];
        end
    end
    fmr_files = fmr_files(ix);
    
    % performing slice timing correction using the time table in the dicom
    for i=1:length(fmr_files)
        fmr = bvqx.OpenDocument(fmr_files{i});
        
        % defining slice timing correction parameters
        interpolation_type = 1;  % 1 - cubic spline interpolation
        
        % performing slice timing correction
        fmr.CorrectSliceTimingUsingTimeTable(interpolation_type);
        
        % deleting the original files without deleting firstvol_as_anat fmr file
        delete(fmr_files{i});
        delete([fmr_files{i}(1:end-3) 'stc']);
    end
end


% FMR motion correction
disp('FMR motion correction')
for s=1:length(subject_names)
    disp(s);
    current_subj = subject_names(s).name;
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    
    % finding relevant fmr files which have undergone slice timing correction
    fmr_files = dir(fullfile(output_dir, '*_SCCTBL.fmr'));
    fmr_files = {fmr_files.name};
    for i=1:length(fmr_files)
        fmr_files{i} = fullfile(output_dir, fmr_files{i});
    end
    
    % performing motion correction
    for i=1:length(fmr_files)
        fmr = bvqx.OpenDocument(fmr_files{i});
        
        % defining process parameters
        TargetVolume = 1;           % the number of volume in the series to align to
        Interpolation_type = 2;     % 0 and 1 - trilinear detection and trilinear interpolation, 2: trilinear detection and sinc interpolation or 3: sinc detection of motion and sinc interpolation
        UseFullDataSet = false;     % false is the default in the GUI
        MaxNumIterations = 100;     % 100 is the default in the GUI
        GenerateMovie = true;
        GenerateLogFile = true;     % creates a log file with the movement parameters
        
        % performing motion correction
        fmr.CorrectMotionEx(TargetVolume, Interpolation_type, UseFullDataSet,...
            MaxNumIterations, GenerateMovie, GenerateLogFile);
        fmr.Remove; % close or remove input FMR
    end
end



% High pass filtering
disp('High pass filtering')
for s=1:length(subject_names)
    disp(s);
    current_subj = subject_names(s).name;
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    
    % finding relevant fmr files
    fmr_files = dir(fullfile(output_dir,'*_3DMCTS.fmr'));   % files which have undergone motion correction
    fmr_files = {fmr_files.name};
    for i=1:length(fmr_files)
        fmr_files{i} = fullfile(output_dir, fmr_files{i});
    end
    
    % performing high pass filtering
    for i=1:length(fmr_files)
        fmr = bvqx.OpenDocument(fmr_files{i});
        
        % defining parameters for filtering
        NumCycles = 2;
        
        % performing filtering
        fmr.TemporalHighPassFilterGLMFourier(NumCycles);
        fmr.Remove; % close or remove input FMR
    end
end


% Spatial smoothing
disp('Smoothing')
for s=1:length(subject_names)
    disp(s);
    current_subj = subject_names(s).name;
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    
    % finding relevant fmr files
    fmr_files = dir(fullfile(output_dir,'*_THPGLMF2c.fmr'));    % files which have undergone high pass filtering
    fmr_files = {fmr_files.name}; 
    for i=1:length(fmr_files)
        fmr_files{i} = fullfile(output_dir, fmr_files{i}); 
    end
    
    for i=1:length(fmr_files)
        fmr = bvqx.OpenDocument(fmr_files{i});
        
        % defining parameters for smoothing
        FWHM = 4;   % smoothing by 4mm FWHM
        
        % performing spatial smoothing
        fmr.SpatialGaussianSmoothing(FWHM, 'mm');
    end
end



%% Anatomical data (VMR) preprocessing

% Inhomogeneity correction and skull stripping
disp('Inhomogeneity correction')
for s=1:length(subject_names)
    disp(s);
    current_subj = subject_names(s).name;
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    
    % identifying the VMR file
    vmr_file = dir(fullfile(output_dir, '*.vmr'));
    vmr_file = fullfile(output_dir, vmr_file(1).name);
    vmr = bvqx.OpenDocument(vmr_file);
    
    % correcting inhomogeneity
    vmr.CorrectIntensityInhomogeneities();
end


% Normalization to MNI space
disp('Normalizing')
for s=1:length(subject_names)
    disp(s);
    current_subj = subject_names(s).name;
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    
    vmr_file = dir(fullfile(output_dir, '*IIHC.vmr'));  % VMR file which has undergone inhomogeneity correction
    vmr_file = fullfile(output_dir, vmr_file(1).name);
    vmr = bvqx.OpenDocument(vmr_file);
    
    % performing spatial normalization
    vmr.NormalizeToMNISpace();
end


%% IMPOTANT - MANUALLY CHECK NORMALIZATION RESULTS BEFORE CLOSING THE WINDOW TO VERIFY NORMALIZATION (USE F8 TO SEE OVERLAY)



%% Co-registration of functional and anatomical data

disp('coregistering...')
for s=1:length(subject_names)
    disp(s);
    current_subj = subject_names(s).name;
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    
    % finding the relevant fmr files
    fmr_files = dir(fullfile(output_dir, '*SD3DSS*mm.fmr'));    % functional files which have been spatially smoothed
    vmr_file = dir(fullfile(output_dir,'*IIHC.vmr'));           % anatomical file which has been corrected for inhomogeneity, BUT NOT YET NORMALIZED TO MNI SPACE (IN NATIVE SPACE)
    vmr_file = fullfile(output_dir, vmr_file(1).name);
    
    % Coregistration using the intensity gradient-based matching algorithm
    useAttachedAMR = 1;
    for i=1:length(fmr_files)
        vmr = bvqx.OpenDocument(vmr_file);
        fmr_filename = fullfile(output_dir, fmr_files(i).name);
        vmr.CoregisterFMRToVMR(fmr_filename, useAttachedAMR);
    end
    
end



%% IMPORTANT - MANUALLY CHECK COREGISTRATION RESULTS BEFORE CLOSING THE WINDOWS, AND CORRECT MANUALLY IN BRAINVOYAGER IF THEY ARE PROBLEMATIC



%% Creating VTC files (files with functional data normalized to MNI space)

disp('Creating VTCs...')
for s=1:length(subject_names)
    current_subj = subject_names(s).name;
    disp(s);
    output_dir = fullfile(subjects_output_data_dir, current_subj);
    
    % Identifying the normalized anatomical data file
    vmr_file = dir(fullfile(output_dir, '*MNI.vmr'));
    vmr_file = fullfile(output_dir, vmr_file(1).name);
    vmr = bvqx.OpenDocument(vmr_file);
    
    % identifying relevant files 
    fmr_files = dir(fullfile(output_dir, '*SD3DSS*mm.fmr'));    % all of the spatially smoothed FMR files
    IA_files = dir(fullfile(output_dir, '*_IIHC_IA.trf'));      % files with parameters of co-registration (output of coregistration step)
    FA_files = dir(fullfile(output_dir, '*_IIHC*_FA.trf'));     % files with parameters of co-registration (output of coregistration step)
    MNI_transform_file = dir(fullfile(output_dir, '*MNI_a12.trf')); % file with the parameters of transformation to MNI space (output of VMR spatial normalization)
    MNI_transform_file = fullfile(output_dir, MNI_transform_file(1).name);
    
    for i=1:length(fmr_files)
        fmr_filename = fullfile(output_dir, fmr_files(i).name);
        IA_filename = fullfile(output_dir, IA_files(i).name);
        FA_filename = fullfile(output_dir, FA_files(i).name);
        
        % defining parameters of VTC creation
        Datatype = 2;       % 1 - integer, 2 - float (GUI default)
        Resolution = 3;     % resolution relative to VMR - 3x3x3
        Interpolation = 1;    % 0 for nearest neighbor , 1 for trilinear (GUI default) , 2 for sinc interpolation
        Intensity_threshold = 100;    % minimum intensity threshold for voxels
        output_VTC_name = [fmr_filename(1:end-4) '_MNI.vtc'];   % name of output VTC file
        
        % creating VTC file
        vmr.CreateVTCInMNISpace(fmr_filename, IA_filename, FA_filename, MNI_transform_file, ...
            output_VTC_name, Datatype, Resolution, Interpolation, Intensity_threshold);
    end
end



%% Manual stage: verify functional coverage for all VTCs (can be done later using the MDM files)



