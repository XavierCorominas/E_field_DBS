%% *************************** Tutorial tDCS electrode optimization 03/2024 ********************************

% Xavier Corominas Teruel, Paris France.
% The tutorial below outlines a systematic approach for tDCS electrode
% optimization.

% All optimizations are performed in the MNI152 template brain.

% Important: if you try to create a leadfield from a headmesh that has been
% previously updated with the electrode mask, this will not work.

%% 1. Compute leadfield with simnibs

% Copy paste head mesh folder from original location to every sub-subject folder.
% Modify paths in your computer
% Source folder
source_folder = (['D:\SOFTWARE\spm12\spm12\toolbox\cat12\templates_MNI152NLin2009cAsym\m2m_MNI152']); % Replace this with the path to your source folder
% Destination folder

subject_id = 'test';
destination_folder =(['D:\MRI\',subject_id]); 

save_folder =(['D:\MRI\',subject_id,'\m2m_MNI152_copy\']); 
if ~exist(save_folder)
      mkdir([save_folder]);
end


% Copy the contents of the source folder to the destination folder
try
    copyfile(source_folder, save_folder);
    disp('Folder copied successfully.');
catch ME
    fprintf('An error occurred: %s\n', ME.message);
end
%%

% Set subject
    subject_id = 'test';
    path_to_mri_nifti = ['D:\MRI\',subject_id];

%create leadfield folder
save_path = ([path_to_mri_nifti,'\m2m_MNI152_copy\']);

save_folder = [save_path,'tdcs_leadfield\']
if ~exist(save_folder)
      mkdir([save_folder]);
end


% place script in the main folder of the example dataset
tdcs_lf = sim_struct('TDCSLEADFIELD');
% subject folder
%tdcs_lf.subpath = 'm2m_MNI152';
tdcs_lf.fnamehead = (['D:\MRI\',subject_id,'\m2m_MNI152_copy\MNI152.msh'])
% Output directory
tdcs_lf.pathfem = [save_path,'\tdcs_leadfield\'];

%Default eeg cap is 10-10 jurak. in case you want to use neuroelectrics cap:
% tdcs_lf.eeg_cap = (['D:\MRI\',subject_id,'\m2m_MNI152\eeg_positions\EEG10-10_Neuroelectrics.csv']);
% Uncomment to use the pardiso solver
tdcs_lf.solver_options = 'pardiso';
% This solver is much faster than the default. However, it requires much more memory (~12 GB)

run_simnibs(tdcs_lf)

%%  Optimization:

% 
cd(['D:\MRI\',subject_id,'\m2m_MNI152_copy\tdcs_leadfield\']);

% Initialize structure
opt = opt_struct('TDCSoptimize');

% Select the leadfield file
opt.leadfield_hdf = (['D:\MRI\',subject_id,'\m2m_MNI152_copy\tdcs_leadfield\MNI152_leadfield_EEG10-10_UI_Jurak_2007.hdf5']);

% Select a name for the optimization
opt.name = 'optimization/single_target_cer';

% Select a maximum total current (in A)
opt.max_total_current = 2e-3;
% Select a maximum current at each electrodes (in A)
opt.max_individual_current = 2e-3;
% Select a maximum number of active electrodes (optional)
opt.max_active_electrodes = 5;

% Define optimization target
% Position of target, in subject space!
% please see tdcs_optimize_mni.m for how to use MNI coordinates
opt.target.positions = [24, -66, -40];
% Intensity of the electric field (in V/m)
opt.target.intensity = 0.25;
% optimize without direction (understanding the e-fild as a mass without
% directions)
opt.target.directions='none';

% Run optimization
run_simnibs(opt);

%%


