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
opt.target.positions = [17, -65, -28];
% Intensity of the electric field (in V/m)
opt.target.intensity = 0.25;
% optimize without direction (understanding the e-fild as a mass without
% directions)
opt.target.directions='none';
opt.avoid.positions = [50, -13, -34]; % tenporal lobe
opt.avoid.positions = [50, -18, -34]; % tenporal lobe

%opt.avoid.positions = [-9, -88, -44]; % cotnralateral cerebell


% Run optimization
run_simnibs(opt);

%% SIMULATION ACCORDING TO OPTIMIZED PARAMETERS

path_to_subject_MNI152_nifti =(['D:\MRI\',subject_id]);

save_path = ([path_to_subject_MNI152_nifti,'\m2m_MNI152_copy\tdcs_leadfield\optimization\']);

save_folder = [save_path,'simulation\']
if ~exist(save_folder)
      mkdir([save_folder]);
end

% RUN SIMULATION
% Initialize  simulation session
s = sim_struct('SESSION');
% path to mesh
s.fnamehead = ([path_to_subject_MNI152_nifti ,'\m2m_MNI152_copy\MNI152.msh']);
% Output folder
s.pathfem = ([path_to_subject_MNI152_nifti,'\m2m_MNI152_copy\tdcs_leadfield\optimization\simulation']);
%EEG cap according to neuroelectrics
%s.eeg_cap = ([path_to_mri_nifti,'\m2m_',subject_id,'\eeg_positions\EEG10-10_Neuroelectrics.csv']); 

% Fields to maps
s.fields = 'eEjJ'; % Save the following results:
                   %  e: Electric field magnitude
                   %  E: Electric field vector
                   %  j: Current density magnitude
                   %  J: Current density vector
% set different output maps                   
s.map_to_surf = true;   %  Map to subject's middle gray matter surface
s.map_to_fsavg = true;  %  Map to FreeSurfer's FSAverage group template
s.map_to_vol = true;    %  Save as nifti volume
s.map_to_MNI = true;    %  Save in MNI space

% Get fields everywhere: 
s.tissues_in_niftis = 'all'
            %s.tissues_in_niftis = [1,2,3]; % Results in the niftis will be masked 
                                            % to only show WM (1), GM (2), CSF(3)
                                            % (standarE: only GM)


%Run init-structure
s.poslist{1} = sim_struct('TDCSLIST');
%Set electrodes (how many and intensities)
s.poslist{1}.currents = [-0.00038, -0.001, -0.00035, 0.002, -0.00027]; % 1e-3=1mA current / 2e-3=2mA current. Here we use one electrode as the anodal and a second as the cathodal
%s.poslist{1}.currents = [2e-3, -0.5e-3, -0.5e-3, -0.5e-3, -0.5e-3]; % 1e-3=1mA current / 2e-3=2mA current. Here we use one electrode as the anodal and a second as the cathodal

% set conductivities for the new tissue 35 (the electrode); THE REST OF THE
% CONDUCTIVITIES ARE SET TO DEFAULT
%s.poslist{1}.cond(36).value = 29.4; % in S/m . We assume that has the same conductivity of a silicon rubber (that is much higher than any brain tissue)/ https://simnibs.github.io/simnibs/build/html/documentation/conductivity.html
%s.poslist{1}.cond(36).name = 'electrodeDBS'; 

%First electrode
% Only valid for rectangular and elliptical electrodes
s.poslist{1}.electrode(1).channelnr = 1;
s.poslist{1}.electrode(1).shape = 'ellipse' % or .ellipse’
s.poslist{1}.electrode(1).dimensions = [10, 10] % Simulate a 10X10  electrode 
s.poslist{1}.electrode(1).thickness = [2, 2] % Electrode 2mm thick  electrode in top od 2 mm gel
s.poslist{1}.electrode(1).centre = 'C6'; % P

%Second electrode
% Only valid for rectangular and elliptical electrodes
s.poslist{1}.electrode(2).channelnr = 2;
s.poslist{1}.electrode(2).shape = 'ellipse' % or .ellipse2
s.poslist{1}.electrode(2).dimensions = [10, 10] % Simulate a 10X10  electrode 
s.poslist{1}.electrode(2).thickness = [2, 2] % Electrode 2mm thick  electrode in top od 2 mm gel
s.poslist{1}.electrode(2).centre = 'P9'; % P

%Third electrode
% Only valid for rectangular and elliptical electrodes
s.poslist{1}.electrode(3).channelnr = 3;
s.poslist{1}.electrode(3).shape = 'ellipse' % or .ellipse2
s.poslist{1}.electrode(3).dimensions = [10, 10] % Simulate a 10X10  electrode 
s.poslist{1}.electrode(3).thickness = [2, 2] % Electrode 2mm thick  electrode in top od 2 mm gel
s.poslist{1}.electrode(3).centre = 'PO4'; % P

%Fourth electrode
% Only valid for rectangular and elliptical electrodes
s.poslist{1}.electrode(4).channelnr = 4;
s.poslist{1}.electrode(4).shape = 'ellipse' % or .ellipse2
s.poslist{1}.electrode(4).dimensions = [10, 10] % Simulate a 10X10  electrode 
s.poslist{1}.electrode(4).thickness = [2, 2] % Electrode 2mm thick  electrode in top od 2 mm gel
s.poslist{1}.electrode(4).centre = 'PO10'; % P

%Fifth electrode
% Only valid for rectangular and elliptical electrodes
s.poslist{1}.electrode(5).channelnr = 5;
s.poslist{1}.electrode(5).shape = 'ellipse' % or .ellipse2
s.poslist{1}.electrode(5).dimensions = [10, 10] % Simulate a 10X10  electrode 
s.poslist{1}.electrode(5).thickness = [2, 2] % Electrode 2mm thick  electrode in top od 2 mm gel
s.poslist{1}.electrode(5).centre = 'Oz'; % P

% Run Simulation
run_simnibs(s);

%%


%% Visualize Simulations and analyze
%load headmesh corrected (with the electrode)

subject_id = (['test'])
head = mesh_load_gmsh4(['D:\MRI\',subject_id,'\m2m_MNI152_copy\tdcs_leadfield\optimization\simulation\MNI152_TDCS_1_scalar.msh']);

%GREY MATTER PLOT and analyses
% Isolate electrode mask
electrode.e =  head.triangles(head.triangle_regions == 1002,:);
electrode.p =  head.nodes;

% Make shure head.elements_data {2,1} corresponds to norme_E
E_field = head.element_data{2,1}.tridata(head.triangle_regions==1002);

% Reference values of Efield impact to the electrode.
Results.E_field_99pct = prctile(E_field,99);
Results.E_field_50pct = prctile(E_field,50);
Results.E_field_mean = mean(E_field);
Results
%visualize electrode alone with E-fields
figure;
trisurf(electrode.e,electrode.p(:,1),electrode.p(:,2),electrode.p(:,3),E_field,'EdgeAlpha',0.5);
colorbar
% Make shure head.elements_data {2,1} corresponds to norme_E
t = E_field > 0.8*prctile(E_field,99); % THIS A CRITICAL POINT, 0.1 is a arbitrary value to select the extension of the passband filter, MODIFY IT TO YOUR NEEDS
figure;
trisurf(electrode.e,electrode.p(:,1),electrode.p(:,2),electrode.p(:,3),double(t),'EdgeAlpha',0.5);
colorbar

% mesh_get_simulation_result
% mesh_get_simulation_result_dbselectrode


% WHite matter plot and analyses
%GREY MATTER PLOT and analyses
% Isolate electrode mask
electrode.e =  head.triangles(head.triangle_regions == 1001,:);
electrode.p =  head.nodes;

% Make shure head.elements_data {2,1} corresponds to norme_E
E_field = head.element_data{2,1}.tridata(head.triangle_regions==1001);

% Reference values of Efield impact to the electrode.
Results.E_field_99pct = prctile(E_field,99);
Results.E_field_50pct = prctile(E_field,50);
Results.E_field_mean = mean(E_field);
Results
%visualize electrode alone with E-fields
figure;
trisurf(electrode.e,electrode.p(:,1),electrode.p(:,2),electrode.p(:,3),E_field,'EdgeAlpha',0.5);
colorbar
% Make shure head.elements_data {2,1} corresponds to norme_E
t = E_field > 0.8*prctile(E_field,99); % THIS A CRITICAL POINT, 0.1 is a arbitrary value to select the extension of the passband filter, MODIFY IT TO YOUR NEEDS
figure;
trisurf(electrode.e,electrode.p(:,1),electrode.p(:,2),electrode.p(:,3),double(t),'EdgeAlpha',0.5);
colorbar
