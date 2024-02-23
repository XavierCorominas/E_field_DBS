
%% *************************** Tutorial DBS_E_FIELD 02/2024 ********************************

% Xavier Corominas Teruel, Paris France.
% The tutorial below outlines a systematic approach for reconstructing the spatial coordinates and trajectory 
% of DBS electrodes implanted in human brains. This process enables the simulation of E-fields with external electrical 
% currents and the estimation of their impact on the electrodes. There are several ways to do so, this is just one of them. 
% To execute this pipeline successfully, ensure the following software is installeE:
% -	Simnibs 4.0 (SimNIBS 4 — SimNIBS g75fd674e') documentation).
% -	SPM12 (SPM12 Software - Statistical Parametric Mapping (ucl.ac.uk))
% -	CAT12 (CAT12 Manual (neuro-jena.github.io))
% -	ITK snap (ITK-SNAP Home (itksnap.org)).
% -	MRIcron (NITRC: MRIcron: Tool/Resource Info) or any other visualized like FSL, freesurfer, etc.
% -	Matlab 2022a or later (or Python if preferred).

% As a general summary, the following precess encompasses the following steps: 
    % (I) MRI head mesh reconstruction), 
    % (II) DBS electrode and mask delineation, 
    % (III) updating head models, 
    % (IV) E-field simulation,      
    % (V) analyses of E-field impact.
% Let’s kick in!


% For the head models, it is highly recommended to complement the head model with T2 or flair images when possible.

%% Step 0: Add paths

clear;
clc;
close all;

% Set up the paths to needed tools --> Modify this according to your location:
% Simnibs :
addpath('C:\Users\XAVIER\SimNIBS-4.0\simnibs_env\Lib\site-packages\simnibs\matlab_tools\')

%Presurfer toolbox for MRI denoising (if necessary, % https://github.com/srikash/presurfer );
addpath('E:\LI_rTMS_project\eeg_analysis\Analysis\External_functions\presurfer-main\func')



% Set subject ID (example):

subject_id = 'test';
sequence = 'S4'

%% Step 1: Generate Head models

%{
       % In case is necessary to denoise the image froma MP2rage T1 sequence

    path_to_INV2 = ['E:\MRI\',subject_id,'\MRI\v_CRETMS_01_001_S5_T1w_mp2rage_INV2.nii'];
    path_to_UNI = ['E:\MRI\',subject_id,'\MRI\v_CRETMS_01_001_S4_T1w_mp2rage_UNI_Images.nii'];
    presurf_MPRAGEise(path_to_INV2,path_to_UNI);

% Alternatively you can manually denoise the T1 with Benoit's Beranger toolbox runing above SPM12 (Cenir ICM engineer) 
% Download https://github.com/benoitberanger/mp2rage
% This works with upon the graphycal user interface of SPM12

%}

%Now, run the SIMNIBS charm protocol to segment the head tissues. Please read the SIMNIBS website for further information:
% Set paths to the MRI that you want to reconstruct in a head model. Change the paths according to your needs.

    path_to_mri_nifti = ['E:\MRI\',subject_id,'\MRI'];
    filename = ['clean_v_CRETMS_01_',subject_id,'_',sequence,'_T1w_mp2rage_UNI_Images.nii'];

cd(path_to_mri_nifti)
clc;
[status,cmdout] = system(['cd ',path_to_mri_nifti],'-echo');
[status,cmdout] = system(['C:\Users\XAVIER\SimNIBS-4.0\bin\charm',' ',subject_id, ' ' ,filename, ' --forceqform --forcerun'], '-echo');

% a folder called m2m_subjectid will be created containing the headmodel


%% Step 2: Manually dealineate the mask and update layers of the head model. 
% This step is not fully done automatically by a script, needs to be done manually.


%2. Here we will do the segmentation manually of the area
%(lesion,electrode...) in to the original MRI with ITK-snap (http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.SNAP3, or similar
% software). Important to do it from the same MRI where we do the
% headmodel. Then the new segmentation layer will be stored in nifti
% format(.nii) and stored in the m2m_subjectID folder of the headmodel (e.g., E:\MRI\test\MRI\m2m_test\)

%Important
% The mask will be stored uncompresses (.nii) with the name: electrode
% The mask needs to be stored in the m2m_subjid forlder

% Then we will actualize the headmodel with the new tissue layer. 
% Here, the new tissue (that we call "lesioned_tissue.nii") will be added with the label 35 : to visualize it later on in the Gmsh software you have to activate the Tools--> Options-->Geometry-->Visibiilty-->Volumes to see the mask)
% Go to the m2m_sibjectID folder and run add_tissues_to_upsampled

cd(path_to_mri_nifti)
% addpath(['m2m_',subject_id]);
clc;
[status,cmdout] = system(['cd ',path_to_mri_nifti,'\m2m_',subject_id,'\'],'-echo');
%[status,cmdout] = system(['dir ',path_to_mri_nifti, '\m2m_',subject_id,'\ -echo']);
[status,cmdout] = system(['C:\Users\XAVIER\SimNIBS-4.0\bin\add_tissues_to_upsampled',' ','-i',' ','m2m_',subject_id,'\electrode.nii',' ','-t',' ' ,'m2m_',subject_id,'\label_prep\tissue_labeling_upsampled.nii.gz',' ','-o',' ','m2m_',subject_id,'\label_prep\tissue_labeling_upsampled.nii.gz',' ','--offset 35'], '-echo');


[status,cmdout] = system(['cd ',path_to_mri_nifti],'-echo');
[status,cmdout] = system(['C:\Users\XAVIER\SimNIBS-4.0\bin\charm',' ',subject_id,' ','--mesh'], '-echo');
% 


%% SIMULATION tDCS cerebellum ( anode right cerebellum, cathode left cerebellum)
% Create Save folder
save_path = ([path_to_mri_nifti ,'\m2m_',subject_id,'\']);

save_folder = [save_path,'simulationtDCScer\']
if ~exist(save_folder)
      mkdir([save_folder]);
end

% RUN SIMULATION
% Initialize  simulation session
s = sim_struct('SESSION');
% path to mesh
s.fnamehead = ([path_to_mri_nifti ,'\m2m_',subject_id,'\',subject_id,'.msh']);
% Output folder
s.pathfem = ([path_to_mri_nifti ,'\m2m_',subject_id,'\simulationtDCScer\']);
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
s.poslist{1}.currents = [2e-3, -2e-3]; % 1e-3=1mA current / 2e-3=2mA current. Here we use one electrode as the anodal and a second as the cathodal
%s.poslist{1}.currents = [2e-3, -0.5e-3, -0.5e-3, -0.5e-3, -0.5e-3]; % 1e-3=1mA current / 2e-3=2mA current. Here we use one electrode as the anodal and a second as the cathodal

% set conductivities for the new tissue 35 (the electrode); THE REST OF THE
% CONDUCTIVITIES ARE SET TO DEFAULT
s.poslist{1}.cond(36).value = 29.4; % in S/m . We assume that has the same conductivity of a silicon rubber (that is much higher than any brain tissue)/ https://simnibs.github.io/simnibs/build/html/documentation/conductivity.html
s.poslist{1}.cond(36).name = 'electrodeDBS'; 

%First electrode
% 3 elements: Electrode with three layers, corresponding silicone rubber electrode inside a sponge with saline solution. 
% Only valid for rectangular and elliptical electrodes
s.poslist{1}.electrode(1).channelnr = 1;
s.poslist{1}.electrode(1).shape = 'rect' % or .ellipse’
s.poslist{1}.electrode(1).dimensions = [45, 65] % Simulate a 45x65mm rubber electrode in a 50X70mm sponge
s.poslist{1}.electrode(1).dimensions_sponge = [50, 70]
s.poslist{1}.electrode(1).thickness = [4, 2, 4] % Electrode 2mm thick rubber electrode in the middle of a 8mm thick sponge
s.poslist{1}.electrode(1).centre = 'PO10'; % Place it over C3
s.poslist{1}.electrode(1).pos_ydir = 'PO9'; % y axis of electrode towards contralateral cerebellum
%{
% Alternatively use MNI coordinates to position the center of theelectrodes. In that case we use the skin-projected coords of crus2;
s.poslist{1}.electrode(1).centre = mni2subject_coords([23, -65, -47], 'm2m_subjectid');
s.poslist{1}.electrode(1).pos_ydir = mni2subject_coords([-23, -65, -47], 'm2m_subjectid');
%}

%Second electrode
s.poslist{1}.electrode(2).channelnr = 2;
s.poslist{1}.electrode(2).shape = 'rect' % or .ellipse2
s.poslist{1}.electrode(2).dimensions = [45, 65] % Simulate a 45x65mm rubber electrode in a 50X70mm sponge
s.poslist{1}.electrode(2).dimensions_sponge = [50, 70]
s.poslist{1}.electrode(2).thickness = [4, 2, 4] % Electrode 2mm thick rubber electrode in the middle of a 8mm thick sponge
s.poslist{1}.electrode(2).centre = 'PO9'; % Place it over C3
s.poslist{1}.electrode(2).pos_ydir = 'PO10'; % y axis of electrode towards contralateral cerebellum
%{
% Alternatively use MNI coordinates to position the center of theelectrodes. In that case we use the skin-projected coords of crus2;
s.poslist{1}.electrode(2).centre = mni2subject_coords([-23, -65, -47], 'm2m_subjectid');
s.poslist{1}.electrode(2).pos_ydir = mni2subject_coords([23, -65, -47], 'm2m_subjectid');
%}

% Run Simulation
run_simnibs(s);




%% Visualize Simulations and analyze
%load headmesh corrected (with the electrode)
head = mesh_load_gmsh4(['E:\MRI\',subject_id,'\MRI\m2m_',subject_id,'\simulationtDCScer\',subject_id,'_tDCS_1_scalar.msh']);

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




%% Simple roi analyses from a MNI coords

head = mesh_load_gmsh4(['E:\MRI\',subject_id,'\MRI\m2m_',subject_id,'\simulationtDCScer\',subject_id,'_tDCS_1_scalar.msh']);

% Reference values of Efield impact to the electrode.
electrode = mesh_extract_regions(head, 'region_idx', 36);

DBSResults.E_field_99pct = prctile(electrode.element_data{2,1}.tetdata,99);
DBSResults.E_field_50pct = prctile(electrode.element_data{2,1}.tetdata,50);
DBSResults.E_field_mean = mean(electrode.element_data{2,1}.tetdata);
DBSResults.E_field_min = min(electrode.element_data{2,1}.tetdata);
DBSResults.E_field_max = max(electrode.element_data{2,1}.tetdata);
DBSResults

faces= electrode.tetrahedra;
vertices= electrode.nodes;
values = electrode.element_data{2, 1}.tetdata ; % Sample values, replace with your actual data

% Plot each face separately and assign colors based on values
for i = 1:size(faces, 1)
    face_vertices = vertices(faces(i, :), :);
    %patch(face_vertices(:, 1), face_vertices(:, 2), face_vertices(:, 3), values(i),'EdgeAlpha',0.5, 'FaceColor', 'flat');
     patch(face_vertices(:, 1), face_vertices(:, 2), face_vertices(:, 3), values(i),'EdgeAlpha',0.2);
    hold on;
end
% lightangle(55,-35)
% h.FaceLighting = 'gouraud';
% h.AmbientStrength = 0.7;
% h.DiffuseStrength = 0.4;
% h.SpecularStrength = 0.9;
% h.SpecularExponent = 25;
% h.BackFaceLighting = 'lit';
% Add a colorbar to represent the values
c= colorbar;
c.Label.String = 'E-field strength (V/m)'; % Set the new legend text
c.FontSize = 15;
axis off



