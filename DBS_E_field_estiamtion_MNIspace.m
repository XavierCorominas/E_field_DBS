%% *************************** Tutorial DBS_E_FIELD 03/2024 ********************************

% Xavier Corominas Teruel, Paris France.

% Tutorial for MNI normalized space.

% The tutorial below outlines a systematic approach for reconstructing the spatial coordinates and trajectory 
% of DBS electrodes implanted in human brains in MNI space. This process enables the simulation of E-fields with external electrical 
% currents and the estimation of their impact on the electrodes. There are several ways to do so, this is just one of them. 
% To execute this pipeline successfully, ensure the following software is installeE:
% -	Simnibs 4.0 (SimNIBS 4 — SimNIBS g75fd674e') documentation).
% -	SPM12 (SPM12 Software - Statistical Parametric Mapping (ucl.ac.uk))
% -	CAT12 (CAT12 Manual (neuro-jena.github.io))
% -	ITK snap (ITK-SNAP Home (itksnap.org)).
% -	MRIcron (NITRC: MRIcron: Tool/Resource Info) or any other visualized like FSL, freesurfer, etc.
% -	Matlab 2022a or later (or Python if preferred).

% As a general summary, the following precess encompasses the following steps: 
    % (I) Native space MRI transformation to MNI space, 
    % (II) DBS electrode and mask delineation (manually with ITK snap), 
    % (III) Updating head MNI model with electrode mask mesh, 
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




%% Step 1: Native space MRI transformation to MNI space

%{
       % In case is necessary to denoise the image froma MP2rage T1 sequence

    path_to_INV2 = ['E:\MRI\',subject_id,'\MRI\v_CRETMS_01_001_S5_T1w_mp2rage_INV2.nii'];
    path_to_UNI = ['E:\MRI\',subject_id,'\MRI\v_CRETMS_01_001_S4_T1w_mp2rage_UNI_Images.nii'];
    presurf_MPRAGEise(path_to_INV2,path_to_UNI);

% Alternatively you can manually denoise the T1 with Benoit's Beranger toolbox runing above SPM12 (Cenir ICM engineer) 
% Download https://github.com/benoitberanger/mp2rage
% This works with upon the graphycal user interface of SPM12

%}

% Now transform to MNi space the subject space MRI

     %   To do so, first we have to segment the subject (step 1.1) space MRI and then
     %   Normalize (step 1.2)

     % Step 1.1 - Segmentation (this will take 20 mins approx)
% List of open inputs
% ------> Set this variables for your subject and computer <--------------------------------
MRI_path = (['D:\MRI\test\MRI\clean_v_CRETMS_01_test_S4_T1w_mp2rage_UNI_Images.nii']); % path to your MRI subject space
TPM_path = (['D:\SOFTWARE\spm12\spm12\tpm\TPM.nii']); % this is the default path after spm installation in my comptuer
Template_path = (['D:\SOFTWARE\spm12\spm12\toolbox\cat12\templates_MNI152NLin2009cAsym\Template_0_GS.nii']); % this is the default location of the template in the cat12

% Job CAT12 segmentation
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
% segmentation function
matlabbatch{1}.spm.tools.cat.estwrite.data = {[MRI_path,',1']};
matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 2;
matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {[TPM_path]};
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasacc = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.optimal = [1 0.3];
matlabbatch{1}.spm.tools.cat.estwrite.extopts.setCOM = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.affmod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.spm_kamap = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASmyostr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.WMHC = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {[Template_path]};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.bb = 12;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.SRP = 22;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.ignoreErrors = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1];
matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 0;

% run segmentation job 
nrun = 1; % enter the number of runs here
jobs = repmat(matlabbatch, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});


% --------> A subfolder called mri\ containing the y / iy transformation will be automatially generated
%%  Step 1.2 - Normalization (this will take a few seconds)

% List of open inputs
% ------> Set this variables for your subject and computer <--------------------------------
y_dimension_mri = (['D:\MRI\test\MRI\mri\y_clean_v_CRETMS_01_test_S4_T1w_mp2rage_UNI_Images.nii']);
% MRI_path = image to normalize is the original MRI previously defined

% Job CAT12 normalization
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
% normalization function

matlabbatchnorm{1}.spm.tools.cat.tools.defs.field1 = {[y_dimension_mri,',1']};
matlabbatchnorm{1}.spm.tools.cat.tools.defs.images = {[MRI_path,',1']};
matlabbatchnorm{1}.spm.tools.cat.tools.defs.bb = [NaN NaN NaN
                                              NaN NaN NaN];
matlabbatchnorm{1}.spm.tools.cat.tools.defs.vox = [NaN NaN NaN];
matlabbatchnorm{1}.spm.tools.cat.tools.defs.interp = 1;
matlabbatchnorm{1}.spm.tools.cat.tools.defs.modulate = 0;

% run normalization job 
nrun = 1; % enter the number of runs here
jobs = repmat(matlabbatchnorm, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

% ----------> a subfile .nii will be automatically generated with the MRI name preceeded witha  w
 % In this example my file is called: 
            % wclean_v_CRETMS_01_test_S4_T1w_mp2rage_UNI_Images

            % Please open the file in MRIcron or any other visualizer to
            % check the quality of the MNI normalization


%% Step 2: Manually dealineate the mask and update layers of the MNI default head model. 

% This step is not fully done automatically by a script, needs to be done manually.

% Then we will create a MNI152 head model from the template MRI we used for the SPM
% normalization step. A subfolder with a MNI152 head model needs to exist
% for every subject

%Now, run the SIMNIBS charm protocol to segment the head tissues. Please read the SIMNIBS website for further information:
% Set paths to the MRI that you want to reconstruct in a head model. Change the paths according to your needs.

% We will do the MNI152 head model once and copy paste inside every subejct
% folder.

path_template_MNI = (['D:\SOFTWARE\spm12\spm12\toolbox\cat12\templates_MNI152NLin2009cAsym']);
filename = ['Template_T1.nii'];

try 
    cd(path_template_MNI)
    clc;

[status,cmdout] = system(['cd ',path_template_MNI],'-echo');
[status,cmdout] = system(['C:\Users\XAVIER\SimNIBS-4.0\bin\charm',' ','MNI152', ' ' ,filename, ' --forceqform --forcerun'], '-echo');

catch ME
    fprintf('An error occurred: %s\n', ME.message);
end


% For tutorial purposes we will create for every of our subejcts a sub
% MNI152 folder and actualize the pertinent MNI152 mesh for every subejct. 
% (The other option for final visualization purposes will be to merge all electrodes in one
% mesh solution. The present code does snot do that).

% Copy paste head mesh folder from original location to every sub-subject folder.
% Modify paths in your computer
% Source folder
source_folder = (['D:\SOFTWARE\spm12\spm12\toolbox\cat12\templates_MNI152NLin2009cAsym\m2m_MNI152']); % Replace this with the path to your source folder
% Destination folder

subject_id = 'test';
destination_folder =(['D:\MRI\',subject_id]); 

save_folder =(['D:\MRI\',subject_id,'\m2m_MNI152\']); 
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
% Here we will do the segmentation manually of the area
%(lesion,electrode...) in to the MNI normalized MRI (wclean_v_CRETMS_01_test_S4_T1w_mp2rage_UNI_Images) 
%  with ITK-snap (http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.SNAP3, or similar
% software). 
% a folder called m2m_subjectid will be created containing the headmodel

% Then the new segmentation layer will be stored in nifti
% format(.nii) and stored in the subfolder m2m_MNI152 folder of the subject

%Important <----------------
% The mask will be stored uncompresses (.nii) with the name: electrode
% The mask needs to be stored in the m2m_MNI152 forlder
% (e.g., D:\MRI\test\MNI152\m2m_MNI152)

% Here, the new tissue (that we call "electrode.nii") will be added with the label 35 : to visualize it later on in the Gmsh software you have to activate the Tools--> Options-->Geometry-->Visibiilty-->Volumes to see the mask)
% Go to the \subjID\m2m_MNI152 folder and run add_tissues_to_upsampled

% Set subject ID (example):

subject_id = 'test';
path_to_subject_MNI152_nifti =(['D:\MRI\',subject_id]);
subject_id_MNI = (['MNI152']);

cd(path_to_subject_MNI152_nifti);

% addpath(['m2m_',subject_id]);
clc;
[status,cmdout] = system(['cd ',path_to_subject_MNI152_nifti],'-echo');
%[status,cmdout] = system(['dir ',path_to_mri_nifti, '\m2m_',subject_id,'\ -echo']);
[status,cmdout] = system(['C:\Users\XAVIER\SimNIBS-4.0\bin\add_tissues_to_upsampled',' ','-i',' ','m2m_',subject_id_MNI,'\electrode.nii',' ','-t',' ' ,'m2m_',subject_id_MNI,'\label_prep\tissue_labeling_upsampled.nii.gz',' ','-o',' ','m2m_',subject_id_MNI,'\label_prep\tissue_labeling_upsampled.nii.gz',' ','--offset 35'], '-echo');

[status,cmdout] = system(['cd ',path_to_subject_MNI152_nifti],'-echo');
[status,cmdout] = system(['C:\Users\XAVIER\SimNIBS-4.0\bin\charm',' ',subject_id_MNI,' ','--mesh'], '-echo');
% 



%% SIMULATION tDCS cerebellum (anode right cerebellum, cathode left cerebellum)
% Create Save folder

path_to_subject_MNI152_nifti =(['D:\MRI\',subject_id]);

save_path = ([path_to_subject_MNI152_nifti,'\m2m_MNI152\']);

save_folder = [save_path,'simulationtDCScer\']
if ~exist(save_folder)
      mkdir([save_folder]);
end

% RUN SIMULATION
% Initialize  simulation session
s = sim_struct('SESSION');
% path to mesh
s.fnamehead = ([path_to_subject_MNI152_nifti ,'\m2m_MNI152\MNI152.msh']);
% Output folder
s.pathfem = ([path_to_subject_MNI152_nifti,'\m2m_MNI152\simulationtDCScer\']);
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
head = mesh_load_gmsh4(['D:\MRI\',subject_id,'\m2m_MNI152\simulationtDCScer\MNI152_TDCS_1_scalar.msh']);

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




%% Simple  analyses of electrode mask

head = mesh_load_gmsh4(['D:\MRI\',subject_id,'\m2m_MNI152\simulationtDCScer\MNI152_TDCS_1_scalar.msh']);

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


