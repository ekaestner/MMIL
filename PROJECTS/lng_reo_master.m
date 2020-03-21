clear; clc;

%  Overall Data Locations
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'LanguageReorganization';

% Data Location
frs_loc = '/home/mmilmcd/data/FSRECONS/';
bld_loc = '/home/mmilmcd/data/MCD_BOLD/subjects';

% File Names
red_fle = 'sbj000_total_2019_03_27.csv';
fsr_sbj_nme = 'mmilmcdRSI_freesurfer_names.csv';
bld_sbj_nme = 'mmilmcdRSI_bold_names.csv';

grp_fle = 'Lng_Reo_subject_list.csv';

%


%% Put together Subjects
fsr_sbj_nme = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' fsr_sbj_nme]);
bld_sbj_nme = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' bld_sbj_nme]);

grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'projects' '/' prj_nme '/' grp_fle]);
    grp_fle(cellfun(@isempty,grp_fle)) = {''};

%% Load Data
% Redcap Data %%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.red_fle = red_fle;
fcfg.sbj_nme = grp_fle(:,1);
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(fcfg);

% Make Subject Groups %%%%%%%%%%%%%%%%%%%%%
grp_fle = [ grp_fle [ {'SideOnset'} ; sbj_sze.sbj_sde_ons ] ];

fcfg = [];
fcfg.grp_fle = grp_fle;
grp_fle = ejk_collate_groups(fcfg);

% Load fMRI Data %%%%%%%%%%%%%%%%%%%%%Fib
lng_reo_fmri_load

lng_reo_tract_load

%% Data Exploration
% Raw Data %%%%%%%%%%%%%%%%%%%%%
% fMRI (# of Voxels)
lng_reo_fmri_check

% Fibers FA
lng_reo_fiber_check

% Fibers MD


% Laterality %%%%%%%%%%%%%%%%%%%%%
% fMRI Laterality (# of Voxels)
lng_reo_fmri_lat_check

% Fibers Laterality FA

% Fibers Laterality MD











