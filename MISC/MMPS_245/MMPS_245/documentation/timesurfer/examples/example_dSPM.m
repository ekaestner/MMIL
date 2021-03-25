function prefix = run_dSPM()

% run_dSPM

%% this is an example script for running ts_dSPM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set variables
subjdir = '/space/pogo1/6/dhagler/subjects';
subjname = 'larap';
session = '060325LP_test';
datamatfile = 'matfiles/avg_data_post.mat';
badchanfile = '060325LP_badchans.txt';
transfile = '060325LP_mri2head.trans';
conditions = [15];
bem_flag = 1;
write_mgh_flag = 1;
sparsesmooth = 0;
postsmooth = 0;
%lh_sourcecov_wfile = 'singecc-lh.w';
%rh_sourcecov_wfile = 'singecc-rh.w';
lh_sourcecov_wfile = [];
rh_sourcecov_wfile = [];
sourcecov_thresh = 1;
useEEG_flag = 0;
usemag_flag = 0;
usegrad_flag = 1;
ncov_type = 1;
SNR=10;
conduct_scalefact = 1;
cen_sph = [0,2,56];
radii = [65.8,71.9,119.6];
nlayers = 3;
conductivities = [0.3 0.012 0.3];
overwrite_forward_flag = 0;
overwrite_inverse_flag = 0;

prefix = 'dSPM';
chanstr = [];
if usegrad_flag, chanstr = [chanstr 'grad']; end;
if usemag_flag, chanstr = [chanstr 'mag']; end;
if useEEG_flag, chanstr = [chanstr 'EEG']; end;
if bem_flag
  bemstr = 'bem';
else
  bemstr = 'sph';
end;

setenv('SUBJECTS_DIR',subjdir);

prefix = sprintf('%s_%s_%s',...
  prefix,chanstr,bemstr);

% load data
load(datamatfile);

% run dSPM
ts_dSPM(...
  avg_data,subjname,...
  'prefix',prefix,...
  'badchanfile',badchanfile,...
  'conditions',conditions,...
  'ncov_type',ncov_type,...
  'lh_sourcecov_wfile',lh_sourcecov_wfile,...
  'rh_sourcecov_wfile',rh_sourcecov_wfile,...
  'sourcecov_thresh',sourcecov_thresh,...
  'bem_flag',bem_flag,...
  'transfile',transfile,...
  'write_mgh_flag',write_mgh_flag,...
  'sparsesmooth',sparsesmooth,...
  'postsmooth',postsmooth,...
  'usegrad_flag',usegrad_flag,...
  'usemag_flag',usemag_flag,...
  'useEEG_flag',useEEG_flag,...
  'SNR',SNR,...
  'cen_sph',cen_sph,...
  'nlayers',nlayers,...
  'conductivities',conductivities,...
  'radii',radii,...
  'conduct_scalefact',conduct_scalefact,...
  'overwrite_forward_flag',overwrite_forward_flag,...
  'overwrite_inverse_flag',overwrite_inverse_flag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
