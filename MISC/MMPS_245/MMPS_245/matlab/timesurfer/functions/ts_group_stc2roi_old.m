function results = ts_group_stc2roi(subjnames,stcfiles,labelfiles,hemi,surfname,subjdir);
%function results = ts_group_stc2roi(subjnames,stcfiles,labelfiles,hemi,[surfname],[subjdir]);
%
% Usage:
%  results = ts_group_stc2roi(subjnames,stcfiles,labelfiles,hemi);
%
% Required input:
%  subjnames: cell array containing list of FreeSurfer subject names
%
%  stcfiles: cell array containing list of stc files
%    (full or relative path)
%    (source time course files output from ts_dSPM)
%
%  labelfiles: struct array containing field called 'fnames'
%    fnames should be a cell array containing a list of file names
%    corresponding to label files or w files defining ROIs
%    (regions of interest)
%    e.g. {'ppc-lh.label', 'sfc-lh.label'} or {'ppc-lh.w', 'sfc-lh.w'}
%    if only one ROI is used, this can be either a cell array or just a string
%
%  hemi: cortical hemisphere
%    should be either 'lh' for left hemisphere or 'rh' for right hemi
%
% Optional input:
%  surfname: name of freesurfer surface to open
%    e.g. white, smoothwm
%    default: 'white'
%
%  subjdir: subjects directory (override SUBJECTS_DIR environment variable)
%    {default = $SUBJECTS_DIR}
%
% Output:
%   results: struct array (one element for each subject) containing fields:
%     time: vector of time points (in msec)
%     wforms: matrix of average source waveforms size = [tpoints,nrois]
%
%  NOTE: subjnames, stcfiles, and labelfiles should all have same number of elements
%
%  created:       09/21/06   by Don Hagler
%  last modified: 07/31/08   by Don Hagler
%

%% todo: use ts_stc2roi (instead of ts_stc2roi_old)
%%       allow stcfiles and labelfiles to be cell matrices
%%       (nsubs x condition for stcfiles and nsubs x roi for labelfiles)

if nargin < 4, help(mfilename); return; end;

if ~exist('surfname','var')
  surfname = 'white';
end;

if ~exist('subjdir','var')
  subjdir = [];
end;

results = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input parameters

% check that input arrays have matching sizes
nsubs = length(subjnames);
nstcfiles = length(stcfiles);
nlabelsets = length(labelfiles);
if nsubs~=nstcfiles
  fprintf('%s: error: # subjects (%d) does not match # stc files (%d)\n',...
    mfilename,nsubs,nstcfiles);
  return;
end;
if nsubs~=nlabelsets
  fprintf('%s: error: # subjects (%d) does not match # label sets (%d)\n',...
    mfilename,nsubs,nlabelsets);
  return;
end;

% check that stc files exist
for i=1:nstcfiles
  if ~exist(stcfiles{i},'file')
    fprintf('%s: error: stcfile %s not found\n',...
      mfilename,stcfiles{i});
    return;
  end;
end;

% check that each subject has same number of label files
% check that label files exist
nrois = length(labelfiles(1).fnames);
for i = 1:nlabelsets
  nlabelfiles = length(labelfiles(i).fnames);
  if nrois~=nlabelfiles
    fprintf('%s: error: for subj %s, # rois (%d) does not match # label files (%d)\n',...
      mfilename,subjnames{i},nrois,nlabelfiles);
    return;
  end;
  for j=1:nrois
    labelfile = labelfiles(i).fnames{j};
    if ~exist(labelfile,'file')
      fprintf('%s: error: labelfile %s not found\n',...
        mfilename,labelfile);
      return;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate average waveforms for each ROI

for i=1:nsubs
  stcfile = stcfiles{i};
  roifilelist = labelfiles(i).fnames;
  subj = subjnames{i};
  [results(i).time,results(i).wforms]=...
    ts_stc2roi_old(stcfile,roifilelist,subj,hemi,surfname,subjdir);
end;

