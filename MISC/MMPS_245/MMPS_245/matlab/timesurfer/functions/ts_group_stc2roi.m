function results = ts_group_stc2roi(subjnames,stcfiles,labelfiles,hemi,surfname,subjdir);
%function results = ts_group_stc2roi(subjnames,stcfiles,labelfiles,hemi,[surfname],[subjdir]);
%
% Usage:
%  results = ts_group_stc2roi(subjnames,stcfiles,labelfiles,hemi);
%
% Required input:
%  subjnames: cell array containing list of FreeSurfer subject names
%
%  stcfiles: cell matrix of stc file names (full or relative path)
%            (source time course files output from ts_dSPM)
%    # rows = # subjects, # columns = # conditions
%
%  labelfiles: cell matrix of file names
%      corresponding to label files or w files defining ROIs
%      (regions of interest)
%      e.g. {'ppc-lh.label', 'sfc-lh.label'} or {'ppc-lh.w', 'sfc-lh.w'}
%    # rows = # subjects, # columns = # ROIs
%    If dSPM for all subjects was run on ico,
%      a cell matrix with 1 column of ROIs may be supplied
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
%       or [nconds,nrois,tpoints] if nconds>1
%
%  created:       09/21/06   by Don Hagler
%  last modified: 01/08/10   by Don Hagler
%

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

if ~iscell(stcfiles)
  error('stcfiles must be a cell matrix');
end;
if nsubs ~= size(stcfiles,1)
  error('# subjects (%d) does not match # of rows in stcfiles (%d)',...
    nsubs,size(stcfiles,1));
end;
nconds = size(stcfiles,2);

if ~iscell(labelfiles)
  error('labelfiles must be a cell matrix');
end;
nlabelsubs = size(labelfiles,1);
nROIs = size(labelfiles,2);
if nlabelsubs==1
  fprintf('%s: will use same set of %d labels for all subjects...\n',...
    mfilename,nROIs);
elseif nsubs ~= nlabelsubs
  error('# subjects (%d) does not match # of rows in labelfiles (%d)',...
    nsubs,nlabelsubs);
end;

% check that stc and label files exist
for s=1:nsubs
  for c=1:nconds
    if ~exist(stcfiles{s,c},'file')
      error('stcfile %s not found',...
        stcfiles{s,c});
    end;
  end;
  if s<=nlabelsubs
    for r=1:nROIs
      if ~exist(labelfiles{s,r},'file')
        error('labelfile %s not found',...
          labelfiles{s,r});
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate average waveforms for each ROI

for s=1:nsubs
  tmp_stcfiles = {stcfiles{s,:}};
  if nlabelsubs>1
    tmp_labelfiles = {labelfiles{s,:}};
  else
    tmp_labelfiles = labelfiles;
  end;
  subj = subjnames{s};
  [results(s).time,results(s).wforms]=...
    ts_stc2roi(tmp_stcfiles,tmp_labelfiles);
end;

