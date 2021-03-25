function [time,wforms] = ts_stc2roi(stcfile,roifilelist,subj,hemi,surfname,subjdir)
%function [time,src_waveform] = ts_stc2roi(stcfile,roifilelist,subj,hemi,[surfname],[subjdir])
%
% Usage:
%  [time,src_waveform] = ts_stc2roi(stcfile,roifilelist,subj,hemi);
%
% Required input:
%  stcfile: stc file (full or relative path)
%    (source time course files output from ts_dSPM)
%    e.g. test-lh.stc
%
%  roifilelist: cell array containing a list of file names corresponding
%    to label files or w files defining ROIs (regions of interest)
%    e.g. {'ppc-lh.label', 'sfc-lh.label'} or {'ppc-lh.w', 'sfc-lh.w'}
%    if only one ROI is used, this can be either a cell array or just a string
%
%  subj: FreeSurfer subject name
%    environment variable SUBJECT_DIR must be set
%    SUBJECT_DIR/subj should contain the freesurfer subject directory
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
%   time: vector of time points (in msec)
%   wforms: matrix of average source waveforms size = [tpoints,nrois]
%
%
%  created:       09/07/06   by Don Hagler
%  last modified: 09/21/06   by Don Hagler
%

if nargin < 4, help(mfilename); return; end;

time = [];
wforms = [];

if ~exist('surfname','var')
  surfname = 'white';
end;

if ~exist('subjdir','var')
  subjdir = [];
end;

% check input
if ~strcmp(hemi,'lh') & ~strcmp(hemi,'rh')
  fprintf('%s: hemi must be ''lh'' or ''rh''...quitting\n',mfilename);
  return;
end;

if ~iscell(roifilelist)
  roifilelist = {roifilelist};
end;
nrois = length(roifilelist);
fprintf('%s: %d ROIs to process...\n',mfilename,nrois);

% extract time course for ROIs
fprintf('%s: reading stcfile %s...\n',mfilename,stcfile);
[starttime,srate,vertices,sol]=ts_read_stc(stcfile);
sampdur = 1000/srate;
vertices = vertices + 1;

fprintf('%s: start time = %0.2f ms\n',mfilename,starttime);
fprintf('%s: sampling rate = %0.2f Hz\n',mfilename,srate);
fprintf('%s: sample duration = %0.2f ms\n',mfilename,sampdur);

nverts = length(vertices);
[ndips,tpoints] = size(sol);

if (nverts ~= ndips)
  fprintf('%s: error: num vertices (%d) does not match num dips (%d)!\n',...
    mfilename,nverts,ndips);
  return;
end;

endtime = starttime+(tpoints-1)*sampdur;
time = starttime:sampdur:endtime;

fprintf('%s: num vertices = %d\n',mfilename,nverts);
fprintf('%s: num time points = %d\n',mfilename,tpoints);

fprintf('%s: extracting ROIs...\n',mfilename);
wforms = zeros(tpoints,nrois);
for i=1:nrois
  roifile = roifilelist{i};
  fprintf('%s: calculating average waveform for ROI %s\n',...
    mfilename,roifile);
  k = findstr(roifile,'.');
  k = k(end);
  if isempty(k)
    fprintf('%s: missing file extension for roifile %s...skipping\n',...
      mfilename,roifile);
    continue;
  end;
  ext = roifile(k+1:end);
  % determine type of roi file    
  try
    switch ext
      case 'label'
        v=fs_read_label(roifile);
      case 'w'
        [w,v]=fs_read_wfile(roifile);
        v = v(find(w)); % make sure v is only non-zero values
      otherwise
        fprintf('%s: unsupported file type for roifile %s...skipping\n',...
          mfilename,roifile);
        continue;
    end;
  catch
    fprintf('%s: error reading ROI file %s...skipping\n',...
      mfilename,roifile);
  end;
  v_stc = intersect(vertices,v); % only include signals from ROI vertices in stc file
  i_v_stc = find(ismember(vertices,v_stc));
  wforms(:,i) = mean(sol(i_v_stc,:),1);
end;

