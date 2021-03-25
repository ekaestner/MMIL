function [time,wforms] = ts_stc2roi(stcfile,roifilelist)
%function [time,wforms] = ts_stc2roi(stcfile,roifilelist)
%
% Usage:
%  [time,src_waveform] = ts_stc2roi(stcfile,roifilelist);
%
% Required input:
%  stcfile: stc file (full or relative path)
%    (source time course files output from ts_dSPM)
%    e.g. test-lh.stc
%    Can be a cell array of stc file names
%
%  roifilelist: cell array containing a list of file names corresponding
%    to label files or w files defining ROIs (regions of interest)
%    e.g. {'ppc-lh.label', 'sfc-lh.label'} or {'ppc-lh.w', 'sfc-lh.w'}
%    if only one ROI is used, this can be either a cell array or just a string
%
% Output:
%   time: vector of time points (in msec)
%   wforms: matrix of average source waveforms size = [nrois,tpoints]
%     If stcfile is cellarray, will be [nstcs,nrois,tpoints]
%
%
%  created:       09/07/06   by Don Hagler
%  last modified: 01/08/10   by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,2)) return; end;

time = [];
wforms = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(roifilelist)
  roifilelist = {roifilelist};
end;
nrois = length(roifilelist);
fprintf('%s: %d ROIs to process...\n',mfilename,nrois);

if iscell(stcfile)
  stcfilelist = stcfile;
else
  stcfilelist = {stcfile};
end;
nstcs = length(stcfilelist);
if nstcs>1
  fprintf('%s: %d stcfiles to process...\n',mfilename,nstcs);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s=1:nstcs
  stcfile = stcfilelist{s};
  if ~exist(stcfile,'file')
    error('file %s not found',stcfile);
  end;

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
    error('num vertices (%d) does not match num dips (%d)',...
      nverts,ndips);
  end;


  endtime = starttime+(tpoints-1)*sampdur;
  time = starttime:sampdur:endtime;

  fprintf('%s: num vertices = %d\n',mfilename,nverts);
  fprintf('%s: num time points = %d\n',mfilename,tpoints);

  fprintf('%s: extracting ROIs...\n',mfilename);
  if isempty(wforms)
    wforms = zeros(nstcs,nrois,tpoints);
  end;
  for i=1:nrois
    roifile = roifilelist{i};
    fprintf('%s: calculating average waveform for ROI %s\n',...
      mfilename,roifile);
    k = findstr(roifile,'.');
    k = k(end);
    if isempty(k)
      error('missing file extension for roifile %s',roifile);
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
          error('unsupported file type for roifile %s',roifile);
      end;
    catch
      error('failed to read roifile %s',roifile);
    end;
    v_stc = intersect(vertices,v); % only include signals from ROI vertices in stc file
    i_v_stc = find(ismember(vertices,v_stc));
    wforms(s,i,:) = mean(sol(i_v_stc,:),1);
  end;
end;

if nstcs==1
  wforms = reshape(wforms,[nrois,tpoints]);
end;

