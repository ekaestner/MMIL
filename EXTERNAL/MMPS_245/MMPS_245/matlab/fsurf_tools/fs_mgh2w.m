function fs_mgh2w(mghfile,wfilestem,hemi,t0,t1,oneframe_flag)
% fs_mgh2w - convert an mgh surface file to one or more w files
% 
% fs_mgh2w(mghfile,wfilestem,hemi,[t0],[t1])
% 
%  Required input:
%    mghfile: full path name of input mgh file
%    wfilestem: file stem (omit hemisphere) for output wfiles
%    hemi:    cortical hemisphere
%
%  Optional parameters:
%    t0: first time sample to extract
%      {default = 1}
%    t1: last time sample to extract
%      {default = last}
%    oneframe_flag: [0|1] single frame only, output name has not tstamp
%      {default = 0}
%
% created:        07/13/06 Don Hagler
% last modified:  03/15/11 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('t0','var') || isempty(t0), t0=0; end;
if ~exist('t1','var') || isempty(t1), t1=0; end;
if ~exist('oneframe_flag','var') || isempty(oneframe_flag), oneframe_flag = 1; end;

if ~exist(mghfile,'file')
  error('file %s not found\n',mghfile);
end;

% read header
[vol,M,mr_parms,volsz]=fs_load_mgh(mghfile,[],[],1);

nverts = volsz(1);
ny = volsz(2);
nz = volsz(3);
nframes = volsz(4);

if ny~=1 | nz~=1
  error('mgh file is not a surface file -- size is [%d,%d,%d,%d]',...
    nverts,ny,nz,nframes);
end;

if t0<1, t0 = 1; end;
if t1<1 | t1>nframes, t1 = nframes; end;
if t0>t1, t0=t1; end;
frames = [t0:t1];

% read data
vol = fs_load_mgh(mghfile,[],frames);
vol(isnan(vol)) = 0;

if isempty(vol)
  error('failed to read mghfile %s',mghfile);
end;

if oneframe_flag
  nframes = 1;
else
  nframes = length(frames);
end;

for t=1:nframes
  vals = vol(:,1,1,t);
  v = find(vals);
  w = vals(v);
  if oneframe_flag
    outfile = sprintf('%s-%s.w',wfilestem,hemi);
  else
    outfile = sprintf('%s-t%d-%s.w',wfilestem,t0+t-1,hemi);
  end;
  fprintf('%s: writing %s...\n',mfilename,outfile);
  fs_write_wfile(outfile,w,v);
end;

fprintf('%s: finished.\n',mfilename);

return
