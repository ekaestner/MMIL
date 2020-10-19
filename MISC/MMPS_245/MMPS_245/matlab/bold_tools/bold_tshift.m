function fname_out = bold_tshift(fname,varargin)
%function fname_out = bold_tshift(fname,[options])
%
% Purpose: perform slice timing correction for BOLD scans using AFNI's 3dTShift
%
% Usage:
%  bold_tshift(fname,'key1', value1,...);
%
% Required Input
%  fname: full path file name of input 4D volume (mgh/mgz format)
%
% Optional Input:
%  'fname_out': output file name of motion corrected 4D volume
%    If empty, will append '_ts' to file stem of fname
%    {default = [])
%  'tpattern': slice time pattern
%    Allowed values: {'alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z'}
%    {default = 'alt+z'}
%  'skipTRs': number of TRs at beginning of scan to be ignored in time shifting
%    If skipTRs >= number of frames in fname, nothing will be done
%    and fname will be returned as fname_out
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Output:
%   fname_out: output file name of slice timing corrected 4D volume
%
% Created:  04/26/10 by Don Hagler
% Last Mod: 08/02/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_out',[],[],...
  'tpattern','alt+z',{'alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z'},...
  'skipTRs',0,[0,Inf],...
  'forceflag',false,[false true],...
...
  'minTRs',5,[],... % hard-coded in 3dTshift (Fatal Error: -ignore value is too large)
  'suffix','ts',[],...
  'tmpdir','tmp_ts',[],...
  'cleanupflag',true,[false true],...
  'verbose',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tpath,tstem,text] = fileparts(fname);
if isempty(tpath), tpath = pwd; end;
if isempty(parms.fname_out)
  fname_out = [tpath '/' tstem '_' parms.suffix text];
else
  fname_out = parms.fname_out;
end;
parms.outdir = fileparts(fname_out);

[M,volsz] = mmil_load_mgh_info(fname,parms.forceflag,parms.outdir);
nframes = volsz(4);

if nframes-parms.minTRs <= max(1,parms.skipTRs)
  if parms.verbose
    fprintf('%s: WARNING: skipping slice timing correction because',mfilename);
    if parms.skipTRs > 1
      fprintf(' skipTRs (%d) > number of frames in %s (%d)\n',...
        parms.skipTRs,fname,nframes);
    else
      fprintf(' of single frame in %s\n',fname);
    end;
  end;
  fname_out = fname;
  return;
end;

if exist(fname_out,'file') & ~parms.forceflag
  return;
end;

if ~exist(fname,'file')
  error('file %s not found',fname);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create tmpdir
if mmil_isrelative(parms.tmpdir)
  parms.tmpdir = [parms.outdir '/' parms.tmpdir];
end;
mmil_mkdir(parms.tmpdir);

% convert input image to nii
fname_tmp = [parms.tmpdir '/' tstem '.nii'];
fs_mri_convert(fname,fname_tmp,'forceflag',parms.forceflag);

% adjust for slice timing with 3dTshift
fname_out_tmp = [parms.tmpdir '/' tstem '_' parms.suffix '.nii'];
% if num TRs < skipTRs, don't do tshift
if ~exist(fname_out_tmp,'file') || parms.forceflag
  if parms.verbose
    fprintf('%s: correcting slice timing for %s\n',mfilename,fname);
  end;
  cmd = sprintf('3dTshift -prefix %s',fname_out_tmp);
  cmd = sprintf('%s -tpattern ''%s''',cmd,parms.tpattern);
  if parms.skipTRs > 0
    cmd = sprintf('%s -ignore %d',cmd,parms.skipTRs);
  end;
  cmd = sprintf('%s %s',cmd,fname_tmp);
  if parms.verbose
    fprintf('%s: cmd = %s\n',mfilename,cmd);
  end;
  [status,result] = unix(cmd);
  if status
    error('3dTshift failed:\n%s\n%s',cmd,result);
  elseif parms.verbose
    disp(result);
  end;
  fname_tmp = fname_out_tmp;
end;

% convert nifti to mgh
fs_mri_convert(fname_out_tmp,fname_out);

% delete temporary files
if parms.cleanupflag && exist(parms.tmpdir,'dir')
  cmd = sprintf('rm -r %s\n',parms.tmpdir);
  [status,result] = unix(cmd);
  if status
    warning('failed to remove tmp dir %s:\n%s',parms.tmpdir,result);
  end;
end;

