function epi_unwarpB0(fname_in,fname_dx,varargin)
%function epi_unwarpB0(fname_in,fname_dx,[options])
%
% Required Parameters:
%   fname_in:    input file name
%   fname_dx:    input file name for displacement field
%
% Optional Parameters:
%  'fname_out': output file name (B0 unwarped)
%     If not supplied, will append 'B0uw' to file stem of fname_in
%    {default = []}
%  'revflag': [0|1] whether to reverse direction of displacements
%    {default = 0}
%  'swapxy_flag': [0|1] whether to swap x and y dimensions before unwarping B0
%    {default = 0}
%  'verbose': [0|1] display status messages
%    {default = 0}
%  'forceflag': [0|1] whether to overwrite existing fname_out
%    {default = 0}
%
% Created:  12/19/08 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'fname_out',[],[],...
  'revflag',false,[false true],...
  'swapxy_flag',false,[false true],...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
...
  'cleanupflag',true,[false true],...
  'outext','.mgz',{'.mgz','.mgh'},...
});

% create fname_B0dx from fname_for and fname_ref stems
if isempty(parms.fname_out)
  [tpath,tstem,text] = fileparts(fname_in);
  parms.fname_out = [tpath '/' tstem '_B0uw' parms.outext];
end;
if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

if ~exist(fname_in,'file')
  error('input file %s not found',fname_in);
end;
if ~exist(fname_dx,'file')
  error('displacement field file %s not found',fname_dx);
end;

% check that input and voldx share dimensions and orientations
[dimsmatch,orientmatch] = fs_vol_match(fname_in,fname_dx);
if ~dimsmatch
  error('mismatch in input volume dimensions for %s and %s',...
    fname_in,fname_dx);
end;
if ~orientmatch
  error('mismatch in volume orientation for %s and %s',...
    fname_in,fname_dx);
end;

% get relative name for output file
[tmp_path,tmp_fstem,tmp_ext] = fileparts(parms.fname_out);
if isempty(tmp_path), tmp_path = pwd; end;

% if swapxy_flag, swap x and y
if parms.swapxy_flag
  orig_orient = fs_read_orient(fname_in);
  pref_orient = orig_orient([2,1,3]);

  tmp_fname_out = 'out.mgz';
  % create tmp directory
  tmpdir = [tmp_path '/tmp_B0uw'];
  [success,msg] = mkdir(tmpdir);
  if ~success
    error('failed to create tmp dir %s:\n%s',tmpdir,msg);
  end;
  if parms.verbose
    fprintf('%s: reorienting input volumes to %s...\n',mfilename,pref_orient);
  end;
  tmp_fname_in = [tmpdir '/in.mgz'];
  tmp_fname_dx = [tmpdir '/dx.mgz'];
  epi_reorient_vol(fname_in,tmp_fname_in,pref_orient);
  epi_reorient_vol(fname_dx,tmp_fname_dx,pref_orient);
  outdir = tmpdir;
else
  tmp_fname_out = [tmp_fstem tmp_ext];
  tmp_fname_in = fname_in;
  tmp_fname_dx = fname_dx;
  outdir = tmp_path;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate command string
cmd = 'applyUnwarpB0';
if parms.revflag
  cmd = sprintf('%s -r %s -ro %s',cmd,tmp_fname_in,tmp_fname_out);  
else
  cmd = sprintf('%s -f %s -fo %s',cmd,tmp_fname_in,tmp_fname_out);  
end;
cmd = sprintf('%s -od %s -d %s',cmd,outdir,tmp_fname_dx);

% run command
if parms.verbose
  fprintf('%s: cmd = %s\n',mfilename,cmd);
end;
[status,result]=unix(cmd);
if status
  error('cmd %s failed:\n%s\n',cmd,result);
elseif parms.verbose
  disp(result);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parms.swapxy_flag
  if parms.verbose
    fprintf('%s: reorienting output volume to %s...\n',mfilename,orig_orient);
  end;
  tmp_fname_out = [tmpdir '/' tmp_fname_out];
  epi_reorient_vol(tmp_fname_out,parms.fname_out,orig_orient);
end;

if parms.cleanupflag && parms.swapxy_flag
  % remove tmp dir
  cmd = sprintf('rm -r %s',tmpdir);
  [status,result]=unix(cmd);
  if status
    fprintf('%s: WARNING: failed to remove tmp dir:\n%s\n',mfilename,result);
  end;
end;

