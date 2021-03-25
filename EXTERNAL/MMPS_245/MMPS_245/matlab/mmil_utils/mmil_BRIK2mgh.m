function mmil_BRIK2mgh(fname_in,fname_out,forceflag)
%function mmil_BRIK2mgh(fname_in,[fname_out],[forceflag])
%
% Required Input:
%   fname_in: full or relative path name of input BRIK file
%   fname_out: full or relative path name of output mgh/mgz file
%
% Optional Input:
%   forceflag: [0|1] toggle overwrite of existing mgh file
%     {default = 0}
%
% Created:  06/15/09 Don Hagler
% Rcnt Mod: 02/25/13 Mark Fiecas
% Last Mod: 02/28/13 Don Hagler
%

if ~mmil_check_nargs(nargin, 1), return; end;

if ~exist('fname_out','var'), fname_out = []; end;
if ~exist('forceflag','var') | isempty(forceflag), forceflag=0; end;

% check input BRIK file exists
if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;
% get input file stem
[in_path,in_fstem,in_ext] = fileparts(fname_in);
if isempty(in_path), in_path = pwd; end;
if ~strcmp(in_ext,'.BRIK')
  error('input file %s lacks proper file extension (BRIK)',fname_in);
end;
k = regexp(in_fstem,'\+(?!.*\+)');
if isempty(k)
  error('input file name %s has wrong pattern for BRIK',fname_in);
end;
in_fstem_name = in_fstem(1:(k-1));
in_fstem_space = in_fstem((k+1):end);

% check input HEAD file exists
fname_head = [in_path '/' in_fstem_name '+' in_fstem_space '.HEAD'];
if ~exist(fname_head,'file')
  error('file %s not found',fname_head);
end;

if isempty(fname_out)
  % create output file name from input file stem
  out_path = in_path;
  out_fstem = in_fstem_name;
  out_ext = '.mgh';
else
  % check that fname_out has .mgh extension
  [out_path,out_fstem,out_ext] = fileparts(fname_out);
  if isempty(out_path), out_path = pwd; end;
  if ~ismember(out_ext,{'.mgh','.mgz'})
    error('output file %s lacks proper file extension (mgh or mgz)',fname_out);
  end;
end;
fname_out = [out_path '/' out_fstem out_ext];
if exist(fname_out,'file') & ~forceflag, return; end;

[tmp,tmp_fstem] = fileparts(tempname);
fstem_tmp = [out_path '/' tmp_fstem];

% use 3dTcat to de-bucketize (e.g. 3dDeconvolve output)
fname_tmp = [fstem_tmp '+orig.BRIK'];
%fprintf('%s: converting %s to 3D+time format in %s...\n',...
%  mfilename,fname_in,fname_tmp);
cmd = sprintf('3dTcat -prefix %s %s',fstem_tmp,fname_in);
[status,msg] = unix(cmd);
if status || ~exist(fname_tmp,'file')
  if regexp(msg,'with only 1 sub-bricks!')
    fname_tmp = fname_in;
  else
    error('failed to concat BRIK:\n%s\n%s',...
      cmd,msg);
  end;
end;

% convert from BRIK to nifti
fname_tmp2 = [fstem_tmp '.nii'];
%fprintf('%s: converting %s to %s...\n',mfilename,fname_tmp,fname_tmp2);
cmd = sprintf('3dAFNItoNIFTI -float -prefix %s %s',fstem_tmp,fname_tmp);
[status,msg] = unix(cmd);
if status || ~exist(fname_tmp2,'file')
  error('failed to convert BRIK to nii:\n%s\n%s',...
    cmd,msg);
end;

% convert to mgh format
%fprintf('%s: converting %s to %s...\n',mfilename,fname_tmp,fname_out);
fs_mri_convert(fname_tmp2,fname_out,'forceflag',forceflag);

% remove temporary file
cmd = sprintf('rm %s*',fstem_tmp);
[status,msg] = unix(cmd);
if status
  error('failed to remove temporary files with stem %s:\n%s\n%s',...
    fstem_tmp,cmd,msg);
end;

