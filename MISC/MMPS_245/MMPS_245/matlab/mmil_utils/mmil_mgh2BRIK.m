function mmil_mgh2BRIK(fname_in,fname_out,frames,out_orient,forceflag)
%function mmil_mgh2BRIK(fname_in,[fname_out],[frames],[out_orient],[forceflag])
%
% Required Input:
%   fname_in: full or relative path name of input mgh file
%   fname_out: full or relative path name of output BRIK file
%
% Optional Input:
%   frames: vector of frame numbers to extract from mgh file
%     if empty or omitted, will convert all frames
%     {default = []}
%   out_orient' - output orientation (e.g. 'LPS', 'RAS', etc.)
%     if empty or omitted, will keep original orientation
%     {default = []}
%   forceflag: [0|1] toggle overwrite of existing BRIK file
%     {default = 0}
%
% Created:  12/01/06 Don Hagler
% Last Mod: 11/01/12 by Don Hagler
%

if ~mmil_check_nargs(nargin, 1), return; end;

status = 0;

if ~exist('fname_out','var'), fname_out = []; end;
if ~exist('frames','var'), frames = []; end;
if ~exist('out_orient','var'), out_orient = []; end;
if ~exist('forceflag','var'), forceflag=0; end;

% check input file exists
if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;

% check that fname_in has .mgh extension
[fpath,fstem,fext] = fileparts(fname_in);
if ~ismember(fext,{'.mgh','.mgz'})
  error('input file %s lacks proper file extension (mgh or mgz)',fname_in);
end;
if isempty(fpath), fpath = pwd; end;

if isempty(fname_out)
  [fname_out_path,fname_out_fstem] = fileparts(fname_in);
else
  [fname_out_path,fname_out_fstem] = fileparts(fname_out);
end;
[tmp,tmp_fstem] = fileparts(tempname);
tmp_fstem = fullfile(fname_out_path,tmp_fstem);

if ~isempty(frames)
  [vol,M,mrparms]=fs_load_mgh(fname_in,[],frames);
  fname_in = sprintf('%s.mgh',tmp_fstem);
  fs_save_mgh(vol,fname_in,M,mrparms);
end;

% if fname_out has +orig.BRIK extension, get file stem
if ~isempty(fname_out)
  k = findstr(fname_out,'+orig.BRIK');
  if isempty(k)
    fstem_out = fname_out;
  else
    fstem_out = fname_out(1:k(end)-1);
  end;
  fname_out = sprintf('%s+orig.BRIK',fstem_out);
else
  fname_out = sprintf('%s/%s+orig.BRIK',fpath,fstem);
  fstem_out = [fpath '/' fstem];
end;

if exist(fname_out,'file')
  if forceflag
    % remove existing BRIK because AFNI will not overwrite
    fprintf('%s: deleting existing BRIK file %s\n',mfilename,fname_out);
    cmd = sprintf('rm %s+orig.*',fstem_out);
    [status,msg] = unix(cmd);
    if status
      error('failed to remove existing output:\n%s\n%s',...
        cmd,msg);
    end;
  else
    return;
  end;
end;

fprintf('%s: converting %s to BRIK...\n',mfilename,fname_in);
% convert to nifti format
niftifile = sprintf('%s.nii',tmp_fstem);
fs_mri_convert(fname_in,niftifile,...
  'out_orient',out_orient,'forceflag',forceflag);

% convert from nifti to BRIK
cmd = sprintf('3dcopy %s %s',niftifile,fstem_out);
[status,msg] = unix(cmd);
if status || ~exist(fname_out,'file')
  error('failed to create BRIK file:\n%s\n%s',...
    cmd,msg);
end;

% remove temporary nifti file
cmd = sprintf('rm %s.*',tmp_fstem);
[status,msg] = unix(cmd);
if status
  error('failed to remove temporary files %s.*:\n%s\n%s',...
    tmp_fstem,cmd,msg);
end;

