function [M,volsz,mr_parms] = mmil_load_mgh_info(fname,forceflag,outdir)
%function [M,volsz,mr_parms] = mmil_load_mgh_info(fname,[forceflag],[outdir])
%
% Purpose: load M and volsz from header of fname
%   and save to mat file, unless mat file already
%   exists, in which case, load from mat file
%
%   Name of mat file will be constructed from file stem
%     of fname, with '_info.mat' appended
%
%   Avoids repeated accessing of compressed mgz files
%
% Required Input:
%   fname: full path file name of mgh or mgz file
%
% Optional Input:
%   forceflag: [0|1] overwrite existing info mat file
%     {default = 0}
%   outdir: output directory
%     {default = path of fname}
%
% Ouptput:
%   M: voxel to RAS transformation matrix
%   volsz: vector of matrix dimensions
%   mr_parms: [tr flipangle te ti fov]
%
% Created:  04/11/11 by Don Hagler
% Prev Mod: 09/08/11 by Don Hagler
% Last Mod: 08/06/17 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
M = []; volsz = []; mr_parms = [];

if ~exist(fname,'file'), error('file %s not found',fname); end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;
if ~exist('outdir','var'), outdir = []; end;

[tpath,tstem,text] = fileparts(fname);
if isempty(outdir), outdir = tpath; end;

if ~ismember(text,{'.mgh','.mgz'})
  error('file %s has wrong extension',fname);
end;

fname_mat = [outdir '/' tstem '_info.mat'];
if ~exist(fname_mat,'file') || forceflag
  [vol,M,mr_parms,volsz] = fs_load_mgh(fname,[],[],1);
  try
    mmil_mkdir(outdir);
    save(fname_mat,'M','volsz','mr_parms');
  catch
    fprintf('%s: WARNING: failed to save info to %s:\n%s\n',...
      mfilename,fname_mat,lasterr);
  end;
else
  load(fname_mat);
end;

