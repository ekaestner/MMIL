function errcode = MMIL_Apply_MRI_Corrections(ContainerPath,varargin)
%function errcode = MMIL_Apply_MRI_Corrections(ContainerPath,[options])
%
% Purpose: correct structural images for B1 intensity inhomogeneities,
%   grad warp distortions, and intensity bias field
%
% Required Input:
%   ContainerPath: full path of MRIPROC container
%
% Optional Input:
%  'gradunwarp_flag': [0|1|2] whether to correct for gradient warping
%     if 2, processing is aborted if gradwarp info is missing
%    {default = 1}
%  'wmbc_flag': [0|1] whether to correct for intensity bias-field
%     using white matter segmentation and sparse smoothing
%    {default = 0}
%  'nu_flag' : [0|1] whether to correct for intensity bias-field
%     using nu (N3) correction
%    {default = 0}
%  'tal_flag': [0|1] register to Talairach space
%    for standardization of white matter intensity values
%    during nu correction
%    {default = 0}
%  'scantypes': cell array of scan types to correct
%    {default={'MPR','FLASHhi','FLASHlo','MEDIChi','MEDIClo','XetaT2'}}
%  'ext': file extension (i.e. '.mgh' or '.mgz')
%    {default = '.mgz'}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Note: nu and wmbc correction are independent
%   wmbc includes nu correction as a pre-processing step (even if nu_flag = 0)
%
% Created:  07/16/08 by Don Hagler
% Prev Mod: 07/17/17 by Don Hagler
% Last Mod: 11/15/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'gradunwarp_flag',1,[0:2],...  
  'wmbc_flag',false,[false true],...
  'nu_flag',false,[false true],...
  'tal_flag',false,[false true],...
  'scantypes',{'MPR','FLASHhi','FLASHlo','MEDIChi','MEDIClo','XetaT2'},[],...
  'ext','.mgz',{'.mgh','.mgz'},...
  'forceflag',false,[false true],...
  ...
  'nu_scantypes',{'MPR','FLASHhi','MEDIChi','XetaT2'},[],...
  'wmbc_scantypes',{'MPR','FLASHhi','MEDIChi','XetaT2'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% B1 correction
errcode = MMIL_B1Corr_Structurals(ContainerPath,...
  'scantypes',parms.scantypes,'ext',parms.ext,'forceflag',parms.forceflag);
if errcode
  fprintf('%s: ERROR: B1 correction failed\n',mfilename);
  return;
end;

% grad warp distortion correction
if parms.gradunwarp_flag
  errcode = MMIL_GradUnwarp_Structurals(ContainerPath,...
    'scantypes',parms.scantypes,'ext',parms.ext,'forceflag',parms.forceflag);
  if errcode
    if parms.gradunwarp_flag==2
      fprintf('%s: ERROR: grad warp distortion correction failed\n',mfilename);
      return;
    else
      fprintf('%s: WARNING: grad warp distortion correction failed\n',mfilename);
    end;
  end;
end;

% nu (N3) bias-field intensity correction
if parms.nu_flag
  % do not allow nu correction on all scan types
  scantypes = intersect(parms.scantypes,parms.nu_scantypes);
  errcode = MMIL_nuCorr_Structurals(ContainerPath,...
    'scantypes',scantypes,'ext',parms.ext,...
    'tal_flag',parms.tal_flag,'forceflag',parms.forceflag);
  if errcode
    fprintf('%s: ERROR: bias-field intensity nu correction failed\n',mfilename);
    return;
  end;
end;

% white matter segmentation-based bias-field intensity correction
if parms.wmbc_flag
  scantypes = intersect(parms.scantypes,parms.wmbc_scantypes);
  errcode = MMIL_wmCorr_Structurals(ContainerPath,...
    'scantypes',scantypes,'ext',parms.ext,...
    'forceflag',parms.forceflag);
  if errcode
    fprintf('%s: ERROR: bias-field intensity wm correction failed\n',mfilename);
    return;
  end;
end;

return;

