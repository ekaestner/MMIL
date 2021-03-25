function errcode = MMIL_Process_Structurals(ContainerPath,varargin)
%function errcode = MMIL_Process_Structurals(ContainerPath,[options])
%
% Required Input:
%   ContainerPath: full path of MRIPROC Container
%
% Optional Input:
%  'refT1flag': which type of T1 series ('MPR' or 'hiFA') to use as reference
%     i.e. which one gets registered to atlas and to which other scans are registered
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default=2}
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
%  'atlasflag': [0|1] whether to resample output in atlas space
%              (rigid-body registration only)
%     {default = 1}
%  'nativeflag': [0|1] if atlasflag=0, whether to keep native resolution
%     otherwise, resample to 1mm, 256^3, LIA
%     {default = 0}
%  'rawQCWflag': [0|1] whether to use manual raw QC info
%     {default = 0}
%  'ProjID': project name (required if rawQCflag=1)
%     {default = []}
%  'minmax': vector of minimum and maximum number of scans required
%     {default = [1 Inf]}
%  'scantypes': cell array of scan types to register
%     {default = {'MPR','XetaT2','FLASHhi','FLASHlo'}}
%  'atlasdir': full path of atlas directory
%     {default =  [getenv('MMPS_DIR') '/atlases']}}
%  'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%     {default =  'T1_Atlas/T1_atlas'}
%  'oldreg_flag': [0|1] use old registration mriRegister
%     otherwise use mmil_reg in mmil_average_volumes
%     {default = 0}
%  'prereg_flag': [0|1] whether to use reg for rigid registration
%     before running dct registration to atlas
%     {default = 1}  
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Created:  06/25/08 by Don Hagler
% Prev Mod: 04/17/16 by Don Hagler
% Last Mod: 05/09/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'refT1flag',2,[0:3],...
  'gradunwarp_flag',1,[0:2],...  
  'wmbc_flag',false,[false true],...
  'nu_flag',false,[false true],...
  'tal_flag',false,[false true],...
  'atlasflag',true,[false true],...
  'nativeflag',false,[false true],...
  'rawQCflag',false,[false true],...
  'ProjID',[],[],...
  'minmax',[1 Inf],[1 Inf],...
  'scantypes',{'MPR','XetaT2','FLASHhi','FLASHlo'},[],...
  'atlasdir',[],[],...
  'atlasname','T1_Atlas/T1_atlas',[],...
  'oldreg_flag',false,[false true],...
  'prereg_flag',true,[false true],...
  'forceflag',false,[false true],...
...
  'getgw_tags',{'gradunwarp_flag','forceflag'},[],...
  'apply_tags',{'gradunwarp_flag','wmbc_flag','nu_flag','tal_flag','forceflag'},[],...
  'reg_tags',{'atlasflag','nativeflag','rawQCflag','ProjID','minmax',...
              'scantypes','refT1flag','gradunwarp_flag',...
              'wmbc_flag','nu_flag','atlasdir',...
              'atlasname','oldreg_flag','prereg_flag','forceflag'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get gradwarp signature from image if ambiguous (GE only)
args = mmil_parms2args(parms,parms.getgw_tags);
errcode = MMIL_Get_Gradwarp_Signature(ContainerPath,args{:});
if errcode ~=0, return; end;

% B1 intensity correction and gradient unwarping (GE and Siemens only)
% Bias-field Intensity correction
args = mmil_parms2args(parms,parms.apply_tags);
errcode = MMIL_Apply_MRI_Corrections(ContainerPath,args{:});
if errcode ~= 0, return; end;

% register structural images, averaging within scan type, register to atlas
args = mmil_parms2args(parms,parms.reg_tags);
errcode = MMIL_Register_and_Resample_MRI_Volumes(ContainerPath,args{:});
if errcode ~= 0, return; end;

MMIL_Calculate_T1Ratio(ContainerPath,parms.forceflag);

