function MMIL_Long_Register_DTI_Exam(dirA,dirB,varargin)
%function MMIL_Long_Register_DTI_Exam(dirA,dirB,[options])
%
% Purpose: Longitudinal DTI registration
%
% Required Input:
%   dirA: full path of baseline directory
%   dirB: full path of followup directory
%
% Optional Input:
%   'resDTI_flag': [0|1] whether to register DTI resolution images
%     otherwise use T1 resolution images
%     {default: 0}
%   'imgtypes': cell array of image types for which to do registrations
%     {default: {'T1','FA'}}
%   'forceflag': run registrations even if output files exist
%     {default = 0}
%
% Created:  08/14/09 by Don Hagler
% Rcnt Mod: 03/13/12 by Vijay Venkatraman
% Last Mod: 09/15/12 by Don Hagler
%

if (~mmil_check_nargs(nargin,2)) return; end;
parms = mmil_args2parms(varargin, { ...
  'resDTI_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'parmdir',[],[],...
  'binfile','mriRegister',[],...
  'imgtypes',{'T1','FA'},{'T1','FA','ADC'},...
  'ext','.mgz',{'.mgh','.mgz'},...
});

if isempty(parms.parmdir)
  parms.parmdir = [getenv('MMPS_PARMS') '/LONGDTI'];
end;

if parms.resDTI_flag
  dir_infix = 'resDTI';
else
  dir_infix = 'resT1';
end;

if ~iscell(parms.imgtypes)
  parms.imgtypes = {parms.imgtypes};
end;

[dirA_path,dirA_stem,dirA_ext] = fileparts(dirA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
for i=1:length(parms.imgtypes)
  imgtype = parms.imgtypes{i};
  fprintf('%s: registering %s images...\n',mfilename,imgtype);
  fname_volA = [dirA '/images_' dir_infix '/' imgtype parms.ext];
  fname_volB = [dirB '/images_' dir_infix '/' imgtype parms.ext];
  fname_maskA = [dirA '/images_' dir_infix '/' imgtype '_mask' parms.ext];
  fname_maskB = [dirB '/images_' dir_infix '/' imgtype '_mask' parms.ext];
  paramfile = [parms.parmdir '/' imgtype '_NonLinReg_Parms.txt'];
  outdir = [dirB '/nonlinreg_' imgtype '_' dir_infix '_'...
            dirA_stem,dirA_ext];

  % if registering FA or ADC, use affine registration from T1
  if ismember('T1',parms.imgtypes)
    switch imgtype
      case {'FA','ADC'}
        T1_outdir = [dirB '/nonlinreg_T1_' dir_infix '_'...
            dirA_stem,dirA_ext];
        fname_affine_mat = [T1_outdir '/affineRegMatrixT2S.txt'];
        if ~exist(fname_affine_mat,'file')
          fprintf('%s: ERROR: affine registration matrix file %s not found\n',...
            mfilename,fname_affine_mat);
          return;
        end;
      otherwise
        fname_affine_mat = [];
    end;
  end;

  mmil_mriRegister(fname_volA,fname_volB,...
    'outdir',outdir,...
    'fname_maskA',fname_maskA,'fname_maskB',fname_maskB,...
    'fname_affine_mat',fname_affine_mat,...
    'paramfile',paramfile,'binfile',parms.binfile,...
    'forceflag',parms.forceflag);
end;

