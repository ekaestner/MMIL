function [pol_stem,ecc_stem,errcode] = BOLD_MMIL_Check_RetFit_Data(ContainerPath,varargin)
%function [pol_stem,ecc_stem,errcode] = BOLD_MMIL_Check_RetFit_Data(ContainerPath,[options])
%
% Purpose: verify that ContainerPath has polar and eccentricity retinotopy
%   data for use with MMIL_RetFit_BOLD_Exam
%
% Required Input:
%   ContainerPath: full path of processed MRI Container (with fMRI retinotopy)
%
% Optional Parameters for Input Data:
%  'multisess_flag': [0|1] indicates source of input analyses
%    0: use results from one or more scans within session
%    1: use results from multiple sessions
%    if 1, pol_dir, ecc_dir, etc. will be ignored
%    {default = 0}
%  'pol_dir': polar angle BOLD analysis subdirectory of ContainerPath
%    if empty, will attempt to find suitable value
%    {default = []}
%  'ecc_dir': eccentricity BOLD analysis subdirectory of ContainerPath
%    if empty, will attempt to find suitable value
%    {default = []}
%  'pol_stem': polar angle file stem
%    if empty, will attempt to find suitable value
%    {default = []}
%  'ecc_stem' eccentricity file stem
%    if empty, will attempt to find suitable value
%    {default = []}
%  'pol_snums': polar angle BOLD scan numbers
%    if empty, will be set to odd scans (for or rev depending on SessInfo.revflag)
%    {default = []}
%  'ecc_snums': eccentricity BOLD scan numbers
%    if empty, will be set to even scans (for or rev depending on SessInfo.revflag)
%    {default = []}
%  'infix': string inside BOLD file names (e.g. 'corr_resBOLD')
%    {default = 'corr_resBOLD'}
%  'fstats_infix': string inside Fourier analysis file names
%    {default = 'fstats_pval'}
%
% Created:  02/22/11 by Don Hagler
% Last Mod: 03/21/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{ ...
  'multisess_flag',false,[false true],...
  'pol_dir',[],[],...
  'ecc_dir',[],[],...
  'pol_stem',[],[],...
  'ecc_stem',[],[],...
  'pol_snums',[],[],...
  'ecc_snums',[],[],...
  'infix','corr_resBOLD',[],...
  'fstats_infix','fstats_pval',[],...
...
  'fnamestem','BOLD',[],...
});

pol_stem = []; ecc_stem = []; errcode = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datatypes = {'pol','ecc'};
for d=1:length(datatypes)
  datatype = datatypes{d};
  [stem,errcode] = set_fstem(ContainerPath,parms,datatype);
  if errcode, return; end;
  flist = dir(sprintf('%s*',stem));
  if isempty(flist)
    errcode = 1;
    return;
  end;
  switch datatype
    case 'pol'
      pol_stem = stem;
    case 'ecc'
      ecc_stem = stem;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stem,errcode] = set_fstem(ContainerPath,parms,datatype)
  stem = [];
  errcode = 0;
  switch datatype
    case 'pol'
      full_datatype = 'polar';
    case 'ecc'
      full_datatype = 'eccen';
  end;

  if parms.multisess_flag
    adir = sprintf('%s_multisess_%s_analysis',...
      parms.fnamestem,full_datatype);
    stem = sprintf('%s_multisess_%s_avg',...
      parms.fnamestem,full_datatype);
    snums = [];
  else
    stem = parms.([datatype '_stem']);
    adir = parms.([datatype '_dir']);
    snums = parms.([datatype '_snums']);
    if isempty(snums)
      errcode = 1;
      return;
    end;
  end;

  [stem,errcode] = BOLD_MMIL_Set_Fourier_Stem(ContainerPath,...
    'infix',parms.infix,'fstats_infix',parms.fstats_infix,...
    'fnamestem',parms.fnamestem,...
    'snums',snums,'stem',stem,'adir',adir);
return;


