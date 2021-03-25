function ABCD_PCQC_REDCap(ProjID,varargin)
%function ABCD_PCQC_REDCap(ProjID,varargin)
%
% Required input:
%   ProjID: project ID  string
%
% Optional input:
%   'infix': input file suffix
%      {default = 'merged_pcqcinfo'}
%
% Created:  02/16/17 by Jose Teruel
% Last Mod: 03/09/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'infix','merged_pcqcinfo',[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms.indir = sprintf('%s/MetaData/%s',getenv('HOME'),ProjID);
parms.outdir = sprintf('%s/MetaData/%s',getenv('HOME'),ProjID);
parms.instem = ProjID;
parms.outstem = ProjID;

if ~exist(parms.indir,'dir')
  error('input directory %s not found',parms.indir);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create pcqc summary
args = mmil_parms2args(parms);
abcd_pcqc_redcap(args{:});

