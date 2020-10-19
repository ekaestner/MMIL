function ABCD_Incomingtopc_Report(ProjID,varargin)
%function ABCD_Incomingtopc_Report(ProjID,varargin)
%
% Required input:
%   ProjID: project ID  string
%
% Optional input:
%  'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  02/16/16 by Jose Teruel
% Prev Mod: 03/08/16 by Don Hagler
% Last Mod: 03/09/16 by Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'forceflag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID);

parms.outdir = sprintf('%s/MetaData/%s',getenv('HOME'),ProjID);
parms.outstem = ProjID;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create pc summary
args = mmil_parms2args(parms);
abcd_incomingtopc_report(args{:});
