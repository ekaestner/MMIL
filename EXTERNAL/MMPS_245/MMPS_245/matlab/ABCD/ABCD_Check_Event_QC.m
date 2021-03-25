function ABCD_Check_Event_QC(ProjID,varargin)
%function ABCD_Check_Event_QC(ProjID,varargin)
%
% Purpose: create report on status of QC completion by event
%
% Required input:
%   ProjID: project ID  string
%
% Optional input:
%   'infix : filename of input info file
%     {default = 'merged_pcqcinfo'}
%   'outfix':  file suffix for output file
%     {default = 'event_qc_info'}
%  'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  03/09/16 by Don Hagler
% Last Mod: 03/09/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'infix','merged_pcqcinfo',[],...
  'outfix','event_qc_info',[],...
  'forceflag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms.indir = sprintf('%s/MetaData/%s',getenv('HOME'),ProjID);
parms.instem = ProjID;
parms.outdir = sprintf('%s/MetaData/%s',getenv('HOME'),ProjID);
parms.outstem = ProjID;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create pc summary
args = mmil_parms2args(parms);
abcd_check_event_qc(args{:});

