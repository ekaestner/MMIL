function ABCD_Check_PC_Exams(ProjID,varargin)
%function ABCD_Check_PC_Exams(ProjID,varargin)
%
% Required input:
%   ProjID: project ID  string
%
% Optional parameters:
%   'batchname': name of batch directory in /home/mmilrec/batchdirs
%     {default = 'ABCD_Check_PC_Exams'}
%   'rootdir': root directory containing project directory
%     with incoming, orig, and pc dirs for ProjID
%     if empty, will get RootDirs from ProjInfo
%     {default = []}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}     
%
% Created:  09/29/16 by Don Hagler
% Prev Mod: 09/29/16 by Don Hagler
% Last Mod: 08/02/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'batchname','ABCD_Check_PC_Exams',[],...
  'rootdir',[],[],...
  'forceflag',false,[false true],...
  ...
  'tags',{'batchname','outdir','forceflag'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.rootdir)
  [ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID);
  indir = RootDirs.orig;
  parms.outdir = RootDirs.pc;
else
  indir = sprintf('%s/%s/orig',parms.rootdir,ProjID);
  parms.outdir = sprintf('%s/%s/pc',parms.rootdir,ProjID);
end;
if isempty(indir)
  error('input directory not specified');
elseif ~exist(indir,'dir')
  error('input directory %s not found',indir);
end;

parms.batchname = [ProjID '_' parms.batchname];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create jobs
args = mmil_parms2args(parms,parms.tags);
abcd_pc_all(indir,args{:});

