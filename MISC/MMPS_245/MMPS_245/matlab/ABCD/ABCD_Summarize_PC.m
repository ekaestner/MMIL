function ABCD_Summarize_PC(ProjID,varargin)
%function ABCD_Summarize_PC(ProjID,varargin)
%
% Required input:
%   ProjID: project ID  string
%
% Optional parameters:
%   'rootdir': root directory containing project directory
%     with incoming, unpack, and orig dirs
%     {default = '/space/syn05/1/data/MMILDB'}
%   'incoming_dir': full path of incoming directory
%     if not set, will be {rootdir}/{ProjID}/incoming
%     {default = []}
%   'outstem': output file stem
%     {default = 'pcinfo'}
%   'col_order': column names in the desired order
%      for output csv file arranged by series
%      additional columns will be moved to the end
%     {default = {'pGUID','VisitID','EventName','SessionType',...
%                 'SeriesType','ABCD_Compliant','SeriesDescription'}}
%   'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  09/29/16 by Don Hagler
% Prev Mod: 12/19/16 by Don Hagler
% Last Mod: 06/08/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'rootdir','/space/syn05/1/data/MMILDB',[],...
  'incoming_dir',[],[],...
  'outstem','pcinfo',[],...
  'col_order',{'pGUID','VisitID','EventName','SessionType','SiteName',...
                      'SeriesType','ABCD_Compliant','SeriesDescription'},[],...
  'forceflag',true,[false true],...
  ...
  'tags',{'auxdir','outdir','outstem','col_order','forceflag'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indir = sprintf('%s/%s/pc',parms.rootdir,ProjID);
if ~exist(indir,'dir')
  error('input directory %s not found',indir);
end;

if isempty(parms.incoming_dir)
  parms.incoming_dir = ...
    sprintf('%s/%s/incoming',parms.rootdir,ProjID);
end;

parms.auxdir = parms.incoming_dir;
parms.outdir = sprintf('%s/MetaData/%s',getenv('HOME'),ProjID);
parms.outstem = sprintf('%s_%s',ProjID,parms.outstem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create pc summary
args = mmil_parms2args(parms,parms.tags);
abcd_summarize_pc(indir,args{:});

