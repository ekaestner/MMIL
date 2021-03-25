function ABCD_PQA_Summarize(ProjID,indir,outdir,varargin)
%function ABCD_PQA_Summarize(ProjID,[options])
%
% Required Input:
%   ProjID: project ID string
%
%
% Created:  03/27/17 by Feng Xue
% Last Mod: 09/28/17 by Feng Xue
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',indir,[],...
  'outdir',outdir,[],...
  'ProjID',ProjID,[],...
  'fname_siteinfo',sprintf('%s/ProjInfo/%s/pqa_sites.mat',getenv('HOME'),ProjID),[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if parms.outdir(1) ~= '/', parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir); end;
  if parms.indir(1) ~= '/', parms.indir= sprintf('%s/%s',getenv('HOME'),parms.indir); end;
  if ~exist(parms.outdir), mmil_mkdir(parms.outdir); end;
  if ~exist(parms.indir,'dir'), error('input directory %s not found',parms.indir); end;

  if ~exist(parms.fname_siteinfo,'file'), error('info file %s not found',parms.fname_siteinfo); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = mmil_parms2args(parms);
abcd_pqa_summarize(args{:});
