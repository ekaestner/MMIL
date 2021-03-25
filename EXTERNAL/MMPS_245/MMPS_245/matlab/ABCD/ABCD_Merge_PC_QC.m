function ABCD_Merge_PC_QC(ProjID,varargin)
%function ABCD_Merge_PC_QC(ProjID,varargin)
%
% Required input:
%   ProjID: project ID  string
%
% Optional parameters:
%   'pcinfo_instem': input file stem ofor protocol compliance
%     file expected to be found here:
%       /home/{user}/MetaData/{ProjID}/{ProjID}_{pcinfo_instem}.csv
%     {default = 'pcinfo'}
%   'qcinfo_instem': input file stem for auto and raw QC summary
%     file expected to be found here:
%       /home/{user}/MetaData/{ProjID}/{ProjID}_{qcinfo_instem}.csv
%     {default = 'qcinfo'}
%   'merged_outstem': output file stem for merger between PC and QC
%     {default = 'merged_pcqcinfo'}
%   'merge_field': column header used to merge PC and QC summaries
%     {default = 'VisitID'}
%   'forceflag': [0|1] overwrite existing output
%     {default = 1}     
%
% Created:  10/22/16 by Don Hagler
% Last Mod: 03/20/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'pcinfo_instem','pcinfo',[],...
  'qcinfo_instem','qcinfo',[],...
  'merged_outstem','merged_pcqcinfo',[],...
  'merge_field','VisitID',[],...                      
  'forceflag',true,[false true],...
  ...
  'merge_flag',3,[0:3],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms.outdir = sprintf('%s/MetaData/%s',getenv('HOME'),ProjID);
parms.pcinfo_instem = sprintf('%s_%s',ProjID,parms.pcinfo_instem);
parms.qcinfo_instem = sprintf('%s_%s',ProjID,parms.qcinfo_instem);
parms.merged_outstem = sprintf('%s_%s',ProjID,parms.merged_outstem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% merge pc and qc summaries
fname_merged = sprintf('%s/%s.csv',parms.outdir,parms.merged_outstem);
if ~exist(fname_merged,'file') || parms.forceflag
  fname_pcinfo = sprintf('%s/%s.csv',parms.outdir,parms.pcinfo_instem);
  fname_qcinfo = sprintf('%s/%s.csv',parms.outdir,parms.qcinfo_instem);
  if exist(fname_pcinfo,'file') && exist(fname_qcinfo,'file')
    mmil_merge_csv(fname_pcinfo,fname_qcinfo,fname_merged,...
      parms.merge_field,parms.merge_flag,parms.forceflag);
  else
    if ~exist(fname_pcinfo,'file')
      fprintf('%s: WARNING: PC summary file %s not found\n',...
        mfilename,fname_pcinfo);
    end;      
    if ~exist(fname_qcinfo,'file')
      fprintf('%s: WARNING: QC summary file %s not found\n',...
        mfilename,fname_qcinfo);
    end;      
  end;
end;

