function ABCD_reg2mni_taskBOLD_Exams(ProjID,varargin)
%function ABCD_reg2mni_taskBOLD_Exams(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject including fields:
%    SubjID, VisitID
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: proc_bold, fsurf (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied, MMIL_ProjInfo.csv is not required
%    {default = []}
%  'user': user name of account with StudyInfo (i.e. ProjInfo, VisitInfo)
%    If not supplied, will use $USER environment variable
%    {default = []}
%  'root_outdir': full path of root output directory
%    if empty, will use RootDirs.proc_bold
%    {default = []}
%  'batchname': name of output batchdir
%    {default = 'ABCD_reg2mni_taskBOLD_Exams'}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Optional Parameters applicable to data selection:
%  'tasknames': cell array of task names
%    {default = {'MID','SST','nBack'}}
%  'infix': input file will be sprintf('BOLD%d_%s.mgz',snum,infix)
%    to specify no infix (i.e. raw data), use 'none'
%    {default = 'corr_resBOLD'}
%  'concat_flag': [0|1|2] analyze concatenated across scans
%    0: analyze each scan individually
%    1: analyze concatenated scans
%    2: analyze individually and concatenated
%    {default = 1}
%
% Created:  11/17/17 by Don Hagler
% Last Mod: 11/17/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'user',[],[],...
  'root_outdir',[],[],...
  'batchname','ABCD_reg2mni_taskBOLD_Exams',[],...
  'forceflag',false,[false true],...
...
  'tasknames',{'MID','SST','nBack'},{'MID','SST','nBack'},...
  'infix','corr_resBOLD',[],...
  'concat_flag',1,[0:2],...
...
  'smooth_fwhm',4,[],...
...
  'info_tags',{'VisitIDs','SubjIDs','StudyInfo','RootDirs','ignore_VisitInfo_flag',...
               'user','numvec_tags'},[],...
  'reg_tags',{'tasknames','infix','concat_flag',...
              'outdir','verbose','forceflag',...
              'fext','orient','brain_thresh','smooth_fwhm',...
              'reg_infix'},[],...
};

parms = mmil_args2parms(varargin,parms_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = mmil_parms2args(parms,parms.info_tags);
[ProjInfo,StudyInfo,RootDirs] = MMIL_Quick_Check_ProjID(ProjID,args{:});

if isempty(StudyInfo), error('empty StudyInfo'); end;
if strcmp(parms.infix,'none'), parms.infix = []; end;
if isempty(parms.root_outdir)
  parms.root_outdir = RootDirs.proc_bold;
end;

% create output batch directory
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s/*\n',batchdir);
  fprintf('cmd = %s',cmd);
  unix(cmd);
end;
mmil_mkdir(batchdir);

fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VisitIDs = {StudyInfo.VisitID};
nvisits = length(VisitIDs);

fprintf('%s: attempting to create jobs for %d visits...\n',...
  mfilename,nvisits);

j = 1;
for i=1:nvisits
  VisitID = VisitIDs{i};
  % find containers for this VisitID
  [ContainerPath,ContainerDir] = MMIL_Get_Container(RootDirs,VisitID,'proc_bold');
  FSContainerPath = MMIL_Get_Container(RootDirs,VisitID,'fsurf');
  if isempty(ContainerPath)
    fprintf('%s: WARNING: missing proc_bold container for %s\n',...
      mfilename,VisitID);
    continue;
  end;
  if isempty(FSContainerPath)
    fprintf('%s: WARNING: missing fsurf container for %s\n',...
      mfilename,VisitID);
    continue;
  end;
  targs = {ContainerPath,FSContainerPath};
  parms.outdir = sprintf('%s/%s/export_reg2mni',parms.root_outdir,ContainerDir);
  
  % create job script
  jstem = regexprep(VisitID,'\^','_');
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem);
  jobfname = [batchdir '/' jobID '.m'];
  mmil_write_script(jobfname,'abcd_reg2mni_taskBOLD',...
    targs,parms.reg_tags,parms);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
  j=j+1;
end;

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs3 %s\n',parms.batchname);

