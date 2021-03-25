function MMIL_Classify_Exams(ProjID,varargin)
% function MMIL_Classify_Exams(ProjID,[options])
%
% Purpose: re-run MRI classification on already unpacked data
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%       (e.g. read from VisitInfo.csv file with MMIL_Read_StudyInfo)
%    Use this option to limit processing to specific subjects
%    Should have a field 'VisitID' that indicates the name of the
%       input data directory (i.e. in orig data directory)
%    If VisitID field is missing, will use SubjID field instead
%    May also specify subject-specific parameter values
%     but note that these will override ProjInfo and command line input
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include these fields: raw, proc
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'DCM_classify_set': string specifying default classify file
%     in MMPS/parms/classify e.g. MMIL_Series_Classify_{set}.csv
%     'All', 'GE', 'Philips', 'Siemens', or 'Strict'
%     Note: 'Strict' is same as 'All', excluding those that
%       rely on SeriesDescription
%     May also use 'none', which will not use classification rules
%       instead rely on function mmil_classify_by_code which
%       does nothing unless overridden with custom-modified version of it
%     {default = 'All'}
%  'DCM_classify_file': full path name of csv file with rules
%     for classifying series of dicoms
%     if empty, and file called <ProjID>_Series_Classify.csv is
%       found in ProjInfo/<ProjID>, will use that
%     {default = []}
%  'batchname': name of output batchdir
%     {default = 'MMIL_Classify_Exams'}
%  'clustflag': [0|1] whether to create cluster jobs or run locally
%    {default = 1}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Created:  09/23/11 by Don Hagler
% Last Mod: 07/31/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'DCM_classify_set','All',{'All','GE','Philips','Siemens','Strict','none'},...
  'DCM_classify_file',[],[],...
  'batchname','MMIL_Classify_Exams',[],...
  'clustflag',true,[false true],...
  'forceflag',false,[false true],...
...
  'DCM_RejectSeries',[],[],... % to be overridden by StudyInfo
...
  'required_rootdirs',{'orig','raw'},[],...
};
parms = mmil_args2parms(varargin,parms_filter);

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if ~isempty(ProjInfo)
  % For arg names present in both varargin and ProjInfo
  % the varargin values will appear in merged_args
  ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
  merged_args = mmil_merge_args(varargin,ProjInfo_args);
  % check that parameters fit allowed range, use defaults if not supplied
  parms = mmil_args2parms(merged_args,parms_filter);
end;

% set series classify file
if isempty(parms.DCM_classify_file) && ~isempty(ProjID)
  fname = sprintf('%s/ProjInfo/%s/%s_Series_Classify.csv',...
    RootDirs.home,ProjID,ProjID);
  if exist(fname,'file'), parms.DCM_classify_file = fname; end
end;
if isempty(parms.DCM_classify_file) &&...
   ~strcmp(lower(parms.DCM_classify_set),'none')
  parms.DCM_classify_file = sprintf('%s/classify/MMIL_Series_Classify_%s.csv',...
    getenv('MMPS_PARMS'),parms.DCM_classify_set);
end;
if ~isempty(parms.DCM_classify_file) && ~exist(parms.DCM_classify_file,'file')
  error('classify file %s not found',parms.DCM_classify_file);
end;

if parms.clustflag
  % create output batch directory
  if isempty(parms.batchname)
    parms.batchname = [mfilename];
  end;
  if ~isempty(ProjID)
    parms.batchname = [ProjID '_' parms.batchname];
  end;
  batchdir = [RootDirs.batch '/' parms.batchname];
  scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
  if exist(batchdir,'file')
    cmd = sprintf('rm -rf %s/*\n',batchdir);
    fprintf('cmd = %s',cmd);
    unix(cmd);
  else
    mmil_mkdir(batchdir);
  end
  fid = fopen(scriptlistfname,'w'); fclose(fid);
  j = 1;
else
  % check available disk space
  MMIL_Check_Usage(RootDirs.raw);
end;

RootDir = RootDirs.raw;
regpat = '^MRIRAW_(?<VisitID>[\w-+^]+)_(?<StudyDate>\d{8})\..+';
dirlist = dir(sprintf('%s/*',RootDir));

for i=1:length(dirlist)
  ContainerDir = char(dirlist(i).name);
  n = regexp(ContainerDir,regpat,'names','once');
  if isempty(n) || ~dirlist(i).isdir
    continue;
  end;
  VisitID = n.VisitID;
  StudyDate = n.StudyDate;

  if ~isempty(StudyInfo)
    % check that VisitID is in StudyInfo  
    ind = find(strcmp(VisitID,{StudyInfo.VisitID}));
    if isempty(ind)
%      fprintf('%s: WARNING: VisitID %s not found in StudyInfo... skipping\n',...
%        mfilename,VisitID);
      continue;
    end;
    if length(ind)>1
      error('%s: VisitID %s is found %d times in StudyInfo - must be unique\n',...
        mfilename,VisitID,length(ind));
      return;
    end;
    RejectSeries = StudyInfo(ind).DCM_RejectSeries;
  else
    RejectSeries = [];
  end;

  ContainerPath = [RootDirs.raw '/' ContainerDir];
  if parms.clustflag
    jstem = regexprep(VisitID,'\^','_');
    jstem = jstem(1:min(20,length(jstem)));
    jobID = sprintf('job_%03d_%s',j,jstem); j = j+1;
    jobfname = sprintf('%s/%s.m',batchdir,jobID);
    args = {ContainerPath,RejectSeries,parms.DCM_classify_file,parms.forceflag};
    mmil_write_script(jobfname,'MMIL_Classify_Dicoms',args,[],[]);
    fid = fopen(scriptlistfname,'a');
    fprintf(fid,'%s\n',jobID);
    fclose(fid);
  else
    MMIL_Classify_Dicoms(ContainerPath,RejectSeries,...
      parms.DCM_classify_file,parms.forceflag);
  end;
end

if parms.clustflag
  % check available disk space
  MMIL_Check_Usage(RootDirs.raw);

  fprintf('%%%% Now login to a cluster and run this:\n',mfilename);
  fprintf('    qmatjobs %s\n',parms.batchname);
end;
