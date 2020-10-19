function MMIL_autoQC_Exams(ProjID,varargin)
% function MMIL_autoQC_Exams(ProjID,[options])
%
% Usage:
%  MMIL_autoQC_Exams(ProjID,'key1', value1,...);
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters that specify study specific information
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%     If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that should contain the following fields:
%      orig, raw, proc, proc_dti, proc_bold
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%  'batchrootdir': top level directory containing output batch job directories
%     {default = /home/$USER/batchdirs}
%  'batchname': name of output batchdir
%     {default = 'MMIL_autoQC_Exams'}
%
% Optional Parameters:
%  'ContainerTypes': cell array of container types to check
%     {default = {'proc' 'proc_dti' 'proc_bold'}}
%  'newflag': [0|1] create jobs for new data only
%    i.e., skip containers that already have final autoQC output file
%    ignored if forceflag=1
%    {default = 1}
%  'qccont_flag': [0|1] put output in separate QC containers
%    if RootDirs.qc not specified, will be set to 0
%    {default = 1}
%  'outdir': output directory
%    relative to ContainerPath
%    if empty, will write output to ContainerPath
%    {default = []}
%  'proc outdir': processing file output directory
%    relative to ContainerPath
%    if empty, will write processed output to ContainerPath
%    {default = []}
%  'infix': file name infix (e.g. [], 'corr')
%    {default = []}
%  'cleanupflag': [0|1] remove temporary files
%    {default = 1}
%  'verbose': [0|1] display status messages and warnings
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Optional File Conversion Parameters:
%  'convert_flag': [0|1] convert files to other format
%    {default = 1}
%  'out_type': output file type
%    supported types: 'nii','mgh','mgz'
%    'nii' format can be used by FSL and AFNI
%    {default = 'nii'}
%  'out_orient': output slice orientation
%    if empty or omitted, keep original orientation
%      e.g. 'LPS', 'RAS', etc.
%    for use with FSL, 'LAS' may be preferred
%    {default = []}
%
% Created:  02/20/16 by Don Hagler
% Last Mod: 05/02/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'batchrootdir',[],[],...
  'batchname','MMIL_autoQC_Exams',[],...
...
  'ContainerTypes',{'proc' 'proc_dti' 'proc_bold'},{'proc' 'proc_dti' 'proc_bold'},...
  'newflag',true,[false true],...
  'qccont_flag',true,[false true],...
  'outdir',[],[],...
  'proc_outdir',[],[],...
  'infix',[],[],...
  'cleanupflag',true,[false true],...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
... % file conversion
  'convert_flag',true,[false true],...
  'out_type','nii',{'nii','nii.gz'},...
  'out_orient',[],[],...
... % hidden
  'DTI_meas_list',{'FA','MD','b0'},[],...
  'DTI_err_list',{'rms_err'},{'err','rms_err','rms'},...
  'BOLD_meas_list',{'mean','std','snr'},[],...
  'required_containers',{'raw'},[],...
  'outdir_group_flag',false,[false true],...
  ...
  'info_tags',{'VisitIDs','SubjIDs','StudyInfo','RootDirs','ignore_VisitInfo_flag',...
               'user','numvec_tags'},[],...
  'qc_tags',{'outdir','proc_outdir','infix',...
             'cleanupflag','verbose','forceflag',...
             'convert_flag','out_type','out_orient','DTI_meas_list',...
             'DTI_err_list','BOLD_meas_list','outdir_group_flag'},[],...
});

if parms.forceflag, parms.newflag = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = mmil_parms2args(parms,parms.info_tags);
[ProjInfo,StudyInfo,RootDirs] = MMIL_Quick_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;
  
% check for qc field in RootDirs
if parms.qccont_flag && isempty(mmil_getfield(RootDirs,'qc'))
  parms.qccont_flag = 0;
  if parms.verbose
    fprintf('%s: WARNING: no qc RootDir specified, so setting qccont_flag = 0\n',...
      mfilename);
  end;
end;  

if parms.qccont_flag
  parms.outdir_group_flag = 1;
end;

if ~isempty(parms.batchrootdir)
  RootDirs.batch = parms.batchrootdir;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output batch directory
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s\n',batchdir);
  fprintf('cmd = %s',cmd);
  unix(cmd);
end;
mmil_mkdir(batchdir);

fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create job scripts

j = 1;
for s=1:length(StudyInfo)
  % select study info for this subject visit
  VisitID = StudyInfo(s).VisitID;
  % loop over Container types
  for k=1:length(parms.ContainerTypes)
	  ContainerType = parms.ContainerTypes{k};
	  [ContainerPath,ContainerDir,ContainerRootDir] = ...
  	  MMIL_Get_Container(RootDirs,VisitID,ContainerType);
    if ~isempty(ContainerDir) && exist(ContainerPath,'dir')
      % check if data files exist
      flist = dir(sprintf('%s/*.mgz',ContainerPath));
      if isempty(flist)
        fprintf('%s: no mgz files found in %s, skpping %s...\n',...
          mfilename,ContainerPath,VisitID);
        continue;
      end;
      if parms.qccont_flag
        parms.outdir = sprintf('%s/%s',RootDirs.qc,...
                                       regexprep(ContainerDir,'PROC','QC'));
      else
        parms.outdir = ContainerPath;
      end;
      % check if autoQC has already been done
      if parms.newflag
        flist = dir(sprintf('%s/*_auto_qcinfo.mat',parms.outdir));
        if ~isempty(flist)
%          fprintf('%s: newflag=1, skipping %s...\n',mfilename,VisitID);
          continue;
        end;
      end;
      % create script
      fprintf('%s: creating job for %s...\n',mfilename,VisitID);
      jstem = regexprep(VisitID,'\^','_');
      jstem = jstem(1:min(20,length(jstem)));
      jobID = sprintf('job_%03d_%s_%s',j,jstem,ContainerType); j = j+1;
      jobfname = [batchdir '/' jobID '.m'];
      mmil_write_script(jobfname,'MMIL_autoQC_Exam',{ContainerPath},...
                        parms.qc_tags,parms);
      % add to list
      fid = fopen(scriptlistfname,'a');
      fprintf(fid,'%s\n',jobID);
      fclose(fid);
    end;
  end;
end

% check available disk space
for k=1:length(parms.ContainerTypes)
	ContainerType = parms.ContainerTypes{k};
  MMIL_Check_Usage(RootDirs.(ContainerType));
end;

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmajjobs4 %s\n',parms.batchname);

