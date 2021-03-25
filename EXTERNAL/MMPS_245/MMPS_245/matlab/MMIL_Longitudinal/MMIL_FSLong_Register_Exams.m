function MMIL_FSLong_Register_Exams(ProjID,varargin)
%function MMIL_FSLong_Register_Exams(ProjID,[options])
%
% Purpose: Registers a FreeSurfer timepoint to the base template for FreeSurfer longitudinal analysis
% Prerequsite: MMIL_FSLong_Base_Exams must be run to create base template prior to this stage.
%
% Required Parameters:
%  'ProjID': project name added to ContainerInfo
%    {default = []}
%
% Optional Parameters:
%  'RootDirs':
%     struct containing the following fields: proc
%     these specify the locations of data
%  'StudyInfo': struct array of study information
%    (e.g. read from csv file with MMIL_Read_StudyInfo)
%    must contain these fields: SubjID, StudyDate, VisitNumber
%    may contain these fields: proc
%    if proc not unspecified, will look for Containers
%      with SubjID and StudyDate (will choose first one if more than one)
%    if empty, use all subjects found in RootDirs.proc
%    {default = []}
%  'batchname': name of output batchdir
%    {default = 'fslong_register'}
%  'FS_version': which version of Freesurfer to use (e.g. 305, 450, 510)
%    {default = 530}
%  'forceflag' - [0|1] whether to create script even if recon is complete
%    {default = 0}
%
% Output:
%   cmd: csh command to run for FreeSurfer recon-all -base
%   err: string containing any error messages
%
% Created:  04/29/16 by Sean Hatton
% Last Mod: 05/18/16 by Sean Hatton
%
% based on code by Don Hagler
cmd = []; err = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{ ...
  'ProjID',[],[],...
  'RootDirs',[],[],...
  'StudyInfo',[],[],...
  'batchname','fslong_register',[],...
  'FS_version',530,[],...
  'forceflag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[StudyInfo,RootDirs] = MMIL_Get_StudyInfo(ProjID);
SubjIDs = {StudyInfo.SubjID};
VisitIDs = {StudyInfo.VisitID};
Subj_List = unique(SubjIDs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output batch directory
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s\n',batchdir);
  fprintf('cmd = %s',cmd);
  [status,result] = unix(cmd);
  if status
    fprintf('%s: WARNING: cmd %s failed:\n%s',mfilename,cmd,result);
  end;
end;
[s,msg] = mkdir(batchdir);
if ~s
  error('failed to create batch job directory %s',batchdir);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create job scripts

scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
fid = fopen(scriptlistfname,'w'); fclose(fid);

for i=1:length(VisitIDs)
  VisitID = VisitIDs{i};
  SubjID = SubjIDs{i};
  
  % Set and check containers exists
  FSContainerPath = MMIL_Get_Container(RootDirs,VisitID,'fsurf');
  if isempty(FSContainerPath) | ~exist(FSContainerPath,'dir')
    fprintf('%s: WARNING: fsurf container %s not found\n',...
      mfilename,FSContainerPath);
  end;
  [tmp,VisitID_Container,ext] = fileparts(FSContainerPath);
  VisitID_Container = [VisitID_Container,ext];
  FSLongBaseContainerPath = MMIL_FSLong_Get_Container(RootDirs,SubjID,'fslong');
  if isempty(FSLongBaseContainerPath) | ~exist(FSLongBaseContainerPath,'dir')
    fprintf('%s: WARNING: fsurf base template container %s not found\n',...
      mfilename,FSLongBaseContainerPath);
  end;
  [tmp,FSLongBase_Container] = fileparts(FSLongBaseContainerPath);
  
  %create script
  fprintf('Creating freesurfer longitudinal script for %s with base template %s...\n',VisitID,FSLongBase_Container);
  tmp_cmd = sprintf('#!/bin/csh\n');
  tmp_cmd = sprintf('%ssource $PUBSH/bin/SetUpFreeSurfer.csh %d\n',...
    tmp_cmd,parms.FS_version);
  tmp_cmd = sprintf('%ssetenv SUBJECTS_DIR %s\n\n\n',tmp_cmd,...
      RootDirs.fslong);
  tmp_cmd = sprintf('%scd %s\n',tmp_cmd,RootDirs.fslong);
  tmp_cmd = sprintf('%srecon-all -long %s %s -all',tmp_cmd,...
      VisitID_Container,FSLongBase_Container);
  tmp_cmd = sprintf('%s\n\n%s\n\n',tmp_cmd);
  cmd = tmp_cmd;

  jobID = sprintf('job_%03d_%s',i,VisitID); i = i+1;
  jobfname = sprintf('%s/%s.csh',batchdir,jobID);
  fid = fopen(jobfname,'w');
  fprintf(fid,'%s',cmd);
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qcshjobs2 %s\n',parms.batchname);

