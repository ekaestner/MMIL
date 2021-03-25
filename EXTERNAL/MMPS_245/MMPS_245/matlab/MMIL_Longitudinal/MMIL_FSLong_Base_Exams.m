function MMIL_FSLong_Base_Exams(ProjID,varargin)
%function MMIL_FSLong_Base_Exams(ProjID,[options])
%
% Purpose: Create the base template for FreeSurfer longitudinal analysis
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
%    {default = 'fslong_base'}
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
% Last Mod: 05/23/16 by Sean Hatton
%
% based on code by Don Hagler
cmd = []; err = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{ ...
  'ProjID',[],[],...
  'RootDirs',[],[],...
  'StudyInfo',[],[],...
  'batchname','fslong_base',[],...
  'FS_version',530,[],...
  'forceflag',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[StudyInfo,RootDirs] = MMIL_Get_StudyInfo(ProjID);
SubjIDs = {StudyInfo.SubjID};
VisitIDs = {StudyInfo.VisitID};
VisitNumbers = [StudyInfo.VisitNumber];
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

for i=1:length(Subj_List)
  SubjID = Subj_List{i};
  ivec = find(strcmp(SubjID,SubjIDs));
  vnums = VisitNumbers(ivec); 
  % check freesurfer container exists
  FSContainerPath = [RootDirs.fsurf '/' StudyInfo(i).fsurf];
  if isempty(StudyInfo(i).fsurf) | ~exist(FSContainerPath,'dir')
    if StudyInfo(i).VisitNumber == 1
      error('baseline fsurf container %s not found',FSContainerPath);
    else
      if ~isempty(StudyInfo(i).fsurf)
        fprintf('%s: WARNING: followup fsurf container %s not found\n',...
          mfilename,FSContainerPath);
      end;
    end;
  end;
  
  % set the fsurf visits for subject  
  nvisits = length(ivec);
  FSPaths = cell(1,nvisits);
  for j=1:nvisits
    FSPaths{j} = MMIL_Get_Container(RootDirs,VisitIDs{ivec(j)},'fsurf');
  end;
  % set base template name and VisitNumbers
  SubjID_base = sprintf('FSLONG_%s_base%s',SubjID,sprintf('_v%d',vnums));
 
  %create script
  fprintf('Creating freesurfer base template script for %s...\n',SubjID)
    tmp_cmd = sprintf('#!/bin/csh\n');
  tmp_cmd = sprintf('%ssource $PUBSH/bin/SetUpFreeSurfer.csh %d\n',...
    tmp_cmd,parms.FS_version);
  tmp_cmd = sprintf('%ssetenv SUBJECTS_DIR %s\n\n\n',tmp_cmd,...
      RootDirs.fslong);
  tmp_cmd = sprintf('%scd %s\n',tmp_cmd,RootDirs.fslong);
  tmp_cmd = sprintf('%sln -s %s/*%s* .\n',tmp_cmd,RootDirs.fsurf,SubjID);
  tmp_cmd = sprintf('%srecon-all -base %s %s-all',tmp_cmd,...
      SubjID_base,sprintf('-tp %s ',FSPaths{:}));
  tmp_cmd = sprintf('%s\n\n%s\n\n',tmp_cmd);
  cmd = tmp_cmd;

  jobID = sprintf('job_%03d_%s',i,SubjID); i = i+1;
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

