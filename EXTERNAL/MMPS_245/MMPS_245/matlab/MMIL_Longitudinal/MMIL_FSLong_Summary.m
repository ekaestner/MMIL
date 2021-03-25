function MMIL_FSLong_Summary(ProjID,varargin)
%function MMIL_FSLong_Summary(ProjID,[options])
%
% Purpose: Provides FreeSurfer longitudinal analysis statistics
% Prerequsite: MMIL_FSLong_Base_Exams must be run to create base template prior to this stage.
% Prerequsite: MMIL_FSLong_Register_Exams must be run to create timepoint-to-base analysis prior to this stage.
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
%  'outdir': output directory for summary csv files
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'FSLong_Summaries'}
%  'outstem': output file stem
%     {default = 'MRI'}
%  'forceflag' - [0|1] whether to create script even if recon is complete
%    {default = 0}
%
% Output:
%   cmd: csh command to run for FreeSurfer recon-all -base
%   err: string containing any error messages
%
% Created:  05/18/16 by Sean Hatton
% Last Mod: 05/25/16 by Sean Hatton
%

cmd = []; err = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{ ...
  'RootDirs',[],[],...
  'StudyInfo',[],[],...
  'batchname','fslong_summary',[],...
  'FS_version',530,[],...
  'outdir',[],[],...
  'forceflag',true,[false true],...
});

% check outdir
if isempty(parms.outdir)
  parms.outdir = [getenv('HOME') '/MetaData'];
  if ~isempty(ProjID)
    parms.outdir = [parms.outdir '/' ProjID];
  end;
  parms.outdir = [parms.outdir '/FSLong_Summaries'];
elseif mmil_isrelative(parms.outdir)
  parms.outdir = [getenv('HOME') '/' parms.outdir];
end;
% create output directory
mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[StudyInfo,RootDirs] = MMIL_Get_StudyInfo(ProjID);
SubjIDs = {StudyInfo.SubjID};
VisitIDs = {StudyInfo.VisitID};
VisitNumbers = [StudyInfo.VisitNumber];
Periods = VisitNumbers - 1;

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
cmd = []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create QDEC Table

fslong_table = sprintf('%s/FSLong_table.txt',parms.outdir);
fid = fopen(fslong_table,'w');
if fid==-1
  error('failed to open longitudinal FreeSurfer table file %s for writing\n',fslong_table);
end;
fclose(fid);
fid = fopen(fslong_table,'a');
fprintf(fid,'fsid fsid-base period\n');

for i=1:length(VisitIDs)
  VisitID = VisitIDs{i};
  SubjID = SubjIDs{i};
  Period = Periods(i);
  % Set FreeSurfer Visit and Base containers
  FSContainerPath = MMIL_Get_Container(RootDirs,VisitID,'fsurf');
  [tmp,VisitID_Container,ext] = fileparts(FSContainerPath);
  VisitID_Container = [VisitID_Container,ext];
  FSLongBaseContainerPath = MMIL_FSLong_Get_Container(RootDirs,SubjID,'fslong');
  [tmp,FSLongBase_Container] = fileparts(FSLongBaseContainerPath);
  fprintf(fid,'%s %s %d\n',VisitID_Container,FSLongBase_Container,Period);
end;
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create script

scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
fid = fopen(scriptlistfname,'w'); fclose(fid);

  tmp_cmd = sprintf('#!/bin/csh\n');
  tmp_cmd = sprintf('%ssource $PUBSH/bin/SetUpFreeSurfer.csh %d\n',...
    tmp_cmd,parms.FS_version);
  tmp_cmd = sprintf('%ssetenv SUBJECTS_DIR %s\n\n',tmp_cmd,...
      RootDirs.fslong);
  tmp_cmd = sprintf('%scd %s\n\n',tmp_cmd,RootDirs.fslong);
  tmp_cmd = sprintf('%slong_stats_slopes --qdec %s --meas volume --do-avg --do-rate --do-pc1 --do-spc --do-stack --time period --sd $SUBJECTS_DIR --stats aseg.stats --stack-avg %s/aseg.avg.csv --stack-rate %s/aseg.rate.csv --stack-pc1 %s/aseg.pc1.csv --stack-spc %s/aseg.spc.csv\n\n',tmp_cmd,fslong_table,parms.outdir,parms.outdir,parms.outdir,parms.outdir);
  tmp_cmd = sprintf('%slong_stats_slopes --qdec %s --meas thickness --do-avg --do-rate --do-pc1 --do-spc --do-stack --time period --sd $SUBJECTS_DIR --stats lh.aparc.stats --stack-avg %s/lh.aparc.avg.csv --stack-rate %s/lh.aparc.rate.csv --stack-pc1 %s/lh.aparc.pc1.csv --stack-spc %s/lh.aparc.spc.csv\n\n',tmp_cmd,fslong_table,parms.outdir,parms.outdir,parms.outdir,parms.outdir);
  tmp_cmd = sprintf('%slong_stats_slopes --qdec %s --meas thickness --do-avg --do-rate --do-pc1 --do-spc --do-stack --time period --sd $SUBJECTS_DIR --stats rh.aparc.stats --stack-avg %s/rh.aparc.avg.csv --stack-rate %s/rh.aparc.rate.csv --stack-pc1 %s/rh.aparc.pc1.csv --stack-spc %s/rh.aparc.spc.csv\n',tmp_cmd,fslong_table,parms.outdir,parms.outdir,parms.outdir,parms.outdir);
  tmp_cmd = sprintf('%s\n\n%s\n\n',tmp_cmd,cmd);
  cmd = tmp_cmd;

  jobfname = sprintf('%s/%s.csh',batchdir,ProjID);
  fid = fopen(jobfname,'w');
  fprintf(fid,'%s',cmd);
  fclose(fid);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',ProjID);
  fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qcshjobs2 %s\n',parms.batchname);

