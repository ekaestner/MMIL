function MMIL_FSReconQC(ProjID,fname_tcl)
%function MMIL_FSReconQC(ProjID,[fname_tcl])
%
% Required Input:
%   ProjID: MMPS Project ID
%
% Optional Input:
%   fname_tcl: name of tcl script
%     if not full path, will look in user's home directory, then MMPS/tcl
%     {default = 'qc_recon.tcl'}
%
% Created:  01/22/2016 by Don Hagler
% Last Mod: 06/15/2016 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('fname_tcl','var') || isempty(fname_tcl)
  fname_tcl = 'qc_recon.tcl';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt_stem = 'nu.mgz lh.white -aux brainmask.mgz -aux-surface rh.white -segmentation aparc+aseg.mgz';

needed_files = {
  '/mri/aparc+aseg.mgz'...
  '/mri/nu.mgz' ...
  '/mri/brainmask.mgz'...
  '/surf/rh.pial'...
  '/surf/lh.pial'...
  '/surf/rh.white'...
  '/surf/lh.white'...
};

headers = '"VisitID","QC","Motion","PialOver","WMUnder","WMOver","Artifact","Notes"';
default_values = '1,0,0,0,0,0,';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for color lut file
fname_lut = [getenv('FREESURFER_HOME') '/FreeSurferColorLUT.txt'];
if ~exist(fname_lut,'file')
  error('file %s not found',fname_lut);
end;

% check tcl file
if mmil_isrelative(fname_tcl)
  tcl_dir_list = {[getenv('HOME') '/tcl'],[getenv('MMPS_DIR') '/tcl']};
  for i=1:length(tcl_dir_list)
    fname_tmp = sprintf('%s/%s',tcl_dir_list{i},fname_tcl);
    if exist(fname_tmp,'file')
      fname_tcl = fname_tmp;
      break;
    end;
  end;
end;
if ~exist(fname_tcl,'file')
  error('file %s not found',fname_tcl);
end;

cmd_options = sprintf('%s %s -tcl %s',opt_stem,fname_lut,fname_tcl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[StudyInfo,RootDirs,ProjInfo] = MMIL_Get_StudyInfo(ProjID);

if isfield(ProjInfo,'FS_version')
  FS_version = ProjInfo.FS_version;
else
  FS_version = str2num(getenv('FREESURFER_VER'));
end;

setenv('SUBJECTS_DIR',RootDirs.fsurf);

fname_qc = sprintf('%s/ProjInfo/%s/%s_FSReconQC.csv',...
  RootDirs.home,ProjID,ProjID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_qc,'file')
  cmd = sprintf('echo %s >> %s',headers,fname_qc);
  unix(cmd);
end;
qc_info = mmil_csv2struct(fname_qc);
qc_VisitIDs = {qc_info.VisitID};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s=1:length(StudyInfo)
  VisitID = StudyInfo(s).VisitID;

  % check whether QC already done for this VisitID
  if ismember(VisitID,qc_VisitIDs)
    fprintf('%s: NOTE: skipping %s (QC already done)\n',...
      mfilename,VisitID);
    continue;
  end

  % check that STRUCT_VisitID is VisitID
  STRUCT_VisitID = StudyInfo(s).STRUCT_VisitID;
  if ~strcmp(VisitID,STRUCT_VisitID)
    fprintf('%s: NOTE: skipping %s (STRUCT_VisitID is %s)\n',...
      mfilename,VisitID,STRUCT_VisitID);
    continue;
  end

  % check whether FS recon exists
  FSContainerDir = StudyInfo(s).fsurf;
  if isempty(FSContainerDir)
    fprintf('%s: WARNING: skipping %s (recon missing)\n',...
      mfilename,VisitID);
    continue;
  end;

  % check that recon is complete
  FSContainerPath = [RootDirs.fsurf '/' FSContainerDir];
  [status,message] = MMIL_Get_FSReconStatus(FSContainerPath,FS_version);
  if status ~= 2 & status ~= 5
    fprintf('%s: WARNING: skipping %s (recon incomplete)\n',...
      mfilename,FSContainerDir);
    continue;
  end;

  % check that required files exist
  all_exist = 1;
  for nf = 1:length(needed_files);
    fname_test = sprintf('%s/%s',FSContainerPath,needed_files{nf});
    if ~exist(fname_test,'file')
      all_exist = 0;
      break;
    end
  end
  if ~all_exist
    fprintf('%s: WARNING: skipping %s (missing required files)\n',...
      mfilename,VisitID);
    continue;
  end

  tic

  % enter VisitID and default values into spreadsheet
  cmd = sprintf('echo "%s",%s >> %s',VisitID,default_values,fname_qc);
  unix(cmd);

  % run tkmedit
  cmd = sprintf('tkmedit %s %s',FSContainerDir,cmd_options);
  mmil_unix(cmd);

  toc

  % prompt user to update the QC spreadsheet
  fprintf('did you UPDATE entries in %s?\n',fname_qc);
  fprintf('press a key to continue\n');
  pause

end

