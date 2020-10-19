function ABCD_RadReview(ProjID,varargin)
%function ABCD_RadReview(ProjID,varargin)
%
% Purpose: identify T1- and T2-weighted images 
%            with acceptable quality and that are protocol compliant
%          and send the corresponding tgz files for review by radiologist
%
% Required input:
%   ProjID: project ID  string
%
% Optional Input: ('key', value pairs)
%   'rootdir': input directory containing project directories
%     {default = '/space/syn05/1/data/MMILDB'}
%   'fname_info': spreadsheet containing PC and QC info for each series
%     if not supplied, will use merged_pcqcinfo.csv in MetaData
%     {default = []}
%   'fname_enrollment': spreadsheet containing pGUIDs for enrolled subjects only
%     if not supplied, will be set to {home}/MetaData/{ProjID}_enrollment.csv'
%     {default = []}
%   'fname_review': name of file created to record date when files were sent
%     for radiologist review
%     if not supplied, will be set to {home}/MetaData/{ProjID}_radreview.csv'
%     {default = []}
%   'logfile': name of file created to log messages
%     {default = []}
%   'local_user': user name on local server
%     {default = 'ABCDRadReview'}
%   'local_ip': IP address of local server
%     {default = '169.228.56.166'}
%   'local_indir': input directory on local server
%     {default = '/data/home/acquisition_sites'}
%   'local_subdir': input subdirectory directory on local server
%       relative to site directory within local_indir
%     {default = 'fiona/outbox/'}
%   'remote_user': user name on remote server
%     {default = 'abcd'}
%   'remote_ip': IP address of remote server
%     {default = 'researchradiology.com'}
%   'remote_outdir': output directory on remote server
%     {default = '/uploads'}
%   'enflag': [0|1] only send data for enrolled subjects
%     if enflag = 1, fname_enrollment must exist
%     {default = 0}
%   'qcflag': [0|1] only send data with good QC
%     {default = 0}
%   'allflag': [0|1] send all image types (T1 and T2) or send nothing
%     if 0, send all available
%     {default = 0}
%
% Created:  11/26/16 by Don Hagler
% Last Mod: 05/19/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'rootdir','/space/syn05/1/data/MMILDB',[],...
  'indir',[],[],...
  'fname_info',[],[],...
  'fname_enrollment',[],[],...
  'fname_review',[],[],...
  'fname_deplck',[],[],...
  'logfile',[],[],...
  'local_user','ABCDRadReview',[],...
  'local_ip','169.228.56.166',[],...
  'local_indir','/data/home/acquisition_sites',[],...
  'local_subdir','fiona/outbox/',[],...
  'remote_user','abcd',[],...
  'remote_ip','researchradiology.com',[],...
  'remote_outdir','/uploads',[],...
  'enflag',false,[false true],...
  'qcflag',false,[false true],...
  'allflag',false,[false true],...
...
  'series_types',{'T1','T2'},[],...
  'date_pat','yyyy-mm-dd',[],...
...
  'tags',{'indir','fname_info','fname_enrollment','fname_review','logfile',...
          'local_user','local_ip','local_indir','local_subdir',...
          'remote_user','remote_ip','remote_outdir','qcflag','allflag',...
          'series_types','date_pat'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID);
if ~isempty(RootDirs.incoming)
  parms.indir = RootDirs.incoming;
else
  parms.indir = sprintf('%s/%s/incoming',parms.rootdir,ProjID);
end;
if ~exist(parms.indir,'dir')
  error('input directory %s not found',parms.indir);
end;
if isempty(parms.fname_info)
  parms.fname_info = sprintf('%s/MetaData/%s/%s_merged_pcqcinfo.csv',...
    RootDirs.home,ProjID,ProjID);
end;
if ~exist(parms.fname_info,'file')
  error('info file %s not found',parms.fname_info);
end;
if parms.enflag
  if isempty(parms.fname_enrollment)
    parms.fname_enrollment = sprintf('%s/MetaData/%s/%s_enrolledsubjects.csv',...
      RootDirs.home,ProjID,ProjID);
  end;
  if ~exist(parms.fname_enrollment,'file')
    error('enrollment file %s not found',parms.fname_enrollment);
  end;
else
  parms.fname_enrollment = [];
end;
if isempty(parms.fname_review)
  parms.fname_review = sprintf('%s/MetaData/%s/%s_radreview.csv',...
    RootDirs.home,ProjID,ProjID);
end;
if isempty(parms.logfile)
  logdir = sprintf('%s/logs',RootDirs.home);
  mmil_mkdir(logdir);
  parms.logfile = sprintf('%s/%s_radreview.log',logdir,ProjID);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run radreview
args = mmil_parms2args(parms,parms.tags);
abcd_radreview(args{:});

