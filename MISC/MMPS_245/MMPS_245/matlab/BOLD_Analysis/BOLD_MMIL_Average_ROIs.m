function BOLD_MMIL_Average_ROIs(ProjID,varargin)
% function BOLD_MMIL_Average_ROIs(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string, used to to load ProjInfo and StudyInfo
%     from user's home directory
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%    e.g. read from csv file with MMIL_Read_StudyInfo
%    If empty, will use ProjID to get StudyInfo
%    May control which subjects are included in average with
%      field 'GroupAvg_ROI'
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    If empty, will use ProjID to get RootDirs
%    {default = []}
%  'qcflag': [0|1] whether to exclude subjects with StudyInfo.QC=0
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0 }
%
% Optional Parameters that specify source of ROIs (from retinotopy fitting):
%  'RF_outdir': output directory
%    full path or relative to ContainerPath
%    {default = 'retfit'}
%  'RF_roi_name': file stem of ROI file (e.g. lh.roi_name.label)
%    {default = 'v123'}
%
% Optional Parameters for averaging ROIs
%  'outdir': output directory
%    may be absolute path or will be relative to /home/{user}/MetaData/{ProjID}
%    {default = 'Average_ROI'}
%  'outstem': output file stem (do not include hemi or extension)
%    {default = 'v123-groupavg'}
%  'outtype': output file type ('mgh' or 'label')
%    {default = 'label'}
%  'smooth_subj': smoothing steps before sampling to ico
%    {default = 0}
%  'smooth_ico': smoothing steps after sampling to ico
%    {default = 0}
%  'thresh': threshold applied to average of binary masks
%    {default =0}
%  'verbose': [0|1] display detailed output
%   {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  03/17/11 by Don Hagler
% Last Mod: 11/15/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin, { ...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'hemilist',{'lh','rh'},[],...
  'qcflag',false,[false true],...
  'forceflag',false,[false true],...
  'StudyInfo',[],[],...
...
  'required_rootdirs',{'proc','fsurf'},[],...
...
  'RF_outdir','retfit',[],...
  'RF_roi_name','v123',[],...
...
  'outdir','Average_ROIs',[],...
  'outstem','v123-groupavg',[],...
  'outtype','label',{'mgh','label'},...
  'hemi',[],{'lh','rh'},...
  'smooth_subj',0,[0 Inf],...
  'smooth_ico',0,[0 Inf],...
  'thresh',0,[0,1],...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'forceflag',false,[false true],...
  'verbose',false,[false true],...
...
  'icolevel',7,[0 7],...
  'fpat','(?<hemi>[l,r]h)\.(?<name>.+)\.label',[],...
  'avg_tag','GroupAvg_ROI',[],...
});

tags_avg = {'outdir','outstem','outtype','hemi','smooth_subj',...
  'smooth_ico','thresh','subjdir','forceflag','verbose','icolevel','fpat'};
tags_RF = {'RF_outdir','RF_roi_name'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = MMIL_Args(parms,'MMIL_Check_ProjID');
[ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
if isempty(StudyInfo), error('empty StudyInfo'); end;
parms.subjdir = RootDirs.fsurf;

if isfield(StudyInfo,parms.avg_tag)
  tmp = {StudyInfo.(parms.avg_tag)};
  ind_avg = find(~cellfun(@isempty,tmp));
  tmp2 = zeros(length(StudyInfo));
  tmp2(ind_avg) = tmp{ind_avg};
  ind_avg = find(tmp2);
  if isempty(ind_avg)
    fprintf('%s: WARNING: no sessions with %s = 1\n',...
      mfilename,parms.avg_tag);
    return;
  end;
  StudyInfo = StudyInfo(ind_avg);
end;

nsubs = length(StudyInfo);

if  mmil_isrelative(parms.outdir)
  parms.outdir = ...
    sprintf('%s/MetaData/%s/%s',...
  RootDirs.home,ProjID,parms.outdir);
end;
mmil_mkdir(parms.outdir);

fprintf('%s: averaging ROIs across subjects...\n',mfilename);
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  switch parms.outtype
    case 'mgh'
      fname_out = [parms.outdir '/' parms.outstem '-' hemi '.' parms.outtype];
    case 'label'
      fname_out = [parms.outdir '/' hemi '.' parms.outstem '.' parms.outtype];
  end;
  if ~exist(fname_out,'file') || parms.forceflag
    fname_log = [parms.outdir '/' hemi '.' parms.outstem '.log'];
    flog = fopen(fname_log,'wt');
    if flog<0, error('failed to open %s for writing',fname_log); end;

    subjnames = {StudyInfo.fsurf};
    labels = cell(1,nsubs);
    for i=1:nsubs
      SubjID = StudyInfo(i).SubjID;
      VisitID = StudyInfo(i).VisitID;
      ContainerPath = sprintf('%s/%s',RootDirs.proc,StudyInfo(i).proc);
      % replace values in parms with (non-empty) values from StudyInfo
      tmp_parms = parms;
      tmp_parms.hemi = parms.hemilist{h};
      for t=1:length(tags_RF)
        tmp_val = mmil_getfield(StudyInfo(i),tags_RF{t},[]);
        if ~isempty(tmp_val), tmp_parms.(tags_RF{t}) = tmp_val; end;
      end;
      fname = sprintf('%s/%s/%s.%s.label',...
        ContainerPath,tmp_parms.RF_outdir,tmp_parms.hemi,tmp_parms.RF_roi_name);
      if ~exist(fname,'file'), error('file %s not found',fname); end;
      labels{i} = fname;
      fprintf(flog,'# SubjID %s\tVisitID %s#\n%s\n',...
        SubjID,VisitID,fname);
    end
    args = mmil_parms2args(tmp_parms,tags_avg);
    fname_out = fs_label_groupavg(labels,subjnames,args{:});
    fprintf(flog,'\n');
    fprintf(flog,'### wrote output file %s\n',fname_out);
  end;
end;

fclose(flog);

