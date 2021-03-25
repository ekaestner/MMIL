function qcinfo = MMIL_Summarize_QCInfo(ProjID,varargin)
% function qcinfo = MMIL_Summarize_QCInfo(ProjID,[options])
%
% Usage:
%  qcinfo = MMIL_Summarize_QCInfo(ProjID,'key1', value1,...);
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
%  'RootDirs': struct that must contain the following fields:
%      fsurf, fsico
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%
% Optional Parameters:
%  'ContainerTypes': cell array of container types to convert
%     {default = {'proc' 'proc_dti' 'proc_bold'}}
%  'qccont_flag': [0|1] put output in separate QC containers
%    if RootDirs.qc not specified, will be set to 0
%    {default = 1}
%  'reviewerIDs': cell array of reviewer ID strings included in qcinfo files
%     allows for multiple, independent reviews
%     {default = []}
%  'fname_out': output file name
%    if not supplied, will be use outdir, ProjID, and outfix
%    {default = []}
%  'outdir': output directory
%    if not supplied, will be /home/{user}/MetaData/{ProjID}
%    {default = []}
%  'outfix': suffix attached to output filename
%    {default = 'qcinfo'}
%  'type_sort_flag': [0|1] review exams sorted by container type
%     otherwise, review exams sorted by visit
%    {default = 0}
%  'verbose': [0|1] display status messages and warnings
%    {default: 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 1}
%
% Output:
%   qcinfo: struct array containing autoQC and rawQC information
%     for each series in each visit
%
% Created:  04/03/16 by Don Hagler
% Last Mod: 05/10/17 by Don Hagler
%

%% todo: add sMRI_ftype_list, etc. to header

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qcinfo = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ProjID,varargin);

if ~exist(parms.fname_out,'file') || parms.forceflag
  if parms.type_sort_flag
    % loop over container types, then visits
    for k=1:parms.ntypes
      for s=1:parms.nvisits
        tmpinfo = get_qcinfo(s,k,parms);
        qcinfo = concat_qcinfo(qcinfo,tmpinfo);
      end;
    end
  else
    % loop over visits, then container types
    for s=1:parms.nvisits
      for k=1:parms.ntypes
        tmpinfo = get_qcinfo(s,k,parms);
        qcinfo = concat_qcinfo(qcinfo,tmpinfo);
      end;
    end
  end;
  if isempty(qcinfo)
    fprintf('%s: WARNING: qcinfo is empty\n',mfilename);
    return;
  end;
  % save compiled qcinfo to output file
  mmil_struct2csv(qcinfo,parms.fname_out);
elseif nargout
  qcinfo = mmil_csv2struct(parms.fname_out);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms = mmil_args2parms(options,{...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
  ...
    'ContainerTypes',{'proc' 'proc_dti' 'proc_bold'},{'proc' 'proc_dti' 'proc_bold'},...
    'qccont_flag',true,[false true],...
    'reviewerIDs',[],[],...
    'fname_out',[],[],...
    'outdir',[],[],...
    'outfix','qcinfo',[],...
    'type_sort_flag',false,[false true],...
    'verbose',true,[false true],...
    'forceflag',true,[false true],...
  ... % hidden
    'thresholds',[0.2 0.3 0.4],[],...
    'info_tags',{'RootDirs','ignore_VisitInfo_flag','user','numvec_tags'},[],...
  ... % for QU
    'sMRI_ftype_list',{'motion'},[],...
    'dMRI_ftype_list',[],[],...
    'fMRI_ftype_list',[],[],...
    'sMRI_stype_list',{'MPR','FLASHhi','XetaT2'},[],...
    'dMRI_stype_list',{'DTI'},[],...
    'fMRI_stype_list',{'BOLD'},[],...
  });
  parms.ntypes = length(parms.ContainerTypes);
  parms.nthresh = length(parms.thresholds);

  % get StudyInfo and RootDirs
  if isempty(parms.StudyInfo)
    if parms.verbose
      fprintf('%s: getting StudyInfo for %s...\n',mfilename,ProjID);
      tic;
    end;
    args = mmil_parms2args(parms,parms.info_tags);
    [parms.StudyInfo,parms.RootDirs] = MMIL_Quick_StudyInfo(ProjID,args{:});
    if parms.verbose, toc; end;
  elseif parms.verbose
    fprintf('%s: getting list of VisitIDs from user supplied StudyInfo...\n',...
      mfilename);
  end;
  if isempty(parms.RootDirs)
    [ProjInfo,parms.RootDirs] = MMIL_Get_ProjInfo(ProjID);
  end;
    
  % check for qc field in RootDirs
  if parms.qccont_flag && isempty(mmil_getfield(parms.RootDirs,'qc'))
    parms.qccont_flag = 0;
    if parms.verbose
      fprintf('%s: WARNING: no qc RootDir specified, so setting qccont_flag = 0\n',...
        mfilename);
    end;
  end;  

  % get list of VisitIDs
  parms.VisitIDs = {parms.StudyInfo.VisitID};
  parms.nvisits = length(parms.VisitIDs);
  if parms.verbose
    fprintf('%s: %d visits and %d types to summarize\n',...
      mfilename,parms.nvisits,parms.ntypes);
  end;

  % set output file name
  if isempty(parms.fname_out)
    if isempty(parms.outdir)
      parms.outdir = sprintf('%s/MetaData/%s',parms.RootDirs.home,ProjID);
    end;
    parms.fname_out = sprintf('%s/%s_%s.csv',...
      parms.outdir,ProjID,parms.outfix);
  else
    parms.outdir = fileparts(parms.fname_out);
    if isempty(parms.outdir), parms.outdir = pwd; end;
  end;
  mmil_mkdir(parms.outdir);
  
  % check reviewerIDs
  if ~iscell(parms.reviewerIDs)
    parms.reviewerIDs = {parms.reviewerIDs};
  end;
  parms.nrev = length(parms.reviewerIDs);
  for r=1:parms.nrev
    if isnumeric(parms.reviewerIDs{r})
      parms.reviewerIDs{r} = num2str(parms.reviewerIDs{r});
    end;
  end;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stypes = get_stypes(ContainerPath,parms)
  stypes = [];
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
  if errcode || isempty(ContainerInfo), return; end;
  stypes = fieldnames(ContainerInfo.ScanInfo);
  ind = find(~cellfun(@isempty,struct2cell(ContainerInfo.ScanInfo)));
  stypes = stypes(ind);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qcinfo = get_qcinfo(s,k,parms)
  qcinfo = [];
  VisitID = parms.StudyInfo(s).VisitID;
  ContainerType = parms.ContainerTypes{k};
  if parms.verbose
    fprintf('%s: checking for %s data for %s...\n',mfilename,ContainerType,VisitID);
    tic;
  end;
  [ContainerPath,ContainerDir,ContainerRootDir] = ...
    MMIL_Get_Container(parms.RootDirs,VisitID,ContainerType);
  if parms.verbose, toc; end;
  if ~isempty(ContainerDir) && exist(ContainerPath,'dir')
    if parms.verbose
      fprintf('%s: checking %s QC info for %s...\n',mfilename,ContainerType,VisitID);
      tic;
    end;
    k = 0; % scan index
    % get scan types from ContainerInfo
    stypes = get_stypes(ContainerPath,parms);
    for i=1:length(stypes)
      stype = stypes{i};
      % set list of qualitative scores
      ftype_list = {};
      switch stype
        case parms.sMRI_stype_list
          ftype_list = parms.sMRI_ftype_list;
        case parms.dMRI_stype_list
          ftype_list = parms.dMRI_ftype_list;
        case parms.fMRI_stype_list
          ftype_list = parms.fMRI_ftype_list;
      end;
      % use separate indir if qccont_flag = 1
      if parms.qccont_flag
        indir = sprintf('%s/%s',parms.RootDirs.qc,...
                        regexprep(ContainerDir,'PROC','QC'));
      else
        indir = ContainerPath;
      end;
      % load auto qcinfo
      auto_qcinfo = load_auto_qcinfo(indir,stype,parms);
      if isempty(auto_qcinfo), continue; end;
      nscans = length(auto_qcinfo);
      % load raw qcinfo
      raw_qcinfo = load_raw_qcinfo(indir,stype,parms,nscans);
      for j=1:nscans
        k = k + 1;
        qcinfo(k).VisitID = VisitID;
        qcinfo(k).fstem = auto_qcinfo(j).fstem;
        % add raw QC info
        qcinfo(k).nrev = 0;
        qcinfo(k).revdisp = 0;
        if ~isempty(raw_qcinfo)
          qcinfo(k).QC = 1;
          QC_last = [];
          for r=1:parms.nrev
            reviewerID = parms.reviewerIDs{r};
            if ~isempty(reviewerID)
              field_infix = ['_' reviewerID];
            else
              field_infix = '';
            end;
            fname = ['QC' field_infix];
            qcinfo(k).(fname) = mmil_getfield(raw_qcinfo(j),fname);
            if ~isempty(qcinfo(k).(fname))
              qcinfo(k).nrev = qcinfo(k).nrev + 1;
              QC = qcinfo(k).(fname);
              if ~QC, qcinfo(k).QC = 0; end;
              if ~isempty(QC_last) && QC~=QC_last
                qcinfo(k).revdisp = 1;
              end;
              QC_last = QC;
            end;            
            fname = ['notes' field_infix];
            qcinfo(k).(fname) = mmil_getfield(raw_qcinfo(j),fname);
            if ~isempty(qcinfo(k).(fname))
              qcinfo(k).(fname)  = regexprep(qcinfo(k).(fname),',',';');
            end;
            % add qualitative scores
            for f=1:length(ftype_list)
              ftype = ftype_list{f};
              fname = (['QU' field_infix '_' ftype]);
              qcinfo(k).(fname) = mmil_getfield(raw_qcinfo(j),fname);
            end;
          end;
        else
          qcinfo(k).QC = [];
          for r=1:parms.nrev
            reviewerID = parms.reviewerIDs{r};
            if ~isempty(reviewerID)
              field_infix = ['_' reviewerID];
            else
              field_infix = '';
            end;
            qcinfo(k).(['QC' field_infix]) = [];
            qcinfo(k).(['notes' field_infix]) = [];
            % add qualitative scores
            for f=1:length(ftype_list)
              ftype = ftype_list{f};
              fname = (['QU' field_infix '_' ftype]);
              qcinfo(k).(fname) = [];
            end;
          end;
        end;
        % copy information from auto_qcinfo
        qcinfo = copy_auto_qcinfo(auto_qcinfo,qcinfo,j,k,parms);
      end;
    end;
    % consolidate duplicate entries (fieldmap for and rev)
    if length(qcinfo)>1
      qcinfo = check_duplicate_qcinfo(qcinfo,parms);
    end;
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qcinfo = load_auto_qcinfo(indir,stype,parms)
  qctype = 'auto';
  qcinfo = []; auto_qcinfo = [];
  fname_qcinfo = sprintf('%s/%s_%s_qcinfo.mat',indir,stype,qctype);
  if ~exist(fname_qcinfo,'file')
    if parms.verbose
      fprintf('%s: WARNING: %s %s qcinfo file for %s not found\n',...
        mfilename,stype,qctype,indir);
    end;
    return;
  end;
  try
    load(fname_qcinfo);
  catch me
    if parms.verbose
      fprintf('%s: WARNING: %s %s qcinfo file %s unreadable\n',...
        mfilename,stype,qctype,fname_qcinfo);
    end;
    return;
  end;
  if ~isempty(auto_qcinfo), qcinfo = auto_qcinfo; end;
  if parms.verbose && isempty(qcinfo)
    fprintf('%s: WARNING: %s %s qcinfo for %s is empty\n',...
      mfilename,stype,qctype,indir);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qcinfo = check_duplicate_qcinfo(qcinfo,parms)
  nscans = length(qcinfo);
  SeUIDs = {qcinfo.SeriesInstanceUID};
  if ~iscellstr(SeUIDs), return; end;
  [tmp,ind] = unique(SeUIDs,'first');
  if length(ind)~=nscans
    for i=1:length(ind)
      j = ind(i);
      k = setdiff(find(strcmp(SeUIDs{j},SeUIDs)),j);
      if isempty(k), continue; end;
      % check QC values for scans with identical Series UID
      for r=1:parms.nrev
        reviewerID = parms.reviewerIDs{r};
        if ~isempty(reviewerID)
          field_infix = ['_' reviewerID];
        else
          field_infix = '';
        end;
        % set QC = 1 if both scans are good
        tag = ['QC' field_infix];
        if ~isempty(qcinfo(j).(tag)) &&...
           ~isempty(qcinfo(k).(tag))
          qcinfo(j).(tag) = (qcinfo(j).(tag)==1 &&...
                             qcinfo(k).(tag)==1);
        else
          qcinfo(j).(tag) = [];
        end;
        % concatenate notes
        tag = ['notes' field_infix];
        if ~isempty(qcinfo(j).(tag)) && ...
           ~isempty(qcinfo(k).(tag))
          qcinfo(j).(tag) = [qcinfo(j).(tag) '; '...
                             qcinfo(k).(tag)];
        elseif ~isempty(qcinfo(k).(tag))
          qcinfo(j).(tag) = qcinfo(k).(tag);
        end;
      end;
    end;
    qcinfo = qcinfo(ind);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qcinfo = load_raw_qcinfo(indir,stype,parms,nscans)
  qctype = 'raw';
  qcinfo = [];
  for r=1:parms.nrev
    reviewerID = parms.reviewerIDs{r};
    if ~isempty(reviewerID)
      field_infix = ['_' reviewerID];
    else
      field_infix = '';
    end;
    outfix = ['raw' field_infix];
    outfix = [outfix '_qcinfo'];
    fname_qcinfo = sprintf('%s/%s_%s.mat',indir,stype,outfix);
    if ~exist(fname_qcinfo,'file')
      if parms.verbose
        fprintf('%s: WARNING: %s %s file for %s not found\n',...
          mfilename,stype,regexprep(outfix,'_',' '),indir);
      end;
      continue;
    end;
    raw_qcinfo = [];
    load(fname_qcinfo);
    if isempty(raw_qcinfo)
      if parms.verbose
        fprintf('%s: WARNING: %s %s for %s is empty (deleting)\n',...
          mfilename,stype,regexprep(outfix,'_',' '),indir);
      end;
      delete(fname_qcinfo);
      continue;
    elseif length(raw_qcinfo) ~= nscans
      % compare length to nscans, delete file if mismatch
      if parms.verbose
        fprintf('%s: WARNING: mismatch between %s auto qcinfo and %s for %s (deleting %s)\n',...
          mfilename,stype,regexprep(outfix,'_',' '),indir,fname_qcinfo);
      end;
      raw_qcinfo = [];
      delete(fname_qcinfo);
      continue;
    end;
    
    if ~isempty(raw_qcinfo)
      % get qualitative scores
      ftype_list = {};
      if isfield(raw_qcinfo,'QU')
        switch stype
          case parms.sMRI_stype_list
            ftype_list = parms.sMRI_ftype_list;
          case parms.dMRI_stype_list
            ftype_list = parms.dMRI_ftype_list;
          case parms.fMRI_stype_list
            ftype_list = parms.fMRI_ftype_list;
        end;
      end; 
    
      for j=1:nscans
        QU = {};
        for f=1:length(ftype_list)
          ftype = ftype_list{f};
          g = find(strcmp(ftype,raw_qcinfo(j).ftype_list));
          if ~isempty(g)
            qu_score = raw_qcinfo(j).QU{g};
            if ischar(qu_score), qu_score = str2num(qu_score); end;
            QU{f} = qu_score;
          end;
        end;
        qcinfo(j).(['QC' field_infix]) = raw_qcinfo(j).QC;
        qcinfo(j).(['notes' field_infix]) = raw_qcinfo(j).notes;
        % add QU scores
        if ~isempty(QU)
          for f=1:length(ftype_list)
            ftype = ftype_list{f};
            qcinfo(j).(['QU' field_infix '_' ftype]) = QU{f};
          end;
        end;
        
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = copy_auto_qcinfo(A,B,a,b,parms)

  % stype-independent metrics
  B(b).fstem = A(a).fstem;
  B(b).stype = A(a).stype;
  B(b).StudyInstanceUID = A(a).seriesinfo.StudyInstanceUID;
  B(b).SeriesInstanceUID = A(a).seriesinfo.SeriesInstanceUID;
  B(b).nvox = A(a).imginfo.volsz(1:3);
  B(b).nreps = A(a).imginfo.nreps;
  B(b).orient = A(a).imginfo.orient;
  B(b).voxvol = A(a).imginfo.voxvol;
  B(b).brainvol = A(a).brain_vol;

  % stype-specific metrics
  switch B(b).stype
    case {'MPR','FLASHhi','XetaT2'}
      % SNR inside brain
      if ~isempty(A(a).brain_SNR)
        B(b).brain_mean = A(a).brain_SNR.mean;
        B(b).brain_std = A(a).brain_SNR.std;
        B(b).brain_SNR = A(a).brain_SNR.SNR;
        B(b).brain_min = A(a).brain_SNR.min;
        B(b).brain_max = A(a).brain_SNR.max;
        B(b).brain_median = A(a).brain_SNR.median;
      end;
      % SNR outside brain
      if ~isempty(A(a).nonbrain_SNR)
        B(b).nonbrain_mean = A(a).nonbrain_SNR.mean;
        B(b).nonbrain_std = A(a).nonbrain_SNR.std;
        B(b).nonbrain_SNR = A(a).nonbrain_SNR.SNR;
      end;
    case 'BOLD'
      % SNR inside brain
      if ~isempty(A(a).brain_SNR)
        B(b).brain_mean = A(a).brain_SNR.mean;
        B(b).brain_std = A(a).brain_SNR.std;
        B(b).brain_SNR = A(a).brain_SNR.SNR;
        B(b).brain_min = A(a).brain_SNR.min;
        B(b).brain_max = A(a).brain_SNR.max;
        B(b).brain_median = A(a).brain_SNR.median;
      end;
      % TR
      B(b).TR = A(a).scaninfo.TR/1000;
      % write average frame-to-frame motion
      if ~isempty(A(a).motion)
        B(b).mean_motion = A(a).motion.motion_stats.mean_motion;
        B(b).mean_motion_nody = A(a).motion.motion_stats_nody.mean_motion;
        B(b).mean_trans = A(a).motion.motion_stats.mean_trans;
        B(b).mean_trans_nody = A(a).motion.motion_stats_nody.mean_trans;
        B(b).mean_rot = A(a).motion.motion_stats.mean_rot;
        B(b).max_dx = A(a).motion.motion_stats.max_dx;
        B(b).max_dy = A(a).motion.motion_stats.max_dy;
        B(b).max_dz = A(a).motion.motion_stats.max_dz;
        B(b).max_rx = A(a).motion.motion_stats.max_rx;
        B(b).max_ry = A(a).motion.motion_stats.max_ry;
        B(b).max_rz = A(a).motion.motion_stats.max_rz;
        for i=1:parms.nthresh
          B(b).(sprintf('subthresh_0%0.0f',10*parms.thresholds(i))) = ...
            A(a).motion.motion_stats.subthresh_nvols(i)*B(b).TR;
          B(b).(sprintf('subthresh_0%0.0f_nody',10*parms.thresholds(i))) = ...
            A(a).motion.motion_stats_nody.subthresh_nvols(i)*B(b).TR;
        end;
      end;
      % write tSNR inside brain
      if ~isempty(A(a).brain_tSNR)
        B(b).brain_tSNR_mean = A(a).brain_tSNR.mean;
        B(b).brain_tSNR_median = A(a).brain_tSNR.median;
        B(b).brain_tSNR_std = A(a).brain_tSNR.std;
      end;
    case 'DTI'
      % SNR inside brain
      if ~isempty(A(a).brain_SNR)
        B(b).brain_mean = A(a).brain_SNR.mean;
        B(b).brain_std = A(a).brain_SNR.std;
        B(b).brain_SNR = A(a).brain_SNR.SNR;
        B(b).brain_min = A(a).brain_SNR.min;
        B(b).brain_max = A(a).brain_SNR.max;
        B(b).brain_median = A(a).brain_SNR.median;
      end;
      % TR
      B(b).TR = A(a).scaninfo.TR/1000;
      % write average frame-to-frame motion
      if ~isempty(A(a).motion)
        B(b).mean_motion = A(a).motion.motion_stats.mean_motion;
        B(b).mean_motion_nody = A(a).motion.motion_stats_nody.mean_motion;
        B(b).mean_trans = A(a).motion.motion_stats.mean_trans;
        B(b).mean_trans_nody = A(a).motion.motion_stats_nody.mean_trans;
        B(b).mean_rot = A(a).motion.motion_stats.mean_rot;
        B(b).max_dx = A(a).motion.motion_stats.max_dx;
        B(b).max_dy = A(a).motion.motion_stats.max_dy;
        B(b).max_dz = A(a).motion.motion_stats.max_dz;
        B(b).max_rx = A(a).motion.motion_stats.max_rx;
        B(b).max_ry = A(a).motion.motion_stats.max_ry;
        B(b).max_rz = A(a).motion.motion_stats.max_rz;
        for i=1:parms.nthresh
          B(b).(sprintf('subthresh_0%0.0f',10*parms.thresholds(i))) = ...
            A(a).motion.motion_stats.subthresh_nvols(i)*B(b).TR;
          B(b).(sprintf('subthresh_0%0.0f_nody',10*parms.thresholds(i))) = ...
            A(a).motion.motion_stats_nody.subthresh_nvols(i)*B(b).TR;
        end;
      end;
      % write number of bad slices
      if ~isempty(A(a).censor)
        B(b).nbad_frame_slices = A(a).censor.nbad_frame_slices;
        B(b).nbad_frames = A(a).censor.nbad_frames;
        B(b).nbad_slices = A(a).censor.nbad_slices;
        B(b).max_nbad_frames_per_slice = A(a).censor.max_nbad_frames_per_slice;
        B(b).max_nbad_slices_per_frame = A(a).censor.max_nbad_slices_per_frame;
      end;
      % write about average brain DTmeas (e.g. FA, MD, b0)
      if ~isempty(A(a).brain_DTmeas)
        meas_list = fieldnames(A(a).brain_DTmeas);
        for m=1:length(meas_list)
          meas = meas_list{m};
          B(b).([meas '_mean']) = A(a).brain_DTmeas.(meas).mean;
          B(b).([meas '_median']) = A(a).brain_DTmeas.(meas).median;
          B(b).([meas '_std']) = A(a).brain_DTmeas.(meas).std;
        end;
      end;
      % write about average brain DTerr
      if ~isempty(A(a).brain_DTerr)
        B(b).DTerr_mean = A(a).brain_DTerr.brain_rms_err_ratio_mean;
        B(b).DTerr_median = A(a).brain_DTerr.brain_rms_err_ratio_median;
        B(b).DTerr_std = A(a).brain_DTerr.brain_rms_err_ratio_std;
      end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qcinfo = concat_qcinfo(qcinfo,tmpinfo)
  if isempty(tmpinfo), return; end;
  if isempty(qcinfo)
    qcinfo = tmpinfo;
  else
    % reconcile fieldnames (add empty fields to match)
    [qcinfo,tmpinfo] = reconcile_fields(qcinfo,tmpinfo);
    qcinfo = cat(2,qcinfo,tmpinfo);    
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,B] = reconcile_fields(A,B)
  fnamesA = fieldnames(A);
  fnamesB = fieldnames(B);
  fnames = setdiff(fnamesB,fnamesA);
  for f=1:length(fnames)
    A(1).(fnames{f}) = [];
  end;
  fnames = setdiff(fnamesA,fnamesB);
  for f=1:length(fnames)
    B(1).(fnames{f}) = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

