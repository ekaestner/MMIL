function BOLD_MMIL_Combine_Fourier_Sessions(ContainerPaths,FSContainerPath,varargin)
%function BOLD_MMIL_Combine_Fourier_Sessions(ContainerPaths,FSContainerPath,[options])
%
% Purpose: Combine Fourier analysis results from multiple sessions
%   Assumes BOLD_MMIL_Fourier_Analysis has been run, including
%   painting to common surface
%
% Usage:
%  BOLD_MMIL_Combine_Fourier_Sessions(ContainerPaths,'key1', value1,...);
%  e.g. BOLD_MMIL_Combine_Fourier_Sessions(ContainerPaths,...
%         'snums',{[3,4],[3,4]},'revflags',{[0,1],[0,1]});
%
% Required Input:
%  ContainerPaths: Cell array containing full paths of BOLDPROC directories
%    for each session
%  FSContainerPath: Full path of FreeSurfer recon container
%   if sphere_flag = 1, FSContainerPath may be empty (will use fsaverage)
%   if sphere_flag = 1 and FSContainerPath is not empty, that recon
%     should be an "average" subject or resampled to ico
%
% Optional Paramters containing information for each session:
%  'snums': cell array containing a vector of scan numbers
%    for each session, to specify which scans to include for each session
%    if empty or omitted, will use all BOLD scans in each ContainerPath
%    {default = []}
%  'snums_valid': cell array containing a vector of scan numbers
%    for each session, to specify which scans were actually processed
%    if empty or omitted, will assume all BOLD scans in each ContainerPath
%    {default = []}
%  'revflags': cell array containing a vector of 0's or 1's for each
%    sesssion, to to specify whether to reverse phase for each scan
%    {default = []} (zero for all)
%  'infix': if empty, will look for files like 'BOLD1.mgz'
%    otherwise, input file will be sprintf('BOLD%d_%s.mgz',snum,infix)
%    example infix = 'corr_resBOLD'
%    {default = []}
%  'multi_flag': [0|1|2] indicates source of input analyses
%    0: use results from individual scans
%    1: use results from multiple scans within session
%    2: use results from multiple sessions
%    NOTE: if multi_flag>0, revflags, phase_offset, phase_offset_postrev
%     are ignored (assumed to be accounted for already)
%    {default = 0}
%
% Optional Parameters specifying how Fourier analysis was done:
%  'fstats_type': [0|1|2] how output Fourier components were scaled
%    0: raw, no scaling
%    1: scaled by sqrt(F-ratio)
%    2: scaled by significance values (-log10(p))
%    {default = 2}
%  'resamp_flag': [0|1] whether Fourier results were resampled to 1x1x1
%    before painting
%    {default = 0}
%  'smoothsteps': smoothing steps on surface after painting
%    {default = 0}
%  'sphere_flag': [0|1] whether data were resampled to icosohedral sphere
%    {default = 0}
%  'sphsmoothsteps': smoothing steps on spherical surface
%    {default = 0}
%
% Optional Parameters for view results:
%  'tksmooth': number of surface smoothing iterations for tksurfer script
%    (relevant only if paint_flag = 1)
%    {default = 2}
%  'fthresh': threshold for scaling amplitudes in tksurfer
%    {default = 0}
%  'fmid': mid point for scaling amplitudes in tksurfer
%    {default = 1.5}
%  'fslope': slope for scaling amplitudes in tksurfer
%    {default = 1.0}
%  'view': view of cortical surface view in tksurfer
%    ('med','ven','lat','pos','dor')
%    {default = 'med'}
%  'surf': surface to display in tksurfer
%    {default = 'inflated'}
%  'datatype': type of complex data (affects color scale in tksurfer)
%    ('polar','eccen')
%    {default = 'polar'}
%
% Other optional parameters:
%  'cxfstatsflag' [0|1] whether to calculate cross-scan complex f-stats
%    (in addition to amplitude f-stats)
%    {default = 1}
%  'phase_offset': fraction of cycles to subtract from phases to account for
%    hemodynamic delay; may be vector with a value for each session
%    {default = 0.1}
%  'phase_offset_postrev': fraction of cycles to subtract from phases
%      to account for stimulus delay; may be vector a value for each session
%    {default = 0}
%  'r_max': maximum eccentricity (applies to eccentricity mapping only)
%     may be vector with a value for each session
%     ignored if datatype = polar
%     {default = 12.5}
%  'r_min': minimum eccentricity (applies to eccentricity mapping only)
%     may be vector with a value for each session
%     ignored if datatype = polar
%    {default = []}
%  'r_min_factor': value multiplied by r_max to get r_min (if empty)
%     ignored if datatype = polar
%    {default = 0.02}
%  'logtrans_flag': whether log transform was used for eccentricity mapping
%    may be vector with a value for each session
%     ignored if datatype = polar
%    {default = 0}
%  'OutContainerPath': output ContainerPath
%    if empty, will use first entry of ContainerPaths
%    {default = []}
%  'outdir': output directory
%     If not a full path, will create a subdir in OutContainerPath
%    {default =  'BOLD_multisess_analysis'}
%  'outstem': output file stem
%    {default =  'BOLD_multisess'}
%  'forceflag': [0|1] whether to run calculations even if output files exist
%    {default = 0}
%
% NOTE: r_max, r_min, r_min_factor, and logtrans_flag are for adjusting
%   phases for varying maximum eccentricity across sessions
%   if all of these are identical for all sessions are used,
%   these values are ignored
%
% Created:  07/05/09 by Don Hagler
% Rcnt Mod: 09/10/12 by Don Hagler
% Last Mod: 09/15/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'snums',[],[],...
  'snums_valid',[],[],...
  'revflags',[],[],...
  'infix',[],[],...
  'multi_flag',0,[0:2],...
  'fstats_type',2,[0:2],...
  'resamp_flag',false,[false true],...
  'smoothsteps',0,[0,1000],...
  'sphsmoothsteps',0,[0,1000],...
  'sphere_flag',false,[false true],...
... % parameters for viewing results
  'tksmooth',2,[0,1000],...
  'fthresh',0,[0,100],...
  'fmid',1.5,[0.1,100],...
  'fslope',1.0,[0.1,100],...
  'view','med',{'lat','med','ven','pos','dor'},...
  'datatype','polar',{'polar','eccen'},...
  'surf','inflated',[],...
...
  'cxfstatsflag',true,[false true],...
  'phase_offset',0.1,[-1,1],...
  'phase_offset_postrev',0,[-1,1],...
  'r_max',12.5,[],...
  'r_min',[],[],...
  'r_min_factor',0.02,[],...
  'logtrans_flag',[],[],...
  'OutContainerPath',[],[],...
  'outdir','BOLD_multisess_analysis',[],...
  'outstem','BOLD_multisess',[],...
  'forceflag',false,[false true],...
...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'suffix_list',{'r','i'},{'r','i'},...
  'out_ext','.mgh',{'.mgh','.mgz'},...
  'fnamestem','BOLD',[],...
...
  'ContainerPaths',ContainerPaths,[],...
  'FSContainerPath',FSContainerPath,[],...
...
  'fsaverage','fsaverage',[],... % used if sphere_flag = 1
  'fsavg_rootdir',[getenv('FREESURFER_HOME') '/subjects'],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check ContainerPaths, snums, and revflags
if ~iscell(parms.ContainerPaths)
  error('ContainerPaths must be a cell array');
end;
nsess = length(parms.ContainerPaths);
if isempty(parms.snums)
  parms.snums = cell(nsess,1);
end;
if isempty(parms.snums_valid)
  parms.snums_valid = cell(nsess,1);
end;
if isempty(parms.revflags)
  parms.revflags = cell(nsess,1);
end;
if ~iscell(parms.snums)
  error('snums must be a cell array with a vector for each ContainerPath');
end;
if ~iscell(parms.revflags)
  error('revflags must be a cell array with a vector for each ContainerPath');
end;

if length(parms.phase_offset)==1
  parms.phase_offset = ones(nsess,1)*parms.phase_offset;
end;
if length(parms.phase_offset)~=nsess
  error('number of phase_offset values must match number of ContainerPaths');
end;

if length(parms.phase_offset_postrev)==1
  parms.phase_offset_postrev = ones(nsess,1)*parms.phase_offset_postrev;
end;
if length(parms.phase_offset_postrev)~=nsess
  error('number of phase_offset_postrev values must match number of ContainerPaths');
end;

parms.eccen_flag = strcmp(parms.datatype,'eccen');
if parms.eccen_flag
  if length(parms.r_max)==1
    parms.r_max = ones(nsess,1)*parms.r_max;
  end;
  if length(parms.r_max)~=nsess
    error('number of r_max values must match number of ContainerPaths');
  end;

  if isempty(parms.r_min)
    parms.r_min = parms.r_min_factor*parms.r_max;  
  elseif length(parms.r_min)==1
    parms.r_min = ones(nsess,1)*parms.r_min;
  elseif length(parms.r_min)~=nsess
    error('number of r_min values must match number of ContainerPaths');
  end;

  if isempty(parms.logtrans_flag)
    parms.logtrans_flag = zeros(nsess,1);
  elseif length(parms.logtrans_flag)==1
    parms.logtrans_flag = ones(nsess,1)*parms.logtrans_flag;
  elseif length(parms.logtrans_flag)~=nsess
    error('number of logtrans_flag values must match number of ContainerPaths');
  end;
else
  parms.r_max = [];
  parms.r_min = [];
  parms.logtrans_flag = [];
end;

for i=1:nsess
  ContainerPath = parms.ContainerPaths{i};
  snums = parms.snums{i};
  snums_valid = parms.snums_valid{i};
  revflags = parms.revflags{i};
  phase_offset = parms.phase_offset(i);
  phase_offset_postrev = parms.phase_offset_postrev(i);
  if parms.eccen_flag
    r_max = parms.r_max(i);
    r_min = parms.r_min(i);
    logtrans_flag = parms.logtrans_flag(i);
  end;
  
  if ~exist(ContainerPath)
    error('ContainerPath %s not found',ContainerPath);
  end;

  % get number of scans, etc.
  [ScanInfo,SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(ContainerPath,...
    'snums',snums_valid,'fnamestem',parms.fnamestem);
  if errcode || isempty(ScanInfo)
    error('no BOLD scans found in %s',ContainerPath);
  end;

  if isempty(snums)
    snums = SessInfo.snums_valid;
  end;
  nscans = length(snums);
  if isempty(revflags), revflags = zeros(1,nscans); end;
  if length(revflags) ~= nscans
    error('length of revflags (%d) does not match number of scans (%d) for %s',...
      length(parms.revflags),nscans,ContainerPath);
  end;
  parms.snums{i} = snums;
  parms.revflags{i} = revflags;
  parms.phase_offsets{i} = phase_offset * ones(1,nscans);
  parms.phase_offset_postrevs{i} = phase_offset_postrev * ones(1,nscans);
  if parms.eccen_flag
    parms.r_max_vecs{i} = r_max * ones(1,nscans);
    parms.r_min_vecs{i} = r_min * ones(1,nscans);
    parms.logtrans_flags{i} = logtrans_flag * ones(1,nscans);
  else
    parms.r_max_vecs = [];
    parms.r_min_vecs = [];
    parms.logtrans_flags = [];
  end;
end;
  
% check FSContainerPath
if isempty(parms.FSContainerPath)
  if parms.sphere_flag
    parms.FSContainerPath = [parms.fsavg_rootdir '/' parms.fsaverage];
  else
    error('FSContainerPath is empty');
  end;
end;
if ~exist(parms.FSContainerPath,'dir')
  error('FSContainerPath %s not found',parms.FSContainerPath);
end;
% get parms.FSRootDir (SUBJECTS_DIR) from FSContainerPath
[parms.FSRootDir,parms.FSDir,ext] = fileparts(parms.FSContainerPath);
parms.FSDir = [parms.FSDir ext]; % in case of .'s and such
setenv('SUBJECTS_DIR',parms.FSRootDir);

if isempty(parms.OutContainerPath)
  parms.OutContainerPath = parms.ContainerPaths{1};
end;
% check if outdir is full path, create output dir
if mmil_isrelative(parms.outdir)
  parms.outdir = [parms.OutContainerPath '/' parms.outdir];
else
  [parms.OutContainerPath,tmp] = fileparts(parms.outdir);
end;
mmil_mkdir(parms.outdir);

switch parms.fstats_type
  case 0
    parms.fstats_infix = 'fstats_raw';
  case 1
    parms.fstats_infix = 'fstats_ratio';
  case 2
    parms.fstats_infix = 'fstats_pval';
end;

parms.full_suffix_list = {};
for j=1:length(parms.suffix_list)
  suffix = parms.suffix_list{j};
  if parms.resamp_flag
    parms.full_suffix_list{j} = ['resT1_' suffix];
  else
    parms.full_suffix_list{j} = suffix;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross-scan average of fourier series

for h=1:length(parms.hemilist)
  % loop over hemis
  hemi = parms.hemilist{h};

  fname_log = [parms.outdir '/' parms.outstem '-' hemi '.log'];
  fid = fopen(fname_log,'wt');
  if fid<0
    error('failed to open %s for writing',fname_log);
  end;
  
  % compile list of files
  fnamelist = [];
  revflags = [];
  phase_offsets = [];
  phase_offset_postrevs = [];
  r_max_vec = [];
  r_min_vec = [];
  logtrans_flags = [];
  k=1;
  for i=1:nsess
    ContainerPath = parms.ContainerPaths{i};
    snums = parms.snums{i};
    snums_valid = parms.snums_valid{i};
    if parms.multi_flag
      if parms.multi_flag==2
        adir = sprintf('%s_multisess_%s_analysis',...
          parms.fnamestem,parms.datatype);
        stem = sprintf('%s_multisess_%s_avg',...
          parms.fnamestem,parms.datatype);
      else
        stem = []; adir = [];
      end;
      tags = {'infix','fstats_infix','fnamestem'};
      args = mmil_parms2args(parms,tags);
      [full_stem,errcode] = BOLD_MMIL_Set_Fourier_Stem(ContainerPath,...
        args{:},'snums_valid',snums_valid,...
        'snums',snums,'stem',stem,'adir',adir);
      if errcode, return; end;
      full_stem_list = {full_stem};
      % assume that this is a multi-session average
      %   and that revflags and phase offsets have already been applied
      tmp_revflags = 0;
      tmp_phase_offsets = 0;
      tmp_phase_offset_postrevs = 0;
      if ~isempty(parms.r_max_vecs)
        tmp_r_max_vec = max(parms.r_max_vecs{i});
      end;
      if ~isempty(parms.r_min_vecs)
        tmp_r_min_vec = min(parms.r_min_vecs{i});
      end;
      if ~isempty(parms.logtrans_flags)
        tmp_logtrans_flags = min(parms.logtrans_flags{i});
      end;
    else
      full_stem_list = cell(length(snums),1);
      for s=1:length(snums)
        tags = {'infix','fstats_infix','fnamestem'};
        args = mmil_parms2args(parms,tags);
        [full_stem_list{s},errcode] = BOLD_MMIL_Set_Fourier_Stem(...
          ContainerPath,args{:},'snums',snums(s));
      end;
      tmp_revflags = parms.revflags{i};
      tmp_phase_offsets = parms.phase_offsets{i};
      tmp_phase_offset_postrevs = parms.phase_offset_postrevs{i};
      if parms.eccen_flag
        tmp_r_max_vec = parms.r_max_vecs{i};
        tmp_r_min_vec = parms.r_min_vecs{i};
        tmp_logtrans_flags = parms.logtrans_flags{i};      
      end;
    end;

    revflags = [revflags,tmp_revflags];
    phase_offsets = [phase_offsets,tmp_phase_offsets];
    phase_offset_postrevs = [phase_offset_postrevs,tmp_phase_offset_postrevs];
    if parms.eccen_flag
      r_max_vec = [r_max_vec,tmp_r_max_vec];
      r_min_vec = [r_min_vec,tmp_r_min_vec];
      logtrans_flags = [logtrans_flags,tmp_logtrans_flags];
    else
      r_max_vec = [];
      r_min_vec = [];
      logtrans_flags = [];
    end;

    for s=1:length(full_stem_list)
      for j=1:length(parms.full_suffix_list)
        full_suffix = parms.full_suffix_list{j};
        tmp_fname = [full_stem_list{s} '_' full_suffix];
        if parms.smoothsteps
          tmp_fname = sprintf('%s-sm%d',tmp_fname,parms.smoothsteps);
        end;
        if parms.sphere_flag
          tmp_fname = sprintf('%s-sphere-sm%d',tmp_fname,parms.sphsmoothsteps);
        end;
        tmp_fname = [tmp_fname '-' hemi parms.out_ext];
        if ~exist(tmp_fname,'file')
          error('file %s not found',tmp_fname);
        end;
        fnamelist{k,j} = tmp_fname;
        fprintf(fid,'%s\n',tmp_fname);
      end;
      k = k + 1;
    end;
  end;
  fclose(fid);

  % decide whether we need to (re)do average
  calc_avg_flag = 0;
  for j=1:length(parms.full_suffix_list)
    fname_out = sprintf('%s/%s_avg_%s_%s-%s.mgh',...
      parms.outdir,parms.outstem,parms.fstats_infix,...
      parms.full_suffix_list{j},hemi);
    if ~exist(fname_out,'file')
      calc_avg_flag = 1;
      break;
    end;
  end;
  if parms.cxfstatsflag
    if parms.resamp_flag
      fname_out = sprintf('%s/%s_cxfstats_resT1-%s.mgh',...
        parms.outdir,parms.outstem,hemi);
    else
      fname_out = sprintf('%s/%s_cxfstats-%s.mgh',...
        parms.outdir,parms.outstem,hemi);
    end;
    if ~exist(fname_out,'file'), calc_avg_flag = 1; end;
  end;

  if calc_avg_flag || parms.forceflag
    % average across scans and save output
    fprintf('%s: calculating cross-scan average for hemi %s...\n',...
      mfilename,hemi);
    results = fs_fourier_avg(fnamelist,...
      'revflags',revflags,...
      'phase_offset',phase_offsets,...
      'cxfstatsflag',parms.cxfstatsflag,...
      'phase_offset_postrev',phase_offset_postrevs,...
      'r_max_vec',r_max_vec,...
      'r_max',parms.r_max,...
      'r_min_vec',r_min_vec,...
      'r_min_factor',parms.r_min_factor,...
      'logtrans_flags',logtrans_flags);
    for j=1:length(parms.full_suffix_list)
      data = getfield(results,sprintf('mean_%s',parms.suffix_list{j}));
      fname_out = sprintf('%s/%s_avg_%s_%s-%s.mgh',...
        parms.outdir,parms.outstem,parms.fstats_infix,...
        parms.full_suffix_list{j},hemi);
      fprintf('%s: saving average to %s...\n',...
        mfilename,fname_out);
      fs_save_mgh(data,fname_out,results.M);
    end;
    % save cross-scan complex fstats
    if parms.cxfstatsflag
      data = getfield(results,sprintf('Fstat',suffix));
      if parms.resamp_flag
        fname_out = sprintf('%s/%s_cxfstats_resT1-%s.mgh',...
          parms.outdir,parms.outstem,hemi);
      else
        fname_out = sprintf('%s/%s_cxfstats-%s.mgh',...
          parms.outdir,parms.outstem,hemi);
      end;
      fprintf('%s: saving complex fstats to %s...\n',...
        mfilename,fname_out);
      fs_save_mgh(data,fname_out,results.M);
    end;
    clear results;
  end;
end;

%% todo: general function for writing fs_surfmgh.csh scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate csh scripts to view results

%% todo: script to view cxfstats

hemi = parms.hemilist{1};
fstem = sprintf('%s_avg_%s',parms.outstem,parms.fstats_infix);
if parms.resamp_flag
  fstem = [fstem '_resT1'];
end;
indir = parms.outdir;
fname_view = sprintf('%s/view_BOLD_Fourier_%s_results_multisess.csh',...
  parms.OutContainerPath,parms.datatype);
revflag = 0;
[fpath,indir] = fileparts(indir);
fid = fopen(fname_view,'wt');
if fid==-1
  error('failed to open %s for writing',fname_view);
end;
fprintf(fid,'#!/bin/csh -f\n');
fprintf(fid,'# automatically generated script for viewing BOLD Fourier results\n');
fprintf(fid,'\n');
fprintf(fid,'if ($#argv == 0) then\n');
fprintf(fid,'  set hemi = %s\n',hemi);
fprintf(fid,'else\n');
fprintf(fid,'  set hemi = $argv[1]\n');
fprintf(fid,'endif\n');
fprintf(fid,'\n');
fprintf(fid,'if ($#argv > 1) then\n');
fprintf(fid,'  set save_flag = $argv[2]\n',hemi);
fprintf(fid,'else\n');
fprintf(fid,'  set save_flag = 0\n');
fprintf(fid,'endif\n');
fprintf(fid,'\n');
fprintf(fid,'setenv SUBJECTS_DIR %s\n',parms.FSRootDir);
fprintf(fid,'set subj = %s\n',parms.FSDir);
fprintf(fid,'set fstem = %s\n',fstem);
fprintf(fid,'set indir = %s\n',indir);
fprintf(fid,'\n');
fprintf(fid,'set realfile = $fstem''_r''\n');
fprintf(fid,'set imagfile = $fstem''_i''\n');
fprintf(fid,'\n');
fprintf(fid,'set smooth = %d\n',parms.tksmooth);
fprintf(fid,'set fthresh = %0.1f\n',parms.fthresh);
fprintf(fid,'set fmid = %0.1f\n',parms.fmid);
fprintf(fid,'set fslope = %0.1f\n',parms.fslope);
fprintf(fid,'set view = %s\n',parms.view);
fprintf(fid,'set surf = %s\n',parms.surf);
fprintf(fid,'\n');
fprintf(fid,'if $save_flag then\n');
fprintf(fid,'  set savestr = "-savetiff -outstem $fstem -offscreen"\n');
fprintf(fid,'else\n');
fprintf(fid,'  set savestr = " "\n');
fprintf(fid,'endif\n');
fprintf(fid,'\n');
fprintf(fid,'fs_surfmgh.csh $subj $realfile $hemi \\\n');
fprintf(fid,'  -imagfile $imagfile \\\n');
if strcmp(parms.datatype,'polar')
  fprintf(fid,'  -polar');
else
  fprintf(fid,'  -eccen');
end;
if revflag
  fprintf(fid,' -revphase \\\n');
else
  fprintf(fid,' \\\n');
end;
fprintf(fid,'  -surf $surf \\\n');
fprintf(fid,'  -indir $indir -smooth $smooth \\\n');
fprintf(fid,'  -fmid $fmid -view $view $savestr\n');
fprintf(fid,'\n');
fprintf(fid,'rm tempsurf.tcl.*\n');
fprintf(fid,'\n');
fclose(fid);

% change permissions so group can write and execute
cmd = sprintf('chmod ug+rwx %s',fname_view);
[status,result] = unix(cmd);
if status
  fprintf('%s: WARNING: failed to set permissions for %s:\n%s\n',...
    mfilename,fname_view,result);
end;
