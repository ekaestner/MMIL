function errcode = BOLD_MMIL_GLM_Analysis(ContainerPath,stim_fnames,varargin)
%function errcode = BOLD_MMIL_GLM_Analysis(ContainerPath,stim_fnames,[options])
%
% Purpose: Run GLM analysis on fMRI data
%
% Usage:
%  BOLD_MMIL_GLM_Analysis(ContainerPath,stim_fnames,'key1', value1,...);
%  e.g. BOLD_MMIL_GLM_Analysis(ContainerPath,stim_fnames,...
%         'FSPath',FSPath,'paint_flag',1,...
%         'snums',[3,4,5,6],'skipTRs',2);
%
% Required Input:
%  ContainerPath: Full path of BOLDPROC directory containing BOLD scans
%  stim_fnames: cell array of stimulus time course 1D files
%    one for each non-baseline condition for each snum
%    containing the names of 1D files
%    {default = []}
%
% Optional Paramters:
%  'mc_flag': [0|1] whether within-scan motion correction was done
%    Will use motion.1D files as nuisance regressors
%    {default = 1}
%  'mc_inter_flag': [0|1] whether between-scan motion correction was done
%    Allows for a single reference scan to be used for registration to FS
%    {default = 1}
%  'FSPath': Full path of FreeSurfer directory containing cortical surface recon
%    (required if paint_flag > 0)
%    {default = []}
%  'snums': vector of scan numbers on which to run GLM analysis
%    if empty or omitted, will use all BOLD scans in ContainerPath
%    {default = []}
%  'infix': BOLD file name infix (e.g. '', 'corr', 'corr_resBOLD')
%    {default = []}
%  'num_conds': vector of numbers of conditions per scan
%    if empty, will assume one condition per scan
%    {default = []}
%  'stim_labels': string or cell array of strings with stimulus condition names
%    If empty, will use labels such as "cond1", "cond2", etc.
%    {default = []}
%  'contrasts_flag': [0|1] calculate glt contrasts between each condition
%    {default = 0}
%  'iresp_flag': [0|1] output impulse response functions
%    {default = 0}
%  'concat_flag': [0|1] concatenate across multiple runs
%    Each scan will also be analyzed individually
%    {default = 1}
%  'skipTRs': number of TRs at beginning of each run to ignore
%    {default = 0}
%  'minlag': number of TRs for minimum "lag" between stimulus and response
%    {default = 0}
%  'maxlag': number of TRs for maximum "lag" between stimulus and response
%    {default = 4}
%  'paint_flag': [0|1|2] whether to paint results to surface
%     if 0, only do GLM analysis, no painting to surface
%     if 1, paint results to surface after averaging across scans
%       (requires motion correction and register.dat file for first scan)
%     if 2, paint results to surface before averaging across scans
%       (requires register.dat file for each scan if not motion corrected)
%    {default = 0}
%  'projdist': dist (mm) to project along normal when painting
%    { default = 1 }
%  'projfrac': fractional dist to project along normal when painting
%    { default = 0.5 }
%  'projfrac_flag': [0|1] whether to use projfrac (1) or projdist (0)
%    { default = 0 }
%  'projfrac_avg': vector of [min max del] for averaging multiple samples
%    If empty, use projfrac instead if projfrac_flag=1
%    {default = []}
%  'projdist_avg': vector of [min max del] for averaging multiple samples
%    If empty, use projdist instead
%    {default = []}
%  'resamp_flag': [0|1] whether to resample GLM results to 1x1x1
%    before painting
%    {default = 0}
%  'force_repaint_flag': [0|1] whether to repaint even if output files exist
%    (do not redo GLM calculations)
%    {default = 0}
%  'sphere_flag': [0|1] whether to sample to icosohedral sphere
%    {default = 0}
%  'smoothsteps': smoothing steps on surface after painting
%    {default = 0}
%  'sphsmoothsteps': smoothing steps on spherical surface
%    {default = 0}
%  'pthresh': probability threshold applied to f-stats for each
%     condition
%    {default = 0.01}
%  'tksmooth': number of surface smoothing iterations for tksurfer script
%    (relevant only if paint_flag = 1)
%    {default = 2}
%  'out_ext': output file extension ('.mgh' or '.mgz')
%    if empty, will use input file extension
%    {default = []}
%  'forceflag': [0|1] whether to run calculations even if output files exist
%    {default = 0}
%
% Created:  09/01/08 by Don Hagler
% Prev Mod: 12/14/16 by Don Hagler
% Last Mod: 08/30/17 by Don Hagler
%

%% todo: accept csv file as input to specify stim_fnames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

%     if 2, resample fseries to T1 before averaging across scans
%           (won't work because of memory limitations)
%  'resamp_flag',0,[0,1,2],...

errcode = [];
if (~mmil_check_nargs(nargin,2)) return; end;
parms = mmil_args2parms(varargin, { ...
  'mc_flag',true,[false true],...
  'mc_inter_flag',true,[false true],...
  'FSPath',[],[],...
  'snums',[],[],...
  'snums_valid',[],[],...
  'infix',[],[],...
  'num_conds',[],[],...
  'stim_labels',[],[],...
  'contrasts_flag',false,[false true],...
  'conds_contrast',[],[],...
  'iresp_flag',false,[false true],...
  'concat_flag',true,[false true],...
  'skipTRs',0,[0,10],...
  'minlag',0,[0,10],...
  'maxlag',4,[0,30],...
  'paint_flag',0,[0,2],...
  'projdist',1,[0,10],...
  'projfrac',0.5,[0,2],...
  'projfrac_flag',false,[false true],...
  'projdist_avg',[],[],...
  'projfrac_avg',[],[],...
  'paint_surf','white',[],...
  'resamp_flag',false,[false true],...
  'pthresh',0.01,[0,1],...
  'force_repaint_flag',false,[false true],...
  'smoothsteps',0,[0,1000],...
  'sphsmoothsteps',0,[0,1000],...
  'sphere_flag',false,[false true],...
  'surfname','white',[],...
  'mask_midbrain_flag',false,[false true],...
... % parameters for viewing results
  'tksmooth',2,[0,1000],...
  'fthresh',0.1,[0,10],...
  'fmid',1,[0.1,10],...
  'fslope',1.0,[0.1,10],...
  'view','med',{'lat','med','ven','pos'},...
  'colscale',1,[0:8],...
...
  'revflag',false,[false true],... % only applys for pep and ipp
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'out_ext',[],{'.mgh','.mgz'},...
  'forceflag',false,[false true],...
...
  'fnamestem','BOLD',[],...
});
parms.overwrite_flag = (parms.forceflag || parms.force_repaint_flag);
paint_tags = {'regfile' 'projfrac_flag' 'projdist' 'projfrac'...
             'projdist_avg' 'projfrac_avg' 'mask_midbrain_flag' ...
             'subjdir' 'surfname' 'overwrite_flag'};
paint_args = mmil_parms2args(parms,paint_tags);

fstem_results = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ScanInfo,SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(ContainerPath,...
  'snums',parms.snums_valid,'fnamestem',parms.fnamestem);
if errcode || isempty(ScanInfo), return; end;

if parms.paint_flag || parms.resamp_flag
  if isempty(parms.FSPath)
    error('FSPath must not be empty if paint_flag > 0');
  end;
  if ~exist(parms.FSPath,'dir')
    error('FSPath %s not found',parms.FSPath);
  end;
  % get parms.subjdir (SUBJECTS_DIR) from FSPath
  [parms.subjdir,parms.subj,ext] = fileparts(parms.FSPath);
  parms.subj = [parms.subj ext]; % in case of .'s and such
  setenv('SUBJECTS_DIR',parms.subjdir);
  fprintf('%s: using scan %d as reference\n',mfilename,SessInfo.regT1_ref);
  fstem_ref = [ContainerPath '/' SessInfo.fstem_regT1_ref];
  if ~isempty(parms.infix), fstem_ref = [fstem_ref '_' parms.infix]; end;
end;

% check input, save input file names
data_fnames = [];
data_fstems = [];
reg_fnames = [];
motion_fnames = [];
for i=1:length(parms.snums)
  s = parms.snums(i);
  % check for bad scan num
  if ~ismember(s,SessInfo.snums_valid)
    fprintf('%s: ERROR: bad BOLD Scan Num (%d)\n',mfilename,s);
    errcode = 1;
    return;
  end;
  % check motion and data files exists
  fstem = sprintf('%s/%s',ContainerPath,ScanInfo(s).fstem);
  if ~isempty(parms.infix)
    fstem = [fstem '_' parms.infix];
  end;
  if parms.mc_flag
    motion_fnames{i} = [fstem '_motion.1D'];
    if ~exist(motion_fnames{i},'file')
      fprintf('%s: ERROR: motion 1D file %s not found\n',...
        mfilename,motion_fnames{i});
      errcode = 1;
      return;
    end;
  end;
  data_fstems{i} = fstem;
  data_fnames{i} = sprintf('%s.mgz',fstem);
  if ~exist(data_fnames{i},'file')
    fprintf('%s: ERROR: data file %s not found\n',mfilename,data_fnames{i});
    errcode =1 ;
    return;
  end;
  if parms.paint_flag || parms.resamp_flag
    % check register.dat exists
    if parms.mc_inter_flag
      reg_fnames{i} = sprintf('%s_register.dat',fstem_ref);
    else
      reg_fnames{i} = sprintf('%s_register.dat',fstem);
    end;
    if ~exist(reg_fnames{i},'file')
      fprintf('%s: WARNING: reg file %s not found (cannot resample to T1 or paint to surface)\n',...
        mfilename,reg_fnames{i});
      parms.resamp_flag = 0;
      parms.paint_flag = 0;
    end;
  end;
end;
SessInfo.nscans = length(parms.snums);
if ~SessInfo.nscans, return; end;


% check num_conds
if length(parms.num_conds) ~= SessInfo.nscans
  error('number of elements in num_conds (%d) does not match number of scans (%d)',...
          length(parms.num_conds),SessInfo.nscans);
end;
if isempty(stim_fnames)
  error('%s: stim_fnames is emtpy');
end;
if ~iscell(stim_fnames), stim_fnames = {stim_fnames}; end;
if numel(stim_fnames) ~= sum(parms.num_conds)
 error('number of elements in stim_fnames (%d) does not match sum of num_conds (%d)',...
   numel(stim_fnames),sum(parms.num_conds));
end;

% make absolute paths for stim_fnames if relative paths were given
for i=1:numel(stim_fnames)
  if isempty(regexp(stim_fnames{i},'^/')) % relative path
    stim_fnames{i} = sprintf('%s/%s',ContainerPath,stim_fnames{i});
  end;
end;

% reorganize into nested cell array
tmp_stim_fnames_nested = cell(SessInfo.nscans,1);
n_cond = 1;
for i=1:SessInfo.nscans
  tmp_stim_fnames = cell(parms.num_conds(i),1);
  for f=1:parms.num_conds(i)
    tmp_stim_fnames{f} = stim_fnames{n_cond};
    n_cond=n_cond+1;
  end;
  tmp_stim_fnames_nested{i} = tmp_stim_fnames;
end;
stim_fnames = tmp_stim_fnames_nested;


if SessInfo.nscans>1 && parms.concat_flag
  % construct cell array of sets of scan numbers
  snum_sets = cell(SessInfo.nscans+1,1);
  for i=1:SessInfo.nscans
    snum_sets{i} = parms.snums(i);
  end;
  snum_sets{SessInfo.nscans+1} = parms.snums;  
else
  % construct cell array of sets of scan numbers
  snum_sets = cell(SessInfo.nscans,1);
  for i=1:SessInfo.nscans
    snum_sets{i} = parms.snums(i);
  end;
end;

% create output dirs
outdirs = cell(size(snum_sets));
outstems = cell(size(snum_sets));
for i=1:length(snum_sets)
  s = snum_sets{i};
  if length(s)>1
    % construct scaninfix for combined output
    outstems{i} = sprintf('%d',parms.snums(1));
    for j=2:SessInfo.nscans
      outstems{i} = sprintf('%s_%d',outstems{i},parms.snums(j));
    end;
    outstems{i} = sprintf('%s_scans_%s',parms.fnamestem,outstems{i});
    if ~isempty(parms.infix)
      outstems{i} = sprintf('%s_%s',outstems{i},parms.infix);
    end;
    outdirs{i} = [ContainerPath '/' outstems{i} '_analysis'];
  else
    [tmp,outstems{i}] = fileparts(data_fstems{i});
    outdirs{i} = sprintf('%s/%s_scan_%d',ContainerPath,parms.fnamestem,s);
    if ~isempty(parms.infix)
      outdirs{i} = sprintf('%s_%s',outdirs{i},parms.infix);
    end;
    outdirs{i} = [outdirs{i} '_analysis'];
  end;
  mmil_mkdir(outdirs{i});
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deconvolution and GLM testing, resample to T1, paint for individual scans

% do analysis for each scan and optionally for all scans combined
for i=1:length(snum_sets)
  s = snum_sets{i};
  [tmp,s_ind] = intersect(parms.snums,s);
  if length(s)>1
    fname_in = data_fnames(s_ind);
    if parms.mc_flag
      fname_motion = motion_fnames(s_ind);
    else
      fname_motion = [];
    end;
    stim1Ds = stim_fnames(s_ind);
  else
    fname_in = data_fnames{s_ind};
    if parms.mc_flag
      fname_motion = motion_fnames{s_ind};
    else
      fname_motion = [];
    end;
    stim1Ds = stim_fnames{s_ind};
  end;
  outdir = outdirs{i};
  outstem = outstems{i};

  [fname_stats,fname_info,fnames_iresp] = mmil_3dDeconv(fname_in,stim1Ds,...
    'fname_motion',fname_motion,...
    'outdir',outdir,...
    'outstem',outstem,...
    'skipTRs',parms.skipTRs,...
    'out_ext',parms.out_ext,...
    'minlag',parms.minlag,...
    'maxlag',parms.maxlag,...
    'stim_labels',parms.stim_labels,...
    'contrasts_flag',parms.contrasts_flag,...
    'conds_contrast',parms.conds_contrast,...
    'iresp_flag',parms.iresp_flag,...
    'forceflag',parms.forceflag);

  if ~exist(fname_stats)
    error('3dDeconv output file %s not found',fname_stats);
  end;
  if ~exist(fname_info)
    error('3dDeconv info file %s not found',fname_info);
  end;

  % resample to T1 space
  if parms.resamp_flag
    fname_regdat = reg_fnames{s_ind(1)};
    fname_ref = sprintf('%s/mri/T1.mgz',parms.FSPath);
    if ~exist(fname_ref,'file')
      error('file %s not found',fname_ref);
    end;

    resamp_args = {...
      'fname_regdat',fname_regdat,...
      'fname_ref',fname_ref,...
      'fname_out',fname_out,...
      'forceflag',(parms.forceflag | parms.force_repaint_flag)};

    % first frame of raw BOLD
    fname_reg = data_fnames{s_ind(1)};
    [fpath,fstem,fext] = fileparts(fname_reg);
    if ~isempty(parms.out_ext)
      fext = parms.out_ext;
    end;
    fname_out = sprintf('%s/%s_resT1_f0%s',fpath,fstem,fext);
    [M_ref2reg,subj,inplane,slicethick] = mmil_resample_by_regdat(fname_reg,...
      resamp_args{:},'frames',1);

    % stats
    fname_reg = fname_stats;
    [fpath,fstem,fext] = fileparts(fname_reg);
    if ~isempty(parms.out_ext)
      fext = parms.out_ext;
    end;
    fname_out = sprintf('%s/%s_resT1%s',fpath,fstem,fext);
    [M_ref2reg,subj,inplane,slicethick] = mmil_resample_by_regdat(fname_reg,...
      resamp_args{:});
    fname_stats = fname_out;

    % iresp
    if parms.iresp_flag
      for c=1:length(fnames_iresp)
        fname_reg = fnames_iresp{c};
        [fpath,fstem,fext] = fileparts(fname_in);
        fname_out = sprintf('%s/%s_resT1%s',fpath,fstem,fext);
        [M_ref2reg,subj,inplane,slicethick] = mmil_resample_by_regdat(fname_reg,...
          resamp_args{:});
        fnames_iresp{c} = fname_out;
      end;
    end;
  end;

  % paint to surface
  if parms.paint_flag
    if parms.resamp_flag
      fname_regdat = [];
    else
      fname_regdat = reg_fnames{s_ind(1)};
    end;

    % paint GLM stats
    fname_in = fname_stats;
    fname_stats = fs_paint(parms.subj,fname_in,...
      'regfile',fname_regdat,paint_args{:});

    % paint impulse response functions
    if parms.iresp_flag
      [fpath,fstem,fext] = fileparts(fname_in);
      if ~isempty(parms.out_ext)
        fext = parms.out_ext;
      end;
      for c=1:length(fnames_iresp)
        fname_in = fnames_iresp{c};
        fname_out = fs_paint(parms.subj,fname_in,...
          'regfile',fname_regdat,paint_args{:});
      end;
    end;
  end;

  % threshold coefs with fstats
  if ~isempty(parms.pthresh) && parms.pthresh~=0 && parms.pthresh~=1
    clear stats_info;
    load(fname_info);
    if ~exist('stats_info','var')
      error('file %s did not contain ''stats_info'' as expected',fname_info);
    end;
    if ~iscell(fname_stats), fname_stats = {fname_stats}; end;
    for f=1:length(fname_stats)
      fname_in = fname_stats{f};
      [fpath,fstem,fext] = fileparts(fname_in);
      n = regexp(fstem,'(?<stem>.+)-(?<hemi>\wh)','names');
      if ~isempty(n)
        hemi = n.hemi;
        fstem = n.stem;
      else
        hemi = [];
      end;
      % a hack to avoid thresholding data on sphere
      n = regexp(fstem,'(?<stem>.+)-sphere','names');
      if ~isempty(n)
        fstem = n.stem; % apply threshold to stats in native space only
      end;
      if ~isempty(hemi)
        fname_in = sprintf('%s/%s-%s%s',...
          fpath,fstem,hemi,fext);
        fname_out = sprintf('%s/%s_pthresh%0.1e-%s%s',...
          fpath,fstem,parms.pthresh,hemi,fext);
      else
        fname_in = sprintf('%s/%s%s',...
          fpath,fstem,fext);
        fname_out = sprintf('%s/%s_pthresh%0.1e%s',...
          fpath,fstem,parms.pthresh,fext);
      end;
      if ~exist(fname_out,'file') || parms.forceflag
        [vol,M] = fs_load_mgh(fname_in,[],[],0,1);
        % check that length of stats_info matches number of frames
        nframes = size(vol,4);
        if nframes ~= length(stats_info)
          error('mismatch between number of frames in %s and stats_info in %s',...
            fname_in,fname_info);
        end;
        % loop over coef frames, thresholding by corresponding fstat frame
        types = {stats_info.type};
        ind = find(strcmp(types,'coef'));
        vol_out = zeros([size(vol,1),size(vol,2),size(vol,3),length(ind)]);
        for j=1:length(ind)
          k = ind(j);
          if k+1>nframes
            error('expected more frames in %s',fname_in);
          end;
          tmp_vol_coef = vol(:,:,:,k);
          % can make assumptions about where fstats are relative to coefs
          tmp_vol_fstat = vol(:,:,:,k+1);
          dofs = stats_info(k+1).dofs;
          if length(dofs)~=2
            error('degrees of freedom for frame %d of %s should have two elements',...
              k+1,fname_info);
          end;
          fthresh = finv(1-parms.pthresh,dofs(1),dofs(2));
          tmp_vol_coef(tmp_vol_fstat<fthresh)=0;
          vol_out(:,:,:,j) = tmp_vol_coef;
        end;
        % save thresholded coefs in multi-frame file
        fs_save_mgh(vol_out,fname_out,M);
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate csh scripts to view results
if parms.paint_flag
  hemi = parms.hemilist{1};
  for i=1:length(snum_sets)
    s = snum_sets{i};
    if length(s)>1
      fname_view = sprintf('%s/view_BOLD_GLM_results.csh',...
        ContainerPath);
    else
      fname_view = sprintf('%s/view_BOLD_GLM_results_scan%d.csh',...
        ContainerPath,s);
    end;
    indir = outdirs{i};
    fstem = [outstems{i} '_3dDeconv'];
    if parms.resamp_flag
      fstem = [fstem '_resT1'];
    end;
    if ~isempty(parms.pthresh) && parms.pthresh~=0 && parms.pthresh~=1
      if parms.smoothsteps
        fstem = sprintf('%s-sm%d_pthresh%0.1e',...
          fstem,parms.smoothsteps,parms.pthresh);
      else
        fstem = sprintf('%s_pthresh%0.1e',fstem,parms.pthresh);
      end;
    end;
    [fpath,indir] = fileparts(indir);
    fid = fopen(fname_view,'wt');
    if fid==-1
      error('failed to open %s for writing',fname_view);
    end;
    fprintf(fid,'#!/bin/csh -f\n');
    fprintf(fid,'# automatically generated script for viewing BOLD GLM results\n');
    fprintf(fid,'\n');
    fprintf(fid,'if ($#argv == 0) then\n');
    fprintf(fid,'  set hemi = %s\n',hemi);
    fprintf(fid,'else\n');
    fprintf(fid,'  set hemi = $argv[1]\n',hemi);
    fprintf(fid,'endif\n');
    fprintf(fid,'\n');
    fprintf(fid,'if ($#argv > 1) then\n');
    fprintf(fid,'  set save_flag = $argv[2]\n',hemi);
    fprintf(fid,'else\n');
    fprintf(fid,'  set save_flag = 0\n');
    fprintf(fid,'endif\n');
    fprintf(fid,'\n');
    fprintf(fid,'if ($#argv > 2) then\n');
    fprintf(fid,'  set view = $argv[3]\n');
    fprintf(fid,'else\n');
    fprintf(fid,'  set view = %s\n',parms.view);
    fprintf(fid,'endif\n');
    fprintf(fid,'\n');
    fprintf(fid,'setenv SUBJECTS_DIR %s\n',parms.subjdir);
    fprintf(fid,'set subj = %s\n',parms.subj);
    fprintf(fid,'set procdir = %s\n',ContainerPath);
    fprintf(fid,'set fstem = %s\n',fstem);
    fprintf(fid,'set indir = %s\n',indir);
    fprintf(fid,'\n');
    fprintf(fid,'set indir = $procdir/$indir\n');
    fprintf(fid,'\n');
    fprintf(fid,'set smooth = %d\n',parms.tksmooth);
    fprintf(fid,'set fthresh = %0.1f\n',parms.fthresh);
    fprintf(fid,'set fmid = %0.1f\n',parms.fmid);
    fprintf(fid,'set fslope = %0.1f\n',parms.fslope);
    fprintf(fid,'set colscale = %d\n',parms.colscale);
    fprintf(fid,'\n');
    fprintf(fid,'if $save_flag then\n');
    fprintf(fid,'  set savestr = "-savetiff -outstem $fstem -offscreen"\n');
    fprintf(fid,'else\n');
    fprintf(fid,'  set savestr = " "\n');
    fprintf(fid,'endif\n');
    fprintf(fid,'\n');
    fprintf(fid,'fs_surfmgh.csh $subj $fstem $hemi \\\n');
    fprintf(fid,'  -indir $indir \\\n');
    fprintf(fid,'  -smooth $smooth \\\n');
    fprintf(fid,'  -fthresh $fthresh \\\n');
    fprintf(fid,'  -fmid $fmid \\\n');
    fprintf(fid,'  -fslope $fslope \\\n');
    fprintf(fid,'  -colscale $colscale \\\n');
    fprintf(fid,'  -view $view $savestr\n');
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
  end;
end;
