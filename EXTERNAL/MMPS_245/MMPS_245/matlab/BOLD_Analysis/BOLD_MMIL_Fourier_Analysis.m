function errcode = BOLD_MMIL_Fourier_Analysis(ContainerPath,varargin)
%function errcode = BOLD_MMIL_Fourier_Analysis(ContainerPath,[options])
%
% Purpose: Fourier analysis for phase-encoded (e.g. retinotopy) fMRI data
%
% Usage:
%  BOLD_MMIL_Fourier_Analysis(ContainerPath,'key1', value1,...);
%  e.g. BOLD_MMIL_Fourier_Analysis(ContainerPath,...
%         'FSPath',FSPath,'paint_flag',1,...
%         'snums',[3,4,5,6],'revflags',[0,1,0,1],'skipTRs',2);
%
% Required Input:
%  ContainerPath: Full path of BOLDPROC directory containing BOLD scans
%
% Optional Paramters:
%  'FSPath': Full path of FreeSurfer directory containing cortical surface recon
%    (required if paint_flag > 0)
%    { default = [] }
%  'mc_flag': [0|1] whether within-scan motion correction was done
%    Will use motion.1D files as nuisance regressors
%    {default = 1}
%  'mc_inter_flag': [0|1] whether between-scan motion correction was done
%    Allows for a single reference scan to be used for registration to FS
%    {default = 1}
%  'snums': vector of scan numbers on which to run Fourier analysis
%    if empty or omitted, will use all BOLD scans in ContainerPath
%    { default = [] }
%  'revflags': vector of 0's or 1's to specify whether to reverse phase for each scan
%    { default = [] } (zero for all)
%  'infix': BOLD file name infix (e.g. '', 'corr', 'corr_resBOLD')
%    { default = [] }
%  'stimfreq': stimulus frequency (cycles per scan)
%    { default = 8 }
%  'skipTRs': number of TRs at beginning of each run to ignore
%    { default = 0 }
%  'phase_offset': fraction of cycles to subtract from phases to account for
%      hemodynamic delay
%    { default = 0.1 }
%  'phase_offset_postrev': fraction of cycles to subtract from phases
%      to account for stimulus delay
%    { default = 0 }
%  'fstats_type': [0|1|2] how output Fourier components should be scaled
%    0: raw, no scaling
%    1: scaled by sqrt(F-ratio)
%    2: scaled by significance values (-log10(p))
%    {default = 2}
%  'cxfstatsflag' [0|1] whether to calculate cross-scan complex f-stats
%    (in addition to amplitude f-stats)
%    { default = 0 }
%  'paint_flag': [0|1|2] whether to paint results to surface
%     if 0, only do Fourier analysis, no painting to surface
%     if 1, paint results to surface after averaging across scans
%       (requires motion correction and register.dat file for reference scan)
%     if 2, paint results to surface before averaging across scans
%       (requires register.dat file for each scan if not motion corrected)
%    { default = 0 }
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
%  'resamp_flag': [0|1] whether to resample fourier results to 1x1x1
%     before painting
%    { default = 0 }
%  'force_repaint_flag': [0|1] whether to repaint even if output files exist
%    (do not redo Fourier calculations)
%    { default = 0 }
%  'sphere_flag': [0|1] whether to sample to icosohedral sphere
%    {default = 0}
%  'smoothsteps': smoothing steps on surface after painting
%    {default = 0}
%  'sphsmoothsteps': smoothing steps on spherical surface
%    {default = 0}
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
%  'forceflag': [0|1] whether to run calculations even if output files exist
%    { default = 0 }
%
% Created:  02/06/07 by Don Hagler
% Last Mod: 08/06/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

%     if 2, resample fseries to T1 before averaging across scans
%           (won't work because of memory limitations)
%  'resamp_flag',0,[0,1,2],...

errcode = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'mc_flag',false,[false true],...
  'mc_inter_flag',true,[false true],...
  'FSPath',[],[],...
  'snums',[],[],...
  'snums_valid',[],[],...
  'infix',[],[],...
  'revflags',[],[],...
  'stimfreq',8,[1,100],...
  'skipTRs',0,[0,10],...
  'fstats_type',2,[0:2],...
  'phase_offset',0.1,[-1,1],...
  'phase_offset_postrev',0,[-1,1],...
  'paint_flag',0,[0,2],...
  'projdist',1,[0,10],...
  'projfrac',0.5,[0,2],...
  'projfrac_flag',false,[false true],...
  'projdist_avg',[],[],...
  'projfrac_avg',[],[],...
  'paint_surf','white',[],...
  'resamp_flag',false,[false true],...
  'force_repaint_flag',false,[false true],...
  'smoothsteps',0,[0,1000],...
  'sphsmoothsteps',0,[0,1000],...
  'sphere_flag',false,[false true],...
  'surfname','white',[],...
  'mask_midbrain_flag',false,[false true],...
... % parameters for viewing results
  'tksmooth',2,[0,1000],...
  'fthresh',0,[0,100],...
  'fmid',1.5,[0.1,100],...
  'fslope',1.0,[0.1,100],...
  'view','med',{'lat','med','ven','pos','dor'},...
  'datatype','polar',{'polar','eccen'},...
  'surf','inflated',[],...
...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'detrend',2,[0:10],...
  'out_ext','.mgh',{'.mgh','.mgz'},...
  'cxfstatsflag',false,[false true],...
  'forceflag',false,[false true],...
...
  'fnamestem','BOLD',[],...
...
  'paint_tags',{'regfile','sphere_flag','sphere_ico','projfrac_flag',...
                'projdist','projfrac','projdist_avg','projfrac_avg',...
                'smoothsteps','sphsmoothsteps','subjdir','surfname'...
                'mask_midbrain_flag','overwrite_flag'},[],...
});
parms.overwrite_flag = (parms.forceflag || parms.force_repaint_flag);
paint_args = mmil_parms2args(parms,parms.paint_tags);

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
    errcode = 1;
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
nscans = length(parms.snums);
if ~nscans, return; end;

% check that length of parms.revflags matches parms.snums
if isempty(parms.revflags), parms.revflags = zeros(nscans,1); end;
if length(parms.revflags) ~= nscans
  error('length of revflags (%d) does not match number of scans (%d)',...
    length(parms.revflags),nscans);
end;

% construct scaninfix for combined output
scaninfix = sprintf('%d',parms.snums(1));
if length(parms.snums)>1
  for i=2:length(parms.snums)
    scaninfix = sprintf('%s_%d',scaninfix,parms.snums(i));
  end;
  scaninfix = sprintf('%s_scans_%s',parms.fnamestem,scaninfix);
end;
if ~isempty(parms.infix), scaninfix = sprintf('%s_%s',scaninfix,parms.infix); end;
if nscans==1
  scaninfix  = [parms.fnamestem scaninfix];
end;

% create output dirs
for i=1:nscans
  s = parms.snums(i);
  outdirs{i} = sprintf('%s/%s%d',ContainerPath,parms.fnamestem,s);
  if ~isempty(parms.infix)
    outdirs{i} = sprintf('%s_%s',outdirs{i},parms.infix);
  end;
  outdirs{i} = [outdirs{i} '_analysis'];
  mmil_mkdir(outdirs{i});
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fourier transform, resample to T1, paint for individual scans

for i=1:nscans
  s = parms.snums(i);
  fname_in = data_fnames{i};
  if parms.mc_flag
    fname_motion = motion_fnames{i};
  else
    fname_motion = [];
  end;
  outdir = outdirs{i};
  tags = {'skipTRs','stimfreq','detrend',...
    'out_ext','fstats_type','forceflag'};
  args = mmil_parms2args(parms,tags);
  [fnames_fstats,fnames_fseries] = fs_fourier(fname_in,...
    'outdir',outdir,'save_fseries_flag',1,...
    'fname_motion',fname_motion,args{:});
    
  % resample to T1 space
  if parms.resamp_flag    
    fname_regdat = reg_fnames{i};
    fname_ref = sprintf('%s/mri/T1.mgz',parms.FSPath);
    if ~exist(fname_ref,'file')
      error('file %s not found',fname_ref);
    end;

    % first frame of raw BOLD
    fname_reg = data_fnames{i};
    [fpath,fstem,fext] = fileparts(fname_reg);
    if ~isempty(parms.out_ext)
      fext = parms.out_ext;
    end;
    fname_out = sprintf('%s/%s_resT1_f0%s',fpath,fstem,fext);
    [M_ref2reg,subj,inplane,slicethick] = mmil_resample_by_regdat(fname_reg,...
      'fname_regdat',fname_regdat,...
      'frames',1,...
      'fname_ref',fname_ref,...
      'fname_out',fname_out,...
      'forceflag',(parms.forceflag | parms.force_repaint_flag));

    % fstats
    for j=1:length(fnames_fstats)
      fname_reg = fnames_fstats{j};
      [fpath,fstem,fext] = fileparts(fname_reg);
      infix = fstem(end); % 'r' or 'i'
      fstem = fstem(1:end-2);
      fname_out = sprintf('%s/%s_resT1_%s%s',fpath,fstem,infix,parms.out_ext);
      [M_ref2reg,subj,inplane,slicethick] = mmil_resample_by_regdat(fname_reg,...
        'fname_regdat',fname_regdat,...
        'fname_ref',fname_ref,...
        'fname_out',fname_out,...
        'forceflag',(parms.forceflag | parms.force_repaint_flag));
      fnames_fstats{j} = fname_out;
    end;

    % fseries
    if parms.resamp_flag==2  % this will run out of memory
      for j=1:length(fnames_fseries)
        fname_reg = fnames_fseries{j};
        [fpath,fstem,fext] = fileparts(fname_reg);
        infix = fstem(end); % 'r' or 'i'
        fstem = fstem(1:end-2);
        fname_out = sprintf('%s/%s_resT1_%s%s',fpath,fstem,infix,parms.out_ext);
        [M_ref2reg,subj,inplane,slicethick] = mmil_resample_by_regdat(fname_reg,...
          'fname_regdat',fname_regdat,...
          'fname_ref',fname_ref,...
          'fname_out',fname_out,...
          'forceflag',(parms.forceflag | parms.force_repaint_flag));
        fnames_fseries{j} = fname_out;
      end;
    end;
  end;

  % paint to surface
  if parms.paint_flag
    if parms.resamp_flag
      fname_regdat = [];
    else
      fname_regdat = reg_fnames{i};
    end;
    if parms.paint_flag==2
      for j=1:length(fnames_fseries)
        fname_in = fnames_fseries{j};
        % paint fourier series (for calculating session average)
        fs_paint(parms.subj,fname_in,...
          'regfile',fname_regdat,paint_args{:});
      end;
    end;
    % paint fourier stats
    for j=1:length(fnames_fstats)
      fname_in = fnames_fstats{j};
      % paint fourier stats (for viewing single scan results)
      fs_paint(parms.subj,fname_in,...
        'regfile',fname_regdat,paint_args{:});
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross-scan average of fourier series

if parms.paint_flag~=2 && parms.mc_inter_flag==0
  fprintf('%s: WARNING: no motion correction, so combining fourier series across scans may be invalid unless paint_flag=2\n',...
    mfilename);
%  return;
end;

suffix_list = {'r','i'};
full_suffix_list = {};
for j=1:length(suffix_list)
  suffix = suffix_list{j};
  if parms.resamp_flag==2
    full_suffix_list{j} = ['resT1_' suffix];
  else
    full_suffix_list{j} = suffix;
  end;
end;

outdir = sprintf('%s/%s_analysis',ContainerPath,scaninfix);
mmil_mkdir(outdir);
calc_avg_flag = 0;
for j=1:length(full_suffix_list)
  full_suffix = full_suffix_list{j};
  if parms.paint_flag==2 % avg on surface
    for k=1:length(hemilist)
      hemi = hemilist{k};
      fname_out = sprintf('%s/%s_avg_fseries_%s-%s.mgh',...
        outdir,scaninfix,full_suffix,hemi);
      if ~exist(fname_out,'file'), calc_avg_flag = 1; end;
    end;
  else % avg in volume
    fname_out = sprintf('%s/%s_avg_fseries_%s%s',...
      outdir,scaninfix,full_suffix,parms.out_ext);
    if ~exist(fname_out,'file'), calc_avg_flag = 1; end;
  end;
end;
if parms.cxfstatsflag
  if parms.paint_flag==2 % avg on surface
    for k=1:length(hemilist)
      hemi = hemilist{k};
      fname_out = sprintf('%s/%s_cxfstats-%s.mgh',...
        outdir,scaninfix,hemi);
      if ~exist(fname_out,'file'), calc_avg_flag = 1; end;
    end;
  else % avg in volume
    fname_out = sprintf('%s/%s_cxfstats%s',...
      outdir,scaninfix,parms.out_ext);
    if ~exist(fname_out,'file'), calc_avg_flag = 1; end;
  end;
end;

if calc_avg_flag || parms.forceflag
  % average across scans and save output
  if parms.paint_flag==2 % avg on surface
    for k=1:length(hemilist)
      hemi = hemilist{k};
      % construct list of files to average
      fnamelist = [];
      for i=1:nscans
        s = parms.snums(i);
        for j=1:length(full_suffix_list) % real or imaginary
          full_suffix = full_suffix_list{j};
          [fpath,fstem] = fileparts(data_fstems{i});
          fnamelist{i,j} = sprintf('%s/%s_fseries_%s-%s.mgh',...
            outdirs{i},fstem,full_suffix,hemi);
        end;
      end;
      fprintf('%s: calculating cross-scan average for hemi %s...\n',...
        mfilename,hemi);
      results = fs_fourier_avg(fnamelist,...
        'revflags',parms.revflags,...
        'phase_offset',parms.phase_offset,...
        'cxfstatsflag',parms.cxfstatsflag,...
        'stimfreq',parms.stimfreq,...
        'phase_offset_postrev',parms.phase_offset_postrev...
      );
      for j=1:length(full_suffix_list)
        full_suffix = full_suffix_list{j};
        data = getfield(results,sprintf('mean_%s',suffix));
        fname_out = sprintf('%s/%s_avg_fseries_%s-%s.mgh',...
          outdir,scaninfix,full_suffix,hemi);
        fprintf('%s: saving average to %s...\n',...
          mfilename,fname_out);
        fs_save_mgh(data,fname_out,results.M);
      end;
      % save cross-scan complex fstats
      if parms.cxfstatsflag
        data = getfield(results,sprintf('Fstat',suffix));
        fname_out = sprintf('%s/%s_cxfstats-%s.mgh',...
          outdir,scaninfix,hemi);
        fprintf('%s: saving complex fstats to %s...\n',...
          mfilename,fname_out);
        fs_save_mgh(data,fname_out,results.M);
      end;
    end;
  else % avg in volume
    % construct list of files to average
    fnamelist = [];
    for i=1:nscans
      s = parms.snums(i);
      for j=1:length(full_suffix_list) % real or imaginary
        full_suffix = full_suffix_list{j};
        [fpath,fstem] = fileparts(data_fstems{i});
        fnamelist{i,j} = sprintf('%s/%s_fseries_%s%s',...
          outdirs{i},fstem,full_suffix,parms.out_ext);
      end;
    end;
    fprintf('%s: calculating cross-scan average...\n',mfilename);
    results = fs_fourier_avg(fnamelist,...
        'revflags',parms.revflags,...
        'phase_offset',parms.phase_offset,...
        'cxfstatsflag',parms.cxfstatsflag,...
        'stimfreq',parms.stimfreq,...
        'phase_offset_postrev',parms.phase_offset_postrev...
      );
    for j=1:length(full_suffix_list)
      suffix = suffix_list{j};
      full_suffix = full_suffix_list{j};
      data = getfield(results,sprintf('mean_%s',suffix));
      fname_out = sprintf('%s/%s_avg_fseries_%s%s',...
        outdir,scaninfix,full_suffix,parms.out_ext);
      fprintf('%s: saving average to %s...\n',...
        mfilename,fname_out);
      fs_save_mgh(data,fname_out,results.M);
    end;
    % save cross-scan complex fstats
    if parms.cxfstatsflag
      data = getfield(results,sprintf('Fstat',suffix));
      fname_out = sprintf('%s/%s_cxfstats%s',...
        outdir,scaninfix,parms.out_ext);
      fprintf('%s: saving complex fstats to %s...\n',...
        mfilename,fname_out);
      fs_save_mgh(data,fname_out,results.M);
    end;
  end;
  clear results;
end;

% calculate F-ratio and pvalue from average fseries
if parms.paint_flag==2 % avg on surface
  fstem = [];
  for k=1:length(hemilist)
    hemi = hemilist{k};
    fnames = [];
    for j=1:length(full_suffix_list)
      full_suffix = full_suffix_list{j};
      fnames{j} = sprintf('%s/%s_avg_fseries_%s-%s.mgh',...
        outdir,scaninfix,full_suffix,hemi);
    end;
    fnames_fstats = fs_fourier(fnames,...
      'outdir',outdir,...
      'input_fseries_flag',1,...
      'stimfreq',parms.stimfreq,...
      'out_ext',parms.out_ext,...
      'fstats_type',parms.fstats_type,...
      'forceflag',parms.forceflag);
    if isempty(fstem)
      [fpath,fstem,fext] = fileparts(fname_in);
      pat = sprintf('(?<stem>\\w+)_%s',component_list{c});
      n = regexp(fstem,pat,'names');
      fstem_results = n.stem;
    end;
  end;
else
  fnames = [];
  for j=1:length(full_suffix_list)
    full_suffix = full_suffix_list{j};
    fnames{j} = sprintf('%s/%s_avg_fseries_%s%s',...
      outdir,scaninfix,full_suffix,parms.out_ext);
  end;
  fnames_fstats = fs_fourier(fnames,...
    'outdir',outdir,...
    'input_fseries_flag',1,...
    'stimfreq',parms.stimfreq,...
    'out_ext',parms.out_ext,...
    'fstats_type',parms.fstats_type,...
    'forceflag',parms.forceflag);
  if parms.resamp_flag
    % resample avg fstats to T1  
    for j=1:length(fnames_fstats)
      fname_regdat = sprintf('%s_register.dat',fstem_ref);
      fname_reg = fnames_fstats{j};
      [fpath,fstem,fext] = fileparts(fname_reg);
      infix = fstem(end); % 'r' or 'i'
      fstem = fstem(1:end-2);
      fname_out = sprintf('%s/%s_resT1_%s%s',fpath,fstem,infix,parms.out_ext);
      [M_ref2reg,subj,inplane,slicethick] = mmil_resample_by_regdat(fname_reg,...
        'fname_regdat',fname_regdat,...
        'fname_ref',fname_ref,...
        'fname_out',fname_out,...
        'forceflag',parms.overwrite_flag);
      fnames_fstats{j} = fname_out;
    end;

    % resample cxfstats to T1
    if parms.cxfstatsflag
      fname_regdat = sprintf('%s_register.dat',fstem_ref);
      fname_in = sprintf('%s/%s_cxfstats%s',...
        outdir,scaninfix,parms.out_ext);
      fname_out = sprintf('%s/%s_cxfstats_resT1%s',...
        outdir,scaninfix,parms.out_ext);
      [M_ref2reg,subj,inplane,slicethick] = mmil_resample_by_regdat(fname_in,...
        'fname_regdat',fname_regdat,...
        'fname_ref',fname_ref,...
        'fname_out',fname_out,...
        'forceflag',parms.overwrite_flag);
    end;
  end;
  if parms.paint_flag
    % paint average fourier stats
    fstem = [];
    if parms.resamp_flag
      fname_regdat = [];
    else
      fname_regdat = sprintf('%s_register.dat',fstem_ref);
    end;
    for j=1:length(fnames_fstats)
      fname_in = fnames_fstats{j};
      fs_paint(parms.subj,fname_in,...
        'regfile',fname_regdat,paint_args{:});
      if isempty(fstem)
        [fpath,fstem,fext] = fileparts(fname_in);
        pat = sprintf('(?<stem>\\w+)_%s',suffix_list{j});
        n = regexp(fstem,pat,'names');
        fstem_results = n.stem;
      end;
    end;

    % paint cross-scan complex fstats
    if parms.cxfstatsflag
      if parms.resamp_flag
        suffix = 'resT1';
      else
        suffix = [];
      end;
      if parms.resamp_flag
        fname_in = sprintf('%s/%s_cxfstats_%s%s',...
          outdir,scaninfix,suffix,parms.out_ext);
      else
        fname_in = sprintf('%s/%s_cxfstats%s',...
          outdir,scaninfix,parms.out_ext);
      end;
      fs_paint(parms.subj,fname_in,...
        'regfile',fname_regdat,paint_args{:});
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate csh scripts to view results
if ~isempty(fstem_results)
  hemi = parms.hemilist{1};
  for i=0:nscans
    if i==0
      indir = outdir;
      fstem = fstem_results;
      fname_view = sprintf('%s/view_BOLD_Fourier_%s_results.csh',...
        ContainerPath,parms.datatype);
      revflag = 0;
    else
      s = parms.snums(i);
      indir = outdirs{i};
      fstem = data_fstems{i};
      [fpath,fstem] = fileparts(data_fstems{i});
      fstem = [fstem '_fstats_pval'];
      if parms.resamp_flag
        fstem = [fstem '_resT1'];
      end;
      fname_view = sprintf('%s/view_BOLD_Fourier_%s_results_scan%d.csh',...
        ContainerPath,parms.datatype,s);
      revflag = parms.revflags(i);
    end;
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
    fprintf(fid,'set realfile = $fstem''_r''\n');
    fprintf(fid,'set imagfile = $fstem''_i''\n');
    fprintf(fid,'\n');
    fprintf(fid,'set smooth = %d\n',parms.tksmooth);
    fprintf(fid,'set fthresh = %0.1f\n',parms.fthresh);
    fprintf(fid,'set fmid = %0.1f\n',parms.fmid);
    fprintf(fid,'set fslope = %0.1f\n',parms.fslope);
    fprintf(fid,'set surf = %s\n',parms.surf);
    fprintf(fid,'\n');
    fprintf(fid,'if $save_flag then\n');
    fprintf(fid,'  set savestr = "-savetiff -outstem $fstem"# -offscreen"\n');
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
    fprintf(fid,'  -indir $indir \\\n');
    fprintf(fid,'  -smooth $smooth \\\n');
    fprintf(fid,'  -fthresh $fthresh \\\n');
    fprintf(fid,'  -fmid $fmid \\\n');
    fprintf(fid,'  -fslope $fslope \\\n');
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

