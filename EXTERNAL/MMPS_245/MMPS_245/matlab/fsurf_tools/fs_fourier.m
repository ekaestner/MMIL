function [fnames_fstats,fnames_fseries] = fs_fourier(fname,varargin)
%function [fnames_fstats,fnames_fseries] = fs_fourier(fname,varargin)
%
% Purpose: run Fourier analysis on input mgh volume
%   and output Fourier stats and/or Fourier series
%
% Usage:
%  fs_fourier(fname,'key1', value1,...); 
%
% Required Parameters:
%   fname: full or relative path name of input mgh file(s)
%     If input_fseries_flag, fname should be a cell array containing
%     the file names of the real and imaginary components
%     (e.g. {'BOLD_avg_fseries_r.mgh','BOLD_avg_fseries_i.mgh'})
%
% Optional parameters:
%  'outdir': output directory
%    If empty, will write output files to directory containing input fname
%    {default: []}
%  'input_fseries_flag': [0|1] whether to treat input mgh file as fseries
%    and only calculate fstats
%    {default: 0}
%  'save_fseries_flag': [0|1] whether to save raw Fourier components for
%    all frequencies (i.e. Fourier series)
%    {default: 0}
%  'save_fstats_flag': [0|1] whether to save Fourier components for stimulus
%    frequency
%    {default: 1}
%  'fstats_type': [0|1|2] how output Fourier components should be scaled
%    0: raw, no scaling
%    1: scaled by sqrt(F-ratio)
%    2: scaled by significance values (-log10(p))
%    {default: 2}
%  'stimfreq': stimulus frequency (cycles per scan)
%    {default: 8}
%  'skipTRs': number of initial repetitions to omit from analysis
%    {default: 0}
%  'freq_low': low cutoff frequency (ignore frequencies lower than this)
%    {default: 3}
%  'freq_high': high cutoff frequency (ignore frequencies higher than this)
%    {default: Inf}
%  'freq_omit': ignore this frequency for fstat calculation
%    Note: 1st and 2nd harmonics of stimfreq are automatically excluded
%    {default: []}
%  'norm_flag': [0|1] whether to normalize input timeseries
%    by mean for each voxel (new mean = 100) before doing Fourier calculations
%    {default: 1}
%  'thresh': when normalize to mean, set voxels with original values
%    less than this to 0
%    {default: 10}
%  'detrend': [0|1|2] whether and how to detrend input timeseries
%    0: no detrend (not recommended)
%    1: linear detrend
%    2: quadratic detrend
%    3: cubic detrend
%    {default: 2}
%  'fname_motion': text file containing 6 columns of motion estimates
%     (e.g. from AFNI's 3dvolreg)
%    If supplied, will regress out motion from timeseries data
%    {default: []}
%  'out_ext': output file extension ('.mgh' or '.mgz')
%    if empty, will use input file extension
%    {default: []}
%  'forceflag': [0|1] toggle overwrite existing output files
%    {default: 0}
%
% Output:
%   fnames_fstats: cell array of output file names (real and imaginary)
%     for Fourier stats
%   fnames_fseries: cell array of output file names (real and imaginary)
%     for Fourier series
%
%   Output files will be created with names based on input fname
%     dictated by save_fstats_flag, save_fseries_flag, etc.
%
% Created:  07/21/08 Don Hagler
% Last Mod: 10/28/14 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

component_list = {'r','i'};
fnames_fstats = [];
fnames_fseries = [];

if (~mmil_check_nargs(nargin, 2)), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir',[],[],...
  'input_fseries_flag',false,[false true],...
  'save_fstats_flag',true,[false true],...
  'save_fseries_flag',false,[false true],...
  'fstats_type',2,[0:2],...
  'stimfreq',8,[1 Inf],...
  'skipTRs',0,[0 Inf],...
  'freq_low',3,[1 Inf],...
  'freq_high',Inf,[1 Inf],...
  'freq_omit',[],[1 Inf],...
  'norm_flag',true,[false true],...
  'thresh',10,[0 Inf],...
  'detrend',2,[0:2],...
  'fname_motion',[],[],...
  'out_ext',[],{'.mgh','.mgz'},...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct output file names, check if they exist

hemi = [];
resamp_flag = 0;
if parms.input_fseries_flag
  if ~iscell(fname) || length(fname)~=2
    error('input fseries fname must be a cell array of file names for real and imaginary components');
  end;
  fstem = [];
  for c=1:length(component_list)
    fname_tmp = fname{c};
    if ~exist(fname_tmp,'file')
      error('file %s not found',fname_tmp);
    end;
    if isempty(fstem) % assume r and i have stem fstem and fext
      [fpath,fstem,fext] = fileparts(fname_tmp);
      if isempty(fpath), fpath = pwd; end;
      if ~ismember(fext,{'.mgh','.mgz'})
        error('input file must be mgh or mgz file type (has %s extension)',...
          fext);
      end;

      % what if input fname has resT1?
      pat = sprintf('(?<stem>\\w+)_fseries_resT1_%s',component_list{c});
      n = regexp(fstem,pat,'names');
      if ~isempty(n)
        resamp_flag = 1;
      else
        pat = sprintf('(?<stem>\\w+)_fseries_%s',component_list{c});
        n = regexp(fstem,pat,'names');
        if isempty(n)
          error('file stem %s does not conform to expected pattern (stem_fseries_%s)',...
            fstem,component_list{c});
        end;
      end;
      fstem = n.stem;

      % what if input fname has -lh or -rh infix?
      pat = '(?<hemi>-\wh$)';
      n = regexp(fstem,pat,'names');
      if ~isempty(n)
        hemi = n.hemi;
      end;
    end;
  end;
else
  if ~exist(fname,'file')
    error('file %s not found',fname);
  end;
  [fpath,fstem,fext] = fileparts(fname);
  if isempty(fpath), fpath = pwd; end;
  if ~ismember(fext,{'.mgh','.mgz'})
    error('input file must be mgh or mgz file type (has %s extension)',...
      fext);
  end;
end;
if isempty(parms.out_ext), parms.out_ext = fext; end;
if ~isempty(hemi), parms.out_ext = [hemi parms.out_ext]; end;
if isempty(parms.outdir), parms.outdir = fpath; end;

if parms.input_fseries_flag && parms.save_fseries_flag
  error('input_fseries_flag and save_fseries_flag cannot both be true');
end;
if ~parms.save_fstats_flag && ~parms.save_fseries_flag
  error('save_fstats_flag and save_fseries_flag cannot both be false');
end;

run_flag = 0;
if parms.save_fstats_flag
  for c=1:length(component_list)
    infix = [];
    switch parms.fstats_type
      case 0
        infix = 'raw';
      case 1
        infix = 'ratio';
      case 2
        infix = 'pval';
    end;
    if resamp_flag
      fnames_fstats{c} = sprintf('%s/%s_fstats_%s_resT1_%s%s',...
        parms.outdir,fstem,infix,component_list{c},parms.out_ext);
    else
      fnames_fstats{c} = sprintf('%s/%s_fstats_%s_%s%s',...
        parms.outdir,fstem,infix,component_list{c},parms.out_ext);
    end;
    if ~exist(fnames_fstats{c},'file')
      run_flag = 1;
    end;
  end;
else
  fnames_fstats = [];
end;
if parms.save_fseries_flag
  for c=1:length(component_list)
    fnames_fseries{c} = sprintf('%s/%s_fseries_%s%s',...
      parms.outdir,fstem,component_list{c},parms.out_ext);
    if ~exist(fnames_fseries{c},'file')
      run_flag = 1;
    end;
  end;
else
  fnames_fseries = [];
end;
if ~run_flag && ~parms.forceflag
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~parms.input_fseries_flag
  % load input mgh volume
  [vol,M] = fs_load_mgh(fname,[],[],0,1);
  [nx,ny,nz,nframes] = size(vol);

  % check that nframes is long enough
  nfreqs = round((nframes-parms.skipTRs)/2);
  if parms.stimfreq>nfreqs-1
    error('stimfreq must be < number of (numTRs/2)-1 (%d)',nfreqs-1);
  end;

  % normalize to mean, remove trend, remove motion
  tags = {'skipTRs','norm_flag','detrend','thresh','fname_motion'};
  args = mmil_parms2args(parms,tags);
  vol = mmil_detrend_vol(vol,args{:});

  % calculate Fourier series
  vol_fseries = fft(vol,[],4);
  % discard the second half of the frequencies
  vol_fseries = vol_fseries(:,:,:,1:nfreqs);

  % save output
  if parms.save_fseries_flag
    for c=1:length(component_list)
      fname_out = fnames_fseries{c};
      if c==1
        vol_tmp = real(vol_fseries);
      else
        vol_tmp = imag(vol_fseries);
      end;
      fprintf('%s: saving Fourier series to %s...\n',mfilename,fname_out);
      fs_save_mgh(vol_tmp,fname_out,M);
    end;
  end;
else
  % load input mgh volumes
  [vol_r,M,mrparms,volsz] = fs_load_mgh(fname{1},[],[],0,1);
  vol_i = fs_load_mgh(fname{2},[],[],0,1);
  vol_fseries = complex(vol_r,vol_i);
  clear vol_r vol_i;
  nfreqs = volsz(4);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate Fourier stats

if parms.save_fstats_flag
  vol_fstats = squeeze(vol_fseries(:,:,:,parms.stimfreq+1));
  % plus 1 because of matlab 1-based indexing and first value is 0 cycles/scan
  if parms.fstats_type>0
    freq_low = round(min(max(parms.freq_low+1,1),nfreqs));
    freq_high = round(min(max(parms.freq_high+1,1),nfreqs));
    noise_freqs = [freq_low:freq_high];
    harmonics = [parms.stimfreq,2*parms.stimfreq,3*parms.stimfreq]+1;
    harmonics = [harmonics-1,harmonics,harmonics+1];
    noise_freqs = setdiff(noise_freqs,[harmonics,parms.freq_omit+1]);
    dof_stim = 2;
    dof_noise = 2*length(noise_freqs);
    vol_stim = abs(vol_fseries(:,:,:,parms.stimfreq+1)).^2;
    vol_noise = sum(abs(vol_fseries(:,:,:,noise_freqs)).^2,4);
    vol_F = (vol_stim/dof_stim)./(eps + vol_noise/dof_noise);
    vol_fstats_ang = angle(vol_fstats);
    switch parms.fstats_type
      case 1 % ratio
        vol_fstats_amp = sqrt(vol_F);
      case 2 % pval
        vol_fstats_amp = -log10(1-fcdf(vol_F,dof_stim,dof_noise));
        vol_fstats_amp(vol_fstats_amp>100) = 100;
    end;
    vol_fstats_r = vol_fstats_amp.*cos(vol_fstats_ang);
    vol_fstats_i = vol_fstats_amp.*sin(vol_fstats_ang);
    vol_fstats = complex(vol_fstats_r,vol_fstats_i);
  end;
  % save output
  for c=1:length(component_list)
    fname_out = fnames_fstats{c};
    if c==1
      vol_tmp = real(vol_fstats);
    else
      vol_tmp = imag(vol_fstats);
    end;
    fs_save_mgh(vol_tmp,fname_out,M);
  end;
end;
