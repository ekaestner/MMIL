function results = fs_fourier_avg(fnamelist,varargin)
%function results = fs_fourier_avg(fnamelist,[options])
%
% Usage:
%  results = fs_fourier_avg(fnamelist,'key1', value1,...);
%
% Required Input:
%   fnamelist: cell matrix of input full path file names
%     Each row should correspond to a single scan / subject
%     Must have two columns corresponding to real and imaginary components
%
%     All input files should be mgh format surface or volume files
%       for a single subject or registered to common space
%       (e.g. spherical ico surface)
%
% Optional Input:
%   'frames': vector of frame numbers (e.g. multiple frequencies)
%     If empty, will loop over all frames
%     {default = []}
%   'revflags': vector of 0's and 1's for each subject (each row of fnamelist)
%     determining whether to set the imaginary components negative
%     (e.g. reverse phase)
%     {default = [0 0 0....]}
%   'phase_offset': phase angle (in cycles) subtracted from
%     complex valued input before applying revflags
%     (e.g. to correct for hemodynamic delay)
%     may be vector with a value for each row of fnamelist
%     {default = 0}
%   'phase_offset_postrev': phase angle (in cycles) subtracted from
%       complex valued input after applying revflags
%       (e.g. to correct for stimulus delay)
%     may be vector with a value for each row of fnamelist
%     {default = 0}
%   'r_max_vec': maximum eccentricity for each scan
%     (i.e. for eccentricity mapping)
%     if not empty, will be used to adjust phases
%     should be a vector of values for each row of fnamelist
%     {default = []}
%   'r_max': maximum eccentricity for each scan
%     if r_max_vec is empty, will use r_max for all
%     {default = 12.5}
%   'r_min_vec': minimum eccentricity for each scan
%     should be a vector of values for each row of fnamelist
%     if empty, will assume to be r_min_factor * r_max_vec
%     {default = []}
%   'r_min_factor': value multiplied by r_max_vec to get r_min_vec (if empty)
%     {default = 0.02}
%   'logtrans_flags': vector of 0's and 1's for each row of fnamelist
%     indicating whether log transform was used for eccentricity mapping
%     if empty, will assume logtrans_flag = 0 for all scans
%     {default = []}
%   'cxfstatsflag': [0|1] whether to calculate complex F-stats
%     {default = 0}
%   'stimfreq': stimulus frequency for cxfstats (cycles per scan)
%     If 0, calculate complex F-stats for all frequencies
%     {default: 8}
%   'verbose_flag': [0|1] whether to print status messages
%     {default = 1}
%
% Output:
%   results: structure containing cross-scan average
%
% Created:  08/18/08 by Don Hagler
% Last Mod: 03/11/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'frames',[],[],...
  'phase_offset',0,[-1,1],...
  'phase_offset_postrev',0,[-1,1],...
  'revflags',[],[],...
  'r_max_vec',[],[],...
  'r_max',12.5,[],...
  'r_min_vec',[],[],...
  'r_min_factor',0.02,[],...
  'logtrans_flags',[],[],...
  'cxfstatsflag',false,[false true],...
  'stimfreq',8,[1 Inf],...
  'verbose_flag',true,[false true],...
});
results = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input

if isempty(fnamelist)
  error('input file name list is empty');
elseif ~iscell(fnamelist)
  fnamelist = {fnamelist};
end;
[nscans,nconds] = size(fnamelist);
if nconds~=2
  error('input file name list must have two columns corresponding to real and imaginary components');
end;

if ~isempty(parms.frames) && min(parms.frames)<1
  error('frame numbers must be > 1');
end;

if parms.verbose_flag
  tic
  fprintf('%s: checking input files...\n',mfilename);
end;
% check files
results.orig_volsz = [];
results.M = [];
results.frames = [];
for s=1:nscans
  for c=1:nconds
    if isempty(fnamelist{s,c})
      error('fnamelist{%d,%d} is empty',s,c);
    end;
    [vol,M,mrparms,volsz] = fs_load_mgh(fnamelist{s,c},[],[],1); % header-only
    if isempty(results.frames)
      if isempty(parms.frames)
        results.frames = [1:volsz(4)];
      else
        results.frames = parms.frames;
      end;
      results.nframes = length(results.frames);
      results.nfreqs = round(results.nframes/2);
    end;
    if (results.nframes>1) & (parms.stimfreq>results.nfreqs-1)
      error('stimfreq must be < (nframes/2)-1 (%d)',results.nfreqs-1);
    end;
    if max(results.frames)>volsz(4)
      error('only %d frames in file %s',volsz(4),fnamelist{s,c});
    end;
    if isempty(results.orig_volsz)
      results.orig_volsz = volsz;
      results.M = M;
    end;
    if any(results.orig_volsz~=volsz)
      error('size of input volume for\n%s (%d,%d,%d,%d)\ndoes not match that of\n%s (%d,%d,%d,%d)',...
        fnamelist{s,c},volsz(1),volsz(2),volsz(3),volsz(4),...
        fnamelist{1,1},results.orig_volsz(1),results.orig_volsz(2),results.orig_volsz(3),results.orig_volsz(4));
    end;
  end;
end;
if parms.verbose_flag, toc; end;

results.nvals = prod(results.orig_volsz(1:3));
results.volsz = [results.orig_volsz(1:3) results.nframes];

if isempty(parms.revflags)
  parms.revflags = zeros(nscans,1);
elseif length(parms.revflags)~=nscans
  error('number of revflags must match number of rows of fnamelist');
end;

if length(parms.phase_offset)==1
  parms.phase_offset = ones(nscans,1)*parms.phase_offset;
end;
if length(parms.phase_offset)~=nscans
  error('number of phase_offset values must match number of rows of fnamelist');
end;

if length(parms.phase_offset_postrev)==1
  parms.phase_offset_postrev = ones(nscans,1)*parms.phase_offset_postrev;
end;
if length(parms.phase_offset_postrev)~=nscans
  error('number of phase_offset_postrev values must match number of rows of fnamelist');
end;

if (~isempty(parms.r_max_vec) && length(unique(parms.r_max_vec))>1) ||...
   (~isempty(parms.r_min_vec) && length(unique(parms.r_min_vec))>1) ||...
   (~isempty(parms.logtrans_flags) && length(unique(parms.logtrans_flags))>1)
  parms.adjust_r_flag = 1;   
  if isempty(parms.r_max_vec)
    parms.r_max_vec = ones(nscans,1)*parms.r_max;
  elseif length(parms.r_max_vec)~=nscans
    error('number of elements in r_max_vec must match number of rows of fnamelist');
  end;
  if isempty(parms.r_min_vec)
    parms.r_min_vec = parms.r_min_factor*parms.r_max_vec;
  elseif length(parms.r_min_vec)~=nscans
    error('number of r_min_vec must match number of rows of fnamelist');
  end;
  if isempty(parms.logtrans_flags)
    parms.logtrans_flags = zeros(nscans,1);
  elseif length(parms.logtrans_flags)~=nscans
    error('number of logtrans_flags must match number of rows of fnamelist');
  end;    
else
  parms.adjust_r_flag = 0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data and calculate sums
results.n = 0;
results.mean_r = zeros(results.nvals,results.nframes,'single');
results.mean_i = zeros(results.nvals,results.nframes,'single');
if parms.cxfstatsflag
  results.var_r = zeros(results.nvals,results.nframes,'single');
  results.var_i = zeros(results.nvals,results.nframes,'single');
end;
for f=1:results.nframes
  frame = results.frames(f);
  if parms.verbose_flag
    tic
    fprintf('%s: summing data for frame %d of %d...\n',mfilename,frame,results.orig_volsz(4));
  end;
  if parms.adjust_r_flag
    final_r_max = max(parms.r_max_vec);
    final_r_min = min(parms.r_min_vec);
    final_logtrans_flag = min(parms.logtrans_flags);
  end;
  for s=1:nscans
    phase_offset = parms.phase_offset(s);
    phase_offset_postrev = parms.phase_offset_postrev(s);  
  
    vec_r = fs_load_mgh(fnamelist{s,1},[],frame);
    vec_r = reshape(vec_r,[prod(size(vec_r)) 1]);
    vec_i = fs_load_mgh(fnamelist{s,2},[],frame);
    vec_i = reshape(vec_i,[prod(size(vec_i)) 1]);
    if phase_offset~=0
      vec_ampl = hypot(vec_r,vec_i);
      vec_phase = atan2(vec_i,vec_r) - phase_offset*2.0*pi;
      vec_r = vec_ampl.*cos(vec_phase);
      vec_i = vec_ampl.*sin(vec_phase);
    end;
    if parms.revflags(s)
      vec_i = -vec_i;
    end;
    if phase_offset_postrev~=0
      vec_ampl = hypot(vec_r,vec_i);
      vec_phase = atan2(vec_i,vec_r) - phase_offset_postrev*2.0*pi;
      vec_r = vec_ampl.*cos(vec_phase);
      vec_i = vec_ampl.*sin(vec_phase);
    end;

    % adjust for differences in r_max, r_min, or logtrans_flags
    if parms.adjust_r_flag
      r_min = parms.r_min_vec(s);
      r_max = parms.r_max_vec(s);
      logtrans_flag = parms.logtrans_flags(s);

      % calculate ampl and phase
      vec_ampl = hypot(vec_r,vec_i);
      vec_phase = atan2(vec_i,vec_r); % in radians
    
      % calculate eccentricity (degrees visual angle) from phase
      vec_phase = vec_phase/(2*pi); % in cycles
      if logtrans_flag
        vec_visang = r_min*exp(vec_phase*log(r_max/r_min));
      else
        vec_visang = r_min + (r_max - r_min)*vec_phase;
      end;
      
      % calculate phase from eccentricity using final_r_max, etc.
      if final_logtrans_flag
        vec_phase = log(vec_visang/(final_r_min*log(final_r_max/final_r_min)));
      else
        vec_phase = (vec_visang - final_r_min)/(final_r_max - final_r_min);
      end;
      vec_phase = 2*pi*vec_phase;
      
      % recalculate real and imag
      vec_r = vec_ampl.*cos(vec_phase);
      vec_i = vec_ampl.*sin(vec_phase);
    end;
    if f==1
      results.n = results.n + 1;
    end;
    results.mean_r(:,frame) = results.mean_r(:,frame) + vec_r;
    results.mean_i(:,frame) = results.mean_i(:,frame) + vec_i;
    if parms.cxfstatsflag
      results.var_r(:,frame) = results.var_r(:,frame) + vec_r.^2;
      results.var_i(:,frame) = results.var_i(:,frame) + vec_i.^2;
    end;
  end;
  if parms.verbose_flag, toc; end;
end;

% calculate mean
if parms.verbose_flag
  tic
  fprintf('%s: calculating means...\n',mfilename);
end;
results.mean_r = reshape(results.mean_r,results.volsz);
results.mean_r = results.mean_r / (eps+results.n);
results.mean_i = reshape(results.mean_i,results.volsz);
results.mean_i = results.mean_i / (eps+results.n);

% reshape var back to volume
if parms.cxfstatsflag
  results.var_r = reshape(results.var_r,results.volsz);
  results.var_i = reshape(results.var_i,results.volsz);
end;

if parms.verbose_flag, toc; end;

if parms.cxfstatsflag
  if parms.verbose_flag
    tic
    fprintf('%s: calculating stats...\n',mfilename);
  end;

  if parms.stimfreq>0 && results.nframes>1
    % calculate variance
    tmp_var_r = results.var_r(:,:,:,parms.stimfreq+1);
    tmp_var_i = results.var_i(:,:,:,parms.stimfreq+1);
    tmp_mean_r = results.mean_r(:,:,:,parms.stimfreq+1);
    tmp_mean_i = results.mean_i(:,:,:,parms.stimfreq+1);
    tmp_volsz = results.volsz(1:3);

    tmp_var_r = (results.n*tmp_var_r - tmp_mean_r.^2)./ ...
                    (eps+results.n*(results.n-1));
    tmp_var_i = (results.n*tmp_var_i - tmp_mean_i.^2)./ ...
                    (eps+results.n*(results.n-1));

    % calculate complex f-stats
    % F ratio is (X2/dof1)/(X2/dof2) with X2 a chi-squared stat
    %   sum of squares is a chi-squared stat
    %   dof is the number of elements in the sum
    %     dof1 = 2 for real and imag
    %     dof2 = 2*N-2 or 2*N (depending on whether zero mean)
    % F ratio = (numer/dof1)/(denom/dof2)
    %   numer equals sum of squared averages
    %   denom equals sum of variance
    results.dof1 = 2;
    results.dof2 = 2*results.n-2;
    numer = (tmp_mean_r.^2 + tmp_mean_i.^2)/results.dof1;
    denom = (tmp_var_r + tmp_var_i)/results.dof2;
    results.Fstat = numer./denom;
  else % do for all frames -- could use lots of memory if many voxels
    % calculate variance
    results.var_r = (results.n*results.var_r - results.mean_r.^2)./ ...
                    (eps+results.n*(results.n-1));
    results.var_i = (results.n*results.var_i - results.mean_i.^2)./ ...
                    (eps+results.n*(results.n-1));
    % calculate complex f-stats
    results.dof1 = 2;
    results.dof2 = 2*results.n-2;
    numer = (results.mean_r.^2 + results.mean_i.^2)/results.dof1;
    denom = (results.var_r + results.var_i)/results.dof2;
    results.Fstat = numer./denom;
  end;

  % calcuate pval
  results.pval = 1-fcdf(results.Fstat,results.dof1,results.dof2);
  results.pval(results.pval<=0) = 1e-40;
  results.pval = -log10(results.pval);

  % calcuate stdev
  results.stdv = sqrt((results.var_r + results.var_i)/2);

  if parms.verbose_flag, toc; end;
end;


% null hypothesis for complex f-stats is zero mean, so could calculate 
%   variance assuming mean is zero with var = sumsq/N
%   and dof2 = 2*N

% another option could be to scale amplitudes to 1
%   so f-stats only depend on phase

