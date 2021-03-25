function rc_mask_pol(subj,varargin)
%function rc_mask_pol(subj,varargin)
%
% Usage:
%  rc_mask_pol(subj,'key1', value1,...);
%
% Required parameters:
%   subj:  a string specifying the subject name (freesurfer recon dir)
%
% Optional parameters:
%  'polstem' - file stem for polar angle retinotopy data
%     e.g. for files polret_r-lh.mgh and polret_i-lh.mgh, polstem = "polret"
%     {default = 'polret'}
%  'poldir' - input directory for polar angle data
%     {default = pwd}
%  'eccstem' - file stem for eccentricity retinotopy data
%     e.g. for files eccret_r-lh.mgh and eccret_i-lh.mgh, eccstem = "eccret"
%     {default = 'eccret'}
%  'eccdir' - input directory for eccentricity data
%     {default = pwd}
%  'maskstem' - file stem for mask data
%     e.g. block-design / event related stats
%     if supplied, phase ecc ('eccstem') data will be ignored
%     {default = []}
%  'maskdir' - input directory for mask data
%     {default = pwd}
%  'outstem'  - output file stem (omit extension, hemi)
%    {default = 'polret-masked'}
%  'outdir' - output directory for fieldsign calculations
%     {default = pwd}
%  'hemilist' - cell array containing list of hemispheres to process
%    {default = {'lh','rh'}}
%  'presmooth' - smoothing steps on surface applied to mask before threshold
%    {default = 0}
%  'thresh' - threshold applied to mask data
%     {default = 1}
%  'postsmooth' - smoothing steps on surface applied to mask after threshold
%    {default = 0}
%  'r0' - starting radius (degrees visual angle)
%     when masking with phase encoded eccentricity data
%  'r1' - ending radius (degrees visual angle)
%     when masking with phase encoded eccentricity data
%  'rmin' - minimum radius (degrees visual angle)
%    for phase encoded eccentricity stimulus
%    determines how phase corresponds to degrees visual angle
%    {default = 0.25}
%  'rmax' - maximum radius (degrees visual angle)
%    for phase encoded eccentricity stimulus
%    {default = 12.5}
%  'logtrans_flag' - [0|1] whether log transform was used with
%    phase encoded eccentricity stimulus
%    {default = 1}
%  'mask_frame' - frame number from mask file
%    ignored if maskstem is empty
%    {default = 1}
%  'surf' - surface to load
%    {default = 'white'}
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'forceflag' - [0|1] whether to overwrite existing output files
%    {default = 0}
%
% Created:  10/09/08 by Don Hagler
% Last Mod: 03/08/11 by Don Hagler
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms( varargin, {...
  'polstem','polret',[],...
  'maskstem',[],[],...
  'eccstem','eccret',[],...
  'outstem','polret-masked',[],...
  'poldir',pwd,[],...
  'maskdir',pwd,[],...
  'eccdir',pwd,[],...
  'outdir',pwd,[],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'intype','mgh',{'mgh','w'},...
  'presmooth',0,[0,1000],...
  'thresh',1,[0,1000],...
  'postsmooth',0,[0,1000],...
  'r0',4.5,[0,100],...
  'r1',5.5,[0,100],...
  'rmin',0.25,[0,100],...
  'rmax',12.5,[0,100],...
  'logtrans_flag',true,[false true],...
  'mask_frame',1,[1,Inf],...
  'surf','white',[],...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'forceflag',false,[false true],...
... % hidden parms
  'intype','mgh',{'mgh'},...
  'infixlist',{'_r','_i'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

if ~isempty(parms.subjdir)
  setenv('SUBJECTS_DIR',parms.subjdir)
else
  error('FreeSurfer SUBJECTS_DIR environment variable not set');
end;

if ~iscell(parms.hemilist), parms.hemilist = {parms.hemilist}; end;

if ~exist(parms.outdir,'dir')
  [success,msg] = mkdir(parms.outdir);
  if ~success
    error('failed to create output dir %s:\n%s',parms.outdir,msg);
  end;
end;

if isempty(parms.maskstem)
  if parms.logtrans_flag
    phase0 = log(parms.r0/parms.rmin)/log(parms.rmax/parms.rmin);
    phase1 = log(parms.r1/parms.rmin)/log(parms.rmax/parms.rmin);
  else
    phase0 = (parms.r0-parms.rmin)/(parms.rmax-parms.rmin);
    phase1 = (parms.r1-parms.rmin)/(parms.rmax-parms.rmin);
  end;
  phase0 = phase0 - floor(phase0); % subtract any cycles > 1
  phase1 = phase1 - floor(phase1);
  if phase0 < 0, phase0 = phase0 + 1; end;
  if phase1 < 0, phase1 = phase1 + 1; end;
  fprintf('%s: will select vertices with phases between %0.2f and %0.2f\n',...
    mfilename,phase0,phase1);
  phase0 = 2*pi*phase0;
  phase1 = 2*pi*phase1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input files
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};

  if ~isempty(parms.maskstem)
    fname = sprintf('%s/%s-%s.%s',...
      parms.maskdir,parms.maskstem,hemi,parms.intype);
    if ~exist(fname,'file'), error('file %s not found',fname); end;
  else
    fname = sprintf('%s/%s%s-%s.%s',...
      parms.eccdir,parms.eccstem,parms.infixlist{1},hemi,parms.intype);
    if ~exist(fname,'file'), error('file %s not found',fname); end;

    fname = sprintf('%s/%s%s-%s.%s',...
      parms.eccdir,parms.eccstem,parms.infixlist{2},hemi,parms.intype);
    if ~exist(fname,'file'), error('file %s not found',fname); end;
  end;

  fname = sprintf('%s/%s%s-%s.%s',...
    parms.poldir,parms.polstem,parms.infixlist{1},hemi,parms.intype);
  if ~exist(fname,'file'), error('file %s not found',fname); end;

  fname = sprintf('%s/%s%s-%s.%s',...
    parms.poldir,parms.polstem,parms.infixlist{2},hemi,parms.intype);
  if ~exist(fname,'file'), error('file %s not found',fname); end;
end;

for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};

  fname_out_r = sprintf('%s/%s%s-%s.%s',...
    parms.outdir,parms.outstem,parms.infixlist{1},hemi,parms.intype);
  fname_out_i = sprintf('%s/%s%s-%s.%s',...
    parms.outdir,parms.outstem,parms.infixlist{2},hemi,parms.intype);
  if ~exist(fname_out_r,'file') || ~exist(fname_out_i,'file') || parms.forceflag
    % load pol data
    fname = sprintf('%s/%s%s-%s.%s',...
      parms.poldir,parms.polstem,parms.infixlist{1},hemi,parms.intype);
    vals_r = fs_load_mgh(fname);

    fname = sprintf('%s/%s%s-%s.%s',...
      parms.poldir,parms.polstem,parms.infixlist{2},hemi,parms.intype);
    if ~exist(fname,'file'), error('file %s not found',fname); end;
    vals_i = fs_load_mgh(fname);

    if ~isempty(parms.maskstem)
      % load mask
      fname = sprintf('%s/%s-%s.%s',...
        parms.maskdir,parms.maskstem,hemi,parms.intype);
      vals_mask = fs_load_mgh(fname,[],parms.mask_frame);

      if length(vals_mask) ~= length(vals_r) ||...
         length(vals_mask) ~= length(vals_i)
        error('input data must have same number of values');
      end;
    else
      % load ecc data
      fname = sprintf('%s/%s%s-%s.%s',...
        parms.eccdir,parms.eccstem,parms.infixlist{1},hemi,parms.intype);
      vals_ecc_r = fs_load_mgh(fname);

      fname = sprintf('%s/%s%s-%s.%s',...
        parms.eccdir,parms.eccstem,parms.infixlist{2},hemi,parms.intype);
      if ~exist(fname,'file'), error('file %s not found',fname); end;
      vals_ecc_i = fs_load_mgh(fname);

      if length(vals_ecc_r) ~= length(vals_r) ||...
         length(vals_ecc_i) ~= length(vals_i) ||...
         length(vals_r) ~= length(vals_i)
        error('input data must have same number of values');
      end;
    end;

    % load surface if going to smooth
    if parms.presmooth>0 || parms.postsmooth>0
      surf = fs_load_subj(subj,hemi,parms.surf,0,parms.subjdir);
    end;

    % presmooth
    if parms.presmooth>0
      if ~isempty(parms.maskstem)
        vals_mask = fs_smooth(surf,vals_mask,parms.presmooth);
      else
        vals_ecc_r =  fs_smooth(surf,vals_ecc_r,parms.presmooth);
        vals_ecc_i =  fs_smooth(surf,vals_ecc_i,parms.presmooth);
      end;
    end;

    if ~isempty(parms.maskstem)
      % threshold
      vals_mask(vals_mask<parms.thresh) = 0;
      vals_mask(vals_mask>=parms.thresh) = 1;
    else
      vals_mask = hypot(vals_ecc_r,vals_ecc_i);
      vals_phase = atan2(vals_ecc_i,vals_ecc_r);
      % threshold      
      vals_mask(vals_mask<parms.thresh) = 0;
      vals_mask(vals_mask>=parms.thresh) = 1;
      % select phases
      vals_phase(vals_phase<0) = vals_phase(vals_phase<0) + 2*pi;
      if phase0 < phase1
        ind = find(vals_phase<phase0 | vals_phase>phase1);
      elseif phase1 > phase0
        ind = find(vals_phase<phase0 & vals_phase>phase1);
      else
        ind = find(vals_phase~=phase0);
      end;
      vals_mask(ind) = 0;
    end;

    % postsmooth
    if parms.postsmooth>0
      vals_mask = fs_smooth(surf,vals_mask,parms.postsmooth);
    end;

    % apply mask to pol data
    vals_r = vals_r.*vals_mask;
    vals_i = vals_i.*vals_mask;

    % save masked pol data
    fs_save_mgh(vals_r,fname_out_r);
    fs_save_mgh(vals_i,fname_out_i);
  end;
end;
