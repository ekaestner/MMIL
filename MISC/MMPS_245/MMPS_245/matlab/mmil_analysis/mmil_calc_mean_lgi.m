function mean_lgi = mmil_calc_mean_lgi(fspath,varargin)
%function mean_lgi = mmil_calc_mean_lgi(fspath,[options])
%
% Purpose: calculate mean cortical lgi (local gyrification index)
%
% Required Input:
%   fspath: FreeSurfer recon container path
%
% Optional Input:
%   'outdir': output directory
%     if not absolute path, will be relative to fspath
%     {default = 'analysis'}
%   'outstem': output file stem
%     if absolute path, outdir will be ignored
%     if empty, will use file stem of fname
%     {default = 'lgi'}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   mean_lgi: local gyrification index averaged over all cortical vertices
%
% Created:  09/29/13 by Don Hagler
% Last Mod: 09/29/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_lgi = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check parameters, check input files and directories
parms = check_input(fspath,varargin);

% convert lgi file from curv to mgz
parms = convert_from_curv(parms);

% calculate mean lgi value for cortical ROIs, average across hemispheres
mean_lgi = calc_mean_lgi(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fspath,options)
  % parse input arguments
  parms = mmil_args2parms(options,{...
    'fspath',fspath,[],...
...
    'outdir','analysis',[],...
    'outstem','lgi',[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
... % hidden parameters
    'hemilist',{'lh','rh'},{'lh','rh'},...
  });
  parms.nhemi = length(parms.hemilist);

  % check FreeSurfer recon exists
  if ~exist(parms.fspath,'dir')
    error('FreeSurfer recon dir %s not found',parms.fspath);
  end;
  [parms.subjdir,parms.subj,text] = fileparts(parms.fspath);
  parms.subj = [parms.subj text];

  % check input files
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    fname = sprintf('%s/surf/%s.pial_lgi',...
      parms.fspath,hemi);
    if ~exist(fname,'file')
      error('file %s not found',fname);
    end;
  end;

  % check cortex labels
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    fname_label = sprintf('%s/label/%s.cortex.label',...
      parms.fspath,hemi);
    if ~exist(fname_label,'file')
      error('file %s not found',fname_label);
    end;
  end;

  % set output directory and create if necessary
  if ~mmil_isrelative(parms.outstem)
    parms.outdir = fileparts(parms.outstem);
  else
    if mmil_isrelative(parms.outdir)
      parms.outdir = [parms.fspath '/' parms.outdir];
    end;
    parms.outstem = [parms.outdir '/' parms.outstem];
  end;
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = convert_from_curv(parms)
  meas = 'pial_lgi';
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    fname_in = sprintf('%s/surf/%s.pial',...
      parms.fspath,hemi);
    fname_out = sprintf('%s-%s.mgz',parms.outstem,hemi);
    if ~exist(fname_out,'file') || parms.forceflag
      cmd = sprintf('setenv SUBJECTS_DIR %s; mri_surf2surf --s %s --hemi %s',...
        parms.subjdir,parms.subj,hemi);
      cmd = sprintf('%s --srcsurfval %s --src_type curv --tval %s',...
        cmd,meas,fname_out);
      [s,r] = unix(cmd);
      if s
        error('cmd %s failed:\n%s',cmd,r);
      elseif parms.verbose
        fprintf('%s\n',r);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mean_lgi = calc_mean_lgi(parms)
  mean_lgi = 0;
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    fname_out = sprintf('%s-%s.mat',parms.outstem,hemi);
    if ~exist(fname_out,'file') || parms.forceflag
      fname_data = sprintf('%s-%s.mgz',parms.outstem,hemi);
      fname_label = sprintf('%s/label/%s.cortex.label',parms.fspath,hemi);
      results = mmil_surf_roi(fname_data,'fname_label',fname_label,...
        'minval',0,'verbose',parms.verbose);
      save(fname_out,'results');
    else
      load(fname_out);
    end;
    mean_lgi = mean_lgi + results.avg;
  end;
  mean_lgi = mean_lgi / parms.nhemi;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

