function fnames_out = mmil_paint(subj,fname_in,varargin);
%function fnames_out = mmil_paint(subj,fname_in,[options]);
%
% Purpose: resample volume data to cortical surface
%
% Usage:
%  mmil_paint(subj,fname_in,'key1', value1,...);
%
% Required parameters:
%  subj: FreeSurfer subject name
%  fname_in: full or relative path of input volume (mgh format)
%
% Optional parameters:
%  'hemi': should be either 'lh' for left hemisphere or 'rh' for right hemi
%    {default = both}
%  'outstem' : output file stem (omit extension, hemi)
%    {default = stem of fname_in}
%  'outdir': output directory
%    if empty, will use pwd or subjdir/analysis if meas supplied
%    {default = []}
%  'outfix' : add extra suffix to outstem before hemi and extension
%    {default = []}
%  'outext': output file extension ('.mgh' or '.mgz')
%    {default = '.mgh'}
%  'pmin': minimum distance (mm) to project along surface vertex normal
%    {default = 1}
%  'pmax': maximum distance (mm) to project along surface vertex normal
%    {default = 1}
%  'pstep': step distance (mm) along surface vertex normal
%    {default = 1}
%  'interpmethod': interpolation method ('nearest','linear','spline','cubic')
%    {default = 'linear'}
%  'avg_flag': [0|1] average values across projection steps
%    otherwise create multi-frame output volume
%    {default = 1}
%  'fname_weights': full path of volume file containing weighting factors
%    if supplied, and avg_flag=1, a weighted average will be calculated
%    assumed to be co-registered with fname_in
%    {deafult = []}
%  'tukey_flag': use Tukey's bisquare function to transform weights
%    {default = []}
%  'tukey_fact': scaling factor used to transform weights with Tukey's bisqaure
%    {default = 0.5}
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'surfname': surface to paint onto
%    {default = white}
%  'frame': frame number selected from fname_in
%    only a single frame will be painted to the surface
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Output:
%   fnames_out: cell array of output file names (e.g. left and right)
%
% Created:  08/21/15 by Don Hagler
% Last Mod: 09/22/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames_out = [];
if ~mmil_check_nargs(nargin,2), return; end;

% parse input parameters
parms = check_input(subj,fname_in,varargin);

% resample data from volume to surface
resample_vol2surf(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(subj,fname_in,options)

  parms = mmil_args2parms(options,{...
    'subj',subj,[],...
    'fname_in',fname_in,[],...
  ...
    'hemi',[],{'lh','rh'},...
    'outstem',[],[],...
    'outdir',[],[],...
    'outfix',[],[],...
    'outext','.mgh',{'.mgh','.mgz'},...
    'regmat',eye(4),[],...
    'pmin',1,[-10,10],...
    'pmax',1,[-10,10],...
    'pstep',1,[0.001,10],...
    'nsmooth',0,[1,10000],...
    'interpmethod','linear',{'nearest','linear','spline','cubic'},...
    'avg_flag',true,[false true],...
    'fname_weights',[],[],...
    'tukey_flag',false,[false true],...
    'tukey_fact',0.5,[1e-10,1e10],...
    'subjdir',[],[],...
    'surfname','white',[],...
    'frame',1,[1,Inf],...
    'forceflag',false,[false true],...
  ...
    'verbose',false,[false true],...
    'hemilist',{'lh','rh'},{'lh','rh'},...
  });

  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  if ~isempty(parms.fname_weights) && ~exist(parms.fname_weights,'file')
    error('file %s not found',parms.fname_weights);
  end;

  if isempty(parms.subjdir)
    parms.subjdir = getenv('SUBJECTS_DIR');
    if isempty(parms.subjdir)
      error('SUBJECTS_DIR not defined as an environment variable');
    end;
  else
    setenv('SUBJECTS_DIR',parms.subjdir);
  end;

  if ~isempty(parms.hemi), parms.hemilist = {parms.hemi}; end;
  parms.nhemi = length(parms.hemilist);
  parms.fnames_surf = cell(parms.nhemi,1);
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    parms.fnames_surf{h} = sprintf('%s/%s/surf/%s.%s',...
      parms.subjdir,parms.subj,hemi,parms.surfname);
    if ~exist(parms.fnames_surf{h},'file')
      error('surface file %s not found',parms.fnames_surf{h});
    end
  end;

  if any(size(parms.regmat)~=[4 4])
    error('regmat is wrong size');
  end;

  if length(parms.frame)>1
    error('frame must have single element');
  end;

  if isempty(parms.outstem)
    [tmp_path,tmp_fstem,tmp_ext] = fileparts(parms.fname_in);
    if isempty(parms.outdir)
      parms.outstem = [tmp_path '/' tmp_fstem];
    else
      parms.outstem = [parms.outdir '/' tmp_fstem];
    end;
  elseif mmil_isrelative(parms.outstem)
    if ~isempty(parms.outdir)
      parms.outstem = [parms.outdir '/' parms.outstem];
    end;
  end;
  if ~isempty(parms.outfix)
    parms.outstem = [parms.outstem parms.outfix];
  end;
  [outdir,tmp_fstem] = fileparts(parms.outstem);
  if isempty(outdir)
    outdir = pwd;
    parms.outstem = [outdir '/' parms.outstem];
  end;
  mmil_mkdir(outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnames_out = resample_vol2surf(parms)
  vol_vals = []; vol_wts = [];
  fnames_out = cell(parms.nhemi,1);
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    fnames_out{h} = sprintf('%s-%s%s',parms.outstem,hemi,parms.outext);
    if ~exist(fnames_out{h},'file') || parms.forceflag
      % load surf
      if parms.verbose
        fprintf('%s: loading surface file %s...\n',...
          mfilename,parms.fnames_surf{h});
      end;
      surf = fs_read_surf(parms.fnames_surf{h});

      % load vol_vals
      if isempty(vol_vals)
        if parms.verbose
          fprintf('%s: loading volume file %s...\n',...
            mfilename,parms.fname_in);
        end;
        [vol,M] = fs_load_mgh(parms.fname_in,[],parms.frame);
        vol_vals = ctx_mgh2ctx(vol,M);
        clear vol;
      end;
      
      % prep surface for use with ctx functions
      [surf,normvecs] = prep_surf(surf,vol_vals.lphcent,parms.regmat);
      %% todo: check orientation of normal vectors? -- need to use mask

      % sample vals from volume onto surface
      if parms.verbose
        fprintf('%s: sampling values...\n',mfilename);
      end;
      surf_vals = get_vals(vol_vals,surf,normvecs,parms);

      % calculate average across steps
      if parms.avg_flag
        % load weights
        if ~isempty(parms.fname_weights)
          if isempty(vol_wts)
            if parms.verbose
              fprintf('%s: loading weights file %s...\n',...
                mfilename,parms.fname_weights);
            end;
            vol_wts = ctx_load_mgh(parms.fname_weights);
            % transform weights using Tukey's bisquare function
            if parms.tukey_flag
              % scale weights by "Tukey" factor
              w = (1 - min(vol_wts.imgs,1)) / parms.tukey_fact;
              % apply Tukey's bisquare weight function
              w = (abs(w)<1) .* (1 - w.^2).^2;
              % replace values in vol_wts
              vol_wts.img = w;
            end;
          end;
          if parms.verbose
            fprintf('%s: sampling weights...\n',mfilename);
          end;
          surf_wts = get_vals(vol_wts,surf,normvecs,parms);
          surf_vals = surf_vals .* surf_wts;
          surf_vals = sum(surf_vals,1) ./ sum(surf_wts,1);
        else
          surf_vals = mean(surf_vals,1);
        end;
      end;

      % save surface vals
      if parms.verbose
        fprintf('%s: saving values to %s...\n',mfilename,fnames_out{h});
      end;
      surf_vals = reshape(surf_vals',[size(surf_vals,2),1,1,size(surf_vals,1)]);
      fs_save_mgh(surf_vals,fnames_out{h});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vals = get_vals(vol,surf,normvecs,parms)
  % sample values along normal vectors onto surface
  vals = getValsOnRay(vol,surf,normvecs,...
                      parms.pmin,parms.pmax,parms.pstep,parms.nsmooth,parms.interpmethod);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [surf,normvecs] = prep_surf(surf,lphcent,M)
  normvecs = [];
  ras2lph = [-1 -1 1];
  lph2ras = ras2lph;
  % add center coordinates
  surf.vertices = bsxfun(@plus,surf.vertices,lph2ras.*mmil_rowvec(lphcent));
  % convert from RAS to LPH
  surf.vertices = bsxfun(@times,surf.vertices,ras2lph);
  % transform surface coordinates with M
  V = cat(2,surf.vertices,ones(surf.nverts,1));
  V = (M * V')';
  surf.vertices = V(:,1:3);
  % preprocess surface into ctx format
  surf = preprocessQ(surf);
  % compute normal vectors
  normvecs = ts_calc_surf_normals(surf);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

