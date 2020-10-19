function mmil_analyze_fuzzy_icoarea(fspath,fsicopath,varargin)
%function mmil_analyze_fuzzy_icoarea(fspath,fsicopath,[options])
%
% Purpose: calculate average ico area for weighted ROIs
%
% Required Input:
%   fspath: freesurfer recon path
%   fsicopath: freesurfer ico resampled path
%
% Optional Parameters:
%   'outdir': where to place output files
%     can be full path, otherwise will be relative to FreeSurfer Container
%     {default = 'analysis'}
%   'verbose': [0|1] display status meassages
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default: 0}
%
%
% Created:  04/03/14 by Don Hagler
% Last Mod: 02/22/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin,{...
  'fspath',fspath,[],...
  'fsicopath',fsicopath,[],...
...
  'outdir','analysis',[],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
... % undocumented
  'hemilist',{'lh','rh'},{'lh' 'rh'},...
  'fuzzy_dir',[],[],...
  'fuzzy_fstem','fuzzy',[],...
  'fuzzy_order',18,[0,2,4,12,18],...
  'fuzzy_thresh',0,[0,Inf],...
  'fuzzy_smooth',2819,[0,Inf],...
  'outtype','mgz',{'mgh','mgz'},...
...
  'fuzzy_name_tags',{'fuzzy_fstem','fuzzy_order'},[],...
  'surf_roi_tags',{'fname_aparc','fname_label','fname_weights','frames',...
                   'minval','scalefact','fname_colorlut','hemi',...
                   'weights_thresh','verbose'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms.nhemi = length(parms.hemilist);

[parms.icosubjdir,parms.icosubj,text] = fileparts(parms.fsicopath);
parms.icosubj = [parms.icosubj text];

if mmil_isrelative(parms.outdir)
  parms.outdir = [parms.fspath '/' parms.outdir];
end;
mmil_mkdir(parms.outdir);

if isempty(parms.fuzzy_dir)
  parms.fuzzy_dir = [getenv('MMPS_DIR') '/atlases/fuzzy_clusters'];
end;
if ~exist(parms.fuzzy_dir,'dir')
  error('fuzzy cluster dir %s not found',parms.fuzzy_dir);
end;
parms.fnames_weights = cell(parms.nhemi,1);
parms.weights_thresh = parms.fuzzy_thresh;
parms.fuzzy_names = [];
if parms.fuzzy_order==0
  parms.fuzzy_outfix = parms.fuzzy_fstem;
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    % get roinames from txt file
    fname_txt = sprintf('%s/%s_%s_roinames.txt',...
      parms.fuzzy_dir,parms.fuzzy_fstem,hemi);
    if ~exist(fname_txt,'file')
      fname_txt = sprintf('%s/%s-%s-roinames.txt',...
        parms.fuzzy_dir,parms.fuzzy_fstem,hemi);
    end;
    if ~exist(fname_txt,'file')
      fprintf('%s: WARNING: roinames file %s not found\n',...
        mfilename,fname_txt);
      parms.fuzzy_names{h} = [];
    else
      parms.fuzzy_names{h} = mmil_readtext(fname_txt);
    end;
  end;
else
  args = mmil_parms2args(parms,parms.fuzzy_name_tags);
  for h=1:parms.nhemi
    parms.fuzzy_names{h} = mmil_fuzzy_names(args{:});
  end;
  parms.fuzzy_outfix = sprintf('%s%d',parms.fuzzy_fstem,parms.fuzzy_order);
end;
for h=1:parms.nhemi
  fname_weights = sprintf('%s/%s-%s.mgz',...
    parms.fuzzy_dir,parms.fuzzy_outfix,parms.hemilist{h});
  if ~exist(fname_weights,'file')
    error('fuzzy cluster file %s not found',fname_weights);
  end;
  parms.fnames_weights{h} = fname_weights;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parms.verbose
  fprintf('%s: calculating weighted averages for ico area...\n',mfilename);
  tic
end;
% convert icoarea to mgz
outstem = [parms.outdir '/icoarea'];
try
  fnames_out = fs_paint(parms.icosubj,[],...
    'meas','icoarea',...
    'outstem',outstem,...
    'outtype',parms.outtype,...
    'subjdir',parms.icosubjdir,...
    'smoothsteps',0,...
    'sphere_flag',0,...
    'sphsmoothsteps',parms.fuzzy_smooth,...
    'mask_midbrain_flag',0,...
    'hemilist',parms.hemilist,...
    'forceflag',parms.forceflag);
catch
  fprintf('\n%s: WARNING: converting icoarea to %s failed:\n%s\n\n',...
    mfilename,parms.outtype,lasterr);
  errcode = 1;
  return;
end;
% calculate weighted averages
fname_mat = sprintf('%s_%s_roi_data.mat',outstem,parms.fuzzy_outfix);
if ~exist(fname_mat,'file') || parms.forceflag
  roi_data = [];
  for h=1:parms.nhemi
    parms.hemi = parms.hemilist{h};
    parms.fname_weights = parms.fnames_weights{h};
    args = mmil_parms2args(parms,parms.surf_roi_tags);
    try
      tmp_roi_data = mmil_surf_roi(fnames_out{h},args{:});
    catch
      fprintf('\n%s: WARNING: fuzzy cluster analysis for icoarea failed:\n%s\n\n',...
        mfilename,lasterr);
      errcode = 1;
      return;
    end;
    for i=1:length(tmp_roi_data)
      tmp_roi_data(i).roiname = sprintf('ctx-%s-%s',...
        parms.hemi,parms.fuzzy_names{h}{i});
    end;
    roi_data = [roi_data,tmp_roi_data];
  end;
  save(fname_mat,'roi_data');
end;
if parms.verbose, toc; end;
