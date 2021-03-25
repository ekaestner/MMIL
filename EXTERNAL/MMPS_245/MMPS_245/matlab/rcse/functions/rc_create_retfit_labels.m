function [label_files,weights_files] = rc_create_retfit_labels(subj,varargin)
%function [label_files,weights_files] = rc_create_retfit_labels(subj,[options])
%
% Purpose: create label (and optional weights) files from retfit results
%   for all or part of visual areas
%
% Required Input:
%   subj: FreeSurfer recon name
%
% Optional Parameters:
%  'retfit_dir':  full path of retfit directory
%     {default = pwd}
%  'retfit_stem': full or absolute path file stem of retfit results
%     If not full path, will be relative to retfit_dir/matfiles
%     {default = 'retfit'}
%  'outdir': output directory
%     {default = pwd}
%  'outstem': output file stem
%     {default = 'retfit'}
%  'weights_flag': [0|1] in addition to creating label files, also save Gaussian
%     weights as mgh files
%    {default = 0}
%  'cond_info' - struct array containing condition information
%    see rc_read_cond_info
%    If cond_info and fname_conds are empty,
%      will create labels for entire visual areas
%    {default: []}
%  'fname_conds' - full path of csv file containing condition information
%    Ignored if cond_info is supplied
%    {default: []}
%  'r_vec' - vector of eccentricities (degrees visual angle)
%    Ignored if cond_info or fname_conds are supplied
%    {default = 7}
%  'th_vec' vector of polar angles (degrees)
%    Ignored if cond_info or fname_conds are supplied
%    {default = [45,135,225,315]}
%  'ecc_width' - eccentricity width of stimuli (deg. vis. ang.)
%    Ignored if cond_info has ecc_width column
%    {default = 10}
%  'theta_width' - polar angle width of stimuli (degrees)
%    Ignored if cond_info has theta_width column
%    will be ignored if fname_conds has theta_width column
%    {default = 90}
%  'r_max': maximum radius (degrees visual angle) used for eccentricity mapping
%    determines phase for a given eccentricity
%    {default = 12.5}
%  'area_names': names of visual areas
%    {default = {'v1','v2','v3'}}
%  'rf_sizes': vector of receptive field sizes (degrees visual angle)
%    -- one for each visual area
%    {default = [0.66,1.03,1.88]}
%  'rf_slopes': vector of slopes of linear trend of receptive field sizes
%    w.r.t. ecc for each visual area
%    Intercept is assumed to be half of r_max
%    {default = [0.06,0.10,0.15]}
%  'vf2ctx_flag': [0|1] precompute mapping from visual field to cortex
%    {default = 1}
%  'stimres': number of pixel elements in stimulus grid
%    a larger number allows smaller receptive field sizes to be used
%    {default = 100}
%  'w_thresh': threshold applied to weights (relative to max across all areas)
%    {default = 0.01}
%  'vfnorm_flag': [0|1] for applying w_thresh for vf2ctx
%     0: normalize to global max
%     1: normalize to max for each cortical location
%     {default = 1}
%  'retfit_thresh': threshold applied to weights (relative to max for one area)
%    {default = 0}
%  'surround_flag': [0|1] whether to model inhibitory surround with difference
%     of Gaussians; Requires vf2ctx_flag = 1
%    {default = 0}
%  'surround_rf_fact': size of surround receptive field relative to center
%    {default = 1.5}
%  'surround_amp_fact': amplitude of surround response relative to center
%    {default = 0.6}
%  'restrict_hemi_flag': [0|1] restrict weights to contralateral hemisphere
%    (otherwise allow weights in both; e.g. vertical meridian)
%    {default = 1}
%  'restrict_uplow_flag': [0|1] restrict weights to upper or lower field
%    portions of visual areas
%    (otherwise allow weights in both; e.g. horizontal meridian)
%    {default = 1}
%  'subjdir': FreeSurfer subject root dir
%    {default = getenv('SUBJECTS_DIR')}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  04/24/11 by Don Hagler
% Last Mod: 12/08/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms( varargin, {...
  'retfit_stem','retfit',[],...
  'retfit_dir',pwd,[],...
  'outstem','retfit',[],...
  'outdir',pwd,[],...
  'weights_flag',false,[false true],...
  'cond_info',[],[],...
  'fname_conds',[],[],...
  'r_vec',7,[],...
  'th_vec',[45,135,225,315],[],...
  'ecc_width',10,[0,100],...
  'theta_width',90,[0,360],...
  'r_max',12.5,[0,Inf],...
  'area_names',{'v1','v2','v3'},[],...
  'rf_sizes',[0.66,1.03,1.88],[],...
  'rf_slopes',[0.06,0.10,0.15],[0,10],...
  'vf2ctx_flag',true,[false true],...
  'stimres',100,[50,1000],...
  'w_thresh',0.01,[0,100],...
  'vfnorm_flag',true,[false true],...
  'retfit_thresh',0,[],...
  'surround_flag',false,[false true],...
  'surround_rf_fact',1.5,[],...
  'surround_amp_fact',0.6,[],...
  'restrict_hemi_flag',true,[false true],...
  'restrict_uplow_flag',true,[false true],...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'forceflag',false,[false,true],...
...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'contrasts',1,[],...
});

label_files = [];
weights_files = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.subjdir), error('subjdir is empty'); end;

if mmil_isrelative(parms.retfit_dir)
  parms.retfit_dir = [pwd '/' parms.retfit_dir];
end;
if mmil_isrelative(parms.retfit_stem)
  parms.retfit_stem = [parms.retfit_dir '/matfiles/' parms.retfit_stem];
end;

retfit_results = rc_load_retfit_results(parms.retfit_stem);

nverts = zeros(1,2);
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  surf = fs_load_subj(subj,hemi,[],1,parms.subjdir);
  nverts(h) = surf.nverts;
end;

% check all output files to avoid unnecessarily running rc_calc_vf2ctx
if exist(parms.outdir,'dir') && ~parms.forceflag
  [all_exist,label_files,weights_files] = check_output_files(parms);
  if all_exist, return; end;
else
  mmil_mkdir(parms.outdir);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define mapping between visual field and cortex
if parms.vf2ctx_flag
  tags = {'stimres' 'rf_sizes' 'rf_slopes' 'surround_flag'...
   'surround_rf_fact' 'surround_amp_fact' 'r_max' 'w_thresh' 'vfnorm_flag'...
   'area_names','hemilist'};
  args = mmil_parms2args(parms,tags);
  [vf2ctx,retfit_results] = rc_calc_vf2ctx(retfit_results,args{:});
else
  vf2ctx = [];
end;

% define patches for each stimulus for each area
tags = {'stimres' 'rf_sizes' 'rf_slopes' 'r_max' 'w_thresh'...
 'cond_info' 'fname_conds' 'r_vec' 'th_vec' 'ecc_width' 'theta_width'...
 'restrict_hemi_flag' 'restrict_uplow_flag' 'area_names','hemilist'};
args = mmil_parms2args(parms,tags);
retmap = rc_define_retmap_from_retfit(retfit_results,...
  'vf2ctx',vf2ctx,args{:});

% create labels for each area and hemi
for c=1:retmap.num_locs
  event_code = mmil_getfield(retmap.cond_info(c),'event_code',[]);
  if isempty(event_code)
    tmp_outstem = sprintf('%s_cond%d',parms.outstem,c);
  else
    tmp_outstem = sprintf('%s_event%d',parms.outstem,event_code);
  end;

  th = retmap.cond_info(c).theta;

  if parms.restrict_uplow_flag
    if th<180
      tmp_uplowlist = {'+'};
    else
      tmp_uplowlist = {'-'};
    end;
  else
    tmp_uplowlist = {'+','-'};
  end;

  if parms.restrict_hemi_flag
    if th>90 && th<270 % left visual field
      tmp_hemilist = {'rh'}; % left visual field, right cortical hemisphere
    else
      tmp_hemilist = {'lh'}; % right visual field, left cortical hemisphere
    end;
  else
    tmp_hemilist = {'lh','rh'};  
  end;

  for j=1:length(tmp_uplowlist)
    uplowstr = tmp_uplowlist{j};
    for k=1:length(tmp_hemilist)
      hemi = tmp_hemilist{k};
      h = find(strcmp(hemi,parms.hemilist));
      for a=1:retmap.num_areas
        fname_label = sprintf('%s/%s.%s_%s%s.label',...
          parms.outdir,hemi,tmp_outstem,retmap.areas(a).name,uplowstr);
        fname_weights = sprintf('%s/%s_%s%s-%s.mgh',...
          parms.outdir,tmp_outstem,retmap.areas(a).name,uplowstr,hemi);
        if ~exist(fname_label,'file') || ...
           (parms.weights_flag && ~exist(fname_weights,'file')) ||...
            parms.forceflag
          % get verts for current stim loc
          v = retmap.areas(a).verts(c).(['v_' hemi]);
          w = retmap.areas(a).verts(c).(['w_' hemi]);
          % normalize to max
          w = w/max(abs(w));
          % apply threshold
          ind = find(abs(w)>=parms.retfit_thresh);
          v = v(ind);
          w = w(ind);
          % save label file
          fs_write_label(v,fname_label,subj);
          if parms.weights_flag
            % save weights
            vals = zeros(nverts(h),1);
            vals(v) = w;
            fs_save_mgh(vals,fname_weights);
          end;
        end;
        label_files{end+1} = fname_label;
        if parms.weights_flag
          weights_files{end+1} = fname_weights;
        end;
      end;
    end;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [all_exist,label_files,weights_files] = check_output_files(parms)
  label_files = [];
  weights_files = [];
  all_exist = 1;

  if ~isempty(parms.cond_info) || ~isempty(parms.fname_conds)
    if ~isempty(parms.cond_info)
      cond_info = parms.cond_info;
    else
      cond_info = rc_read_cond_info(parms.fname_conds);
    end;
    if ~isempty(parms.contrasts) && isfield(cond_info,'contrast')
      contrasts = cell2mat({cond_info.contrast});
      cond_order = find(ismember(contrasts,parms.contrasts));
      cond_info = cond_info(cond_order);
    end;
  else
    c=1;
    for i=1:length(parms.r_vec)
      r = parms.r_vec(i);
      for j=1:length(parms.th_vec)
        th = parms.th_vec(j);
        if th<0, th=th+360; end;
        cond_info(c).ecc = r;
        cond_info(c).theta = th;
        c=c+1;
      end;
    end;
  end;
  num_locs = length(cond_info);
  num_areas = length(parms.area_names);
  
  for c=1:num_locs
    event_code = mmil_getfield(cond_info(c),'event_code',[]);
    if isempty(event_code)
      tmp_outstem = sprintf('%s_cond%d',parms.outstem,c);
    else
      tmp_outstem = sprintf('%s_event%d',parms.outstem,event_code);
    end;

    th = cond_info(c).theta;

    if parms.restrict_uplow_flag
      if th<180
        tmp_uplowlist = {'+'};
      else
        tmp_uplowlist = {'-'};
      end;
    else
      tmp_uplowlist = {'+','-'};
    end;

    if parms.restrict_hemi_flag
      if th>90 && th<270 % left visual field
        tmp_hemilist = {'rh'}; % left visual field, right cortical hemisphere
      else
        tmp_hemilist = {'lh'}; % right visual field, left cortical hemisphere
      end;
    else
      tmp_hemilist = {'lh','rh'};  
    end;

    for j=1:length(tmp_uplowlist)
      uplowstr = tmp_uplowlist{j};
      for k=1:length(tmp_hemilist)
        hemi = tmp_hemilist{k};
        h = find(strcmp(hemi,parms.hemilist));
        for a=1:num_areas
          fname_label = sprintf('%s/%s.%s_%s%s.label',...
            parms.outdir,hemi,tmp_outstem,parms.area_names{a},uplowstr);
          fname_weights = sprintf('%s/%s_%s%s-%s.mgh',...
            parms.outdir,tmp_outstem,parms.area_names{a},uplowstr,hemi);
          if ~exist(fname_label,'file')
            all_exist = 0;
            return;
          end;
          label_files{end+1} = fname_label;
          if parms.weights_flag
            if ~exist(fname_weights,'file')
              all_exist = 0;
              return;
            end;
            weights_files{end+1} = fname_weights;
          end;
        end;
      end;
    end;
  end;

return;

