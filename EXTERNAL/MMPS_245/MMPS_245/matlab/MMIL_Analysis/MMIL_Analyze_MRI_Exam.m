function errcode = MMIL_Analyze_MRI_Exam(FSContainerPath,varargin)
%function errcode = MMIL_Analyze_MRI_Exam(FSContainerPath,[options])
%
% Purpose:
%  Smooth surface stats such as thickness, area, T1w, extract ROI averages
%
% Usage:
%  MMIL_Analyze_MRI_Exam(FSContainerPath,'key1', value1,...)
%
% Required Parameters
%   FSContainerPath: full path of directory containing FreeSurfer recon
%
% Optional Parameters controlling which steps to run:
%   'run_all_flag': run all steps, overriding other control flags
%     {default = 0}
%   'aseg_flag': read freesurfer's aseg.stats, save in mat file
%     {default = 1}
%   'aparc_flag': read freesurfer's aparc.stats, save in mat files
%     {default = 1}
%   'thick_flag': resample thickness to sphere, smooth
%     if fname_weights supplied, also get weighted average of thickness
%     {default = 1}
%   'sulc_flag': resample sulc to sphere, smooth,
%       calculate aparc ROI averages
%     {default = 1}
%   'area_flag': resample area to sphere (with Jacobian correction), smooth
%     if fname_weights supplied, also get weighted average of area
%     {default = 1}
%   'cortvol_flag': calculate cortical volume from thickness and area
%     {default = 1}
%   'T1w_aseg_flag': extract average T1-weighted values (nu) in aseg ROIs
%     {default = 1}
%   'T1w_surf_flag': sample T1-weighted volume to surface,
%       resample to sphere, smooth, compute gray-white contrast,
%       calculate aparc ROI averages
%     {default = 1}
%   'proc_aseg_flag': extract average values in aseg ROIs
%       with input files determined by 'proc_inputlist'
%     {default = 1}
%   'proc_surf_flag': sample volume data to surface,
%       resample to sphere, smooth, compute gray-white contrast,
%       calculate aparc ROI averages
%       with input files determined by 'proc_inputlist'
%     {default = 1}
%   'fuzzy_flag': [0|1] use fuzzy cluster ROIs for thickness, area, and T1w
%     {default = 1}
%   'aseg_roigroups_flag': [0|1] create masks for groups of aseg roi codes
%       defined by aseg_roigroups
%     {default = 0}
%   'subhippo_flag': create subdivided hippocampal ROIs
%     {default = 0}
%   'erode_flag': create and use eroded ROIs (aseg, groups, hippo)
%     {default = 0}
%
% Other Optional Parameters:
%   'outdir': where to place output files
%     can be full path, otherwise will be relative to FreeSurfer Container
%     {default = 'analysis'}
%   'proc_inputlist': list of input file stems
%       relative to proc_indir e.g. 'MPR_res', 'T2w_res'
%       input files should be in register with orig.mgz
%         and have identical dimensions and vox2ras matrices
%       if empty, proc_aseg_flag and proc_surf_flag will be set to 0
%     {default = []}
%   'proc_outputlist': list of output file stems
%       corresponding to elements of proc_inputlist
%       if empty, will use file stems from proc_inputlist
%         e.g. {'T1w','T2w'}
%       if proc_outputlist contains 'T1w'
%         T1w_aseg_flag and T1w_surf_flag will be set to 0
%     {default = []}
%   'proc_indir': full path of directory containing processed data
%       if empty, proc_aseg_flag and proc_surf_flag will be set to 0
%     {default = []}
%   'proc_scalefacts': vector of scaling factors for
%       each element of proc_inputlist
%       if empty, will be set to 1 for each
%     {default = []}
%   'aparc_infix': string in aparc annot and stats file
%     e.g. 'aparc', 'aparc.a2009s'
%     {default = 'aparc'}
%   'fuzzy_dir': input directory for fuzzy cluster ROIs
%     if empty, will use [getenv('MMPS_dir') '/atlases/fuzzy_clusters']
%     {default = []}
%   'fuzzy_fstem': input file stem for fuzzy cluster ROIs
%     with expected names like {fstem}{order}-{hemi}.mgz
%     {default = 'fuzzy'}
%   'fuzzy_order': [0|2|4|12|18] number of fuzzy cluster ROIs
%     note: set of 18 includes combined sets of 2, 4, and 12
%     if order=0, names are like {fstem}-{hemi}.mgz
%     {default = 18}
%   'fuzzy_thresh': threshold applied to fuzzy cluster ROIs
%     {default = 0}
%   'fuzzy_smooth': smoothing applied to surface maps before extracting
%     values from fuzzy clusters
%     {default = 2819}
%   'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     {default = [-0.2,0.2]}
%   'smoothsteps': number of smoothing iterations on individual subject surface
%     {default = 0}
%   'sphere_flag': [0|1] whether to resample to spherical atlas
%     {default = 0}
%   'sphsmoothsteps': number of smoothing iterations on sphere (can be vector)
%     {default = 0}
%   'mask_midbrain_flag', [0|1] whether to mask out mid brain and other
%     cortical regions marked "unknown" (masking done before smoothing)
%     {default = 0}
%   'aseg_roigroups': struct array containing the following fields:
%      roiname: name of new ROI
%      roicode: new ROI code number
%      roicodes: vector of aseg ROI code numbers
%     {default: aseg_roigroups = fs_define_aseg_roigroups}
%       (includes 'WholeBrain', 'LatVentricles', and 'AllVentricles')
%   'erode_nvoxels': number of voxels to erode
%     {default = 1}
%   'check_stale_flag': check creation date of output directory
%     and delete if stale (created prior to fs.finish.all.touch)
%     {default = 0}
%   'check_complete_flag': [0|1] whether to require that recon is complete
%     {default = 1}
%   'FS_version': which version of Freesurfer used (e.g. 305, 450, 510, 530)
%     for checking whether recon is complete
%     {default = 530}
%   'verbose': [0|1] display status meassages
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default: 0}
%
% Created:  03/31/09 by Don Hagler
% Prev Mod: 12/02/16 by Don Hagler
% Last Mod: 07/11/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: aparc_roigroups for lobar analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

errcode = 0;

[parms,errcode] = check_input(FSContainerPath,varargin);
if errcode, return; end;

mmil_mkdir(parms.outdir);

% read aseg stats, save to mat file
if parms.aseg_flag
  errcode = errcode + read_aseg_stats(parms);
end;

% read aparc stats, save to mat file
if parms.aparc_flag
  errcode = errcode + read_aparc_stats(parms);
end;

% resample thickness to ico sphere surface, smooth
if parms.thick_flag
  errcode = errcode + smooth_surfmeas(parms,'thickness');
  % calculate average thickness for weighted ROIs
  if parms.fuzzy_flag
    errcode = errcode + analyze_fuzzy_surfmeas(parms,'thickness');
  end;
end;

% resample sulc to ico sphere surface, smooth
if parms.sulc_flag
  errcode = errcode + smooth_surfmeas(parms,'sulc');
  % calculate average sulc for aparc ROIs
  errcode = errcode + analyze_aparc_surfmeas(parms,'sulc');
  % calculate average sulc for weighted ROIs
  if parms.fuzzy_flag
    errcode = errcode + analyze_fuzzy_surfmeas(parms,'sulc');
  end;
end;

% resample area to ico sphere surface, smooth
if parms.area_flag
  errcode = errcode + smooth_surfmeas(parms,'area');
  % calculate average area for weighted ROIs
  if parms.fuzzy_flag
    errcode = errcode + analyze_fuzzy_surfmeas(parms,'area');
  end;
end;

% calculate cortical volume on ico, smooth
if parms.cortvol_flag
  errcode = errcode + smooth_cortvol(parms);
  % calculate average cortvol for weighted ROIs
  if parms.fuzzy_flag
    errcode = errcode + analyze_fuzzy_cortvol(parms);
  end;
end;

% sample T1 to surface, calculate gray-white contrast, smooth, ROI averages
if parms.T1w_surf_flag
  errcode = errcode + analyze_cortsurf_T1w(parms);
  % calculate average T1w for weighted ROIs
  if parms.fuzzy_flag
    errcode = errcode + analyze_fuzzy_T1w(parms);
  end;
end;

% sample values from proc files to surface, calculate gray-white contrast, smooth, ROI averages
if parms.proc_surf_flag
  errcode = errcode + analyze_cortsurf_proc(parms);
  % calculate average values for weighted ROIs
  if parms.fuzzy_flag
    errcode = errcode + analyze_fuzzy_proc(parms);
  end;
end;

% create subdivided hippocampal ROIs
if parms.subhippo_flag
  [tmp_errcode,parms.flist_hippo_masks] = subdiv_hippo(parms);
  errcode = errcode + tmp_errcode;
end;

% create aseg ROI group masks (e.g. WholeBrain, AllVentricles)
if parms.aseg_roigroups_flag
  [tmp_errcode,parms.aseg_roigroups] = create_roigroups(parms);
  errcode = errcode + tmp_errcode;
end;

% erode ROIs
if parms.erode_flag
  [tmp_errcode,parms.fname_aseg_erode] = erode_aseg(parms);
  if tmp_errcode, parms.T1w_aseg_flag = 0; end;
  if tmp_errcode, parms.proc_aseg_flag = 0; end;
  errcode = errcode + tmp_errcode;
end;

% extract average T1-weighted values for aseg ROIs
if parms.T1w_aseg_flag
  errcode = errcode + analyze_aseg_T1w(parms);
end;

% extract average values from proc files for aseg ROIs
if parms.proc_aseg_flag
  errcode = errcode + analyze_aseg_proc(parms);
end;

% display message about errors and/or finished
final_status(errcode,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_input(FSContainerPath,options)
  errcode = 0;
  parms = mmil_args2parms(options,{...
    'fspath',FSContainerPath,[],...
  ... % control flags:
    'run_all_flag',false,[false true],...
    'aseg_flag',true,[false true],...
    'aparc_flag',true,[false true],...
    'thick_flag',true,[false true],...
    'sulc_flag',true,[false true],...
    'area_flag',true,[false true],...
    'cortvol_flag',true,[false true],...
    'T1w_aseg_flag',true,[false true],...
    'T1w_surf_flag',true,[false true],...
    'proc_aseg_flag',true,[false true],...
    'proc_surf_flag',true,[false true],...
    'fuzzy_flag',true,[false true],...
    'aseg_roigroups_flag',false,[false true],...
    'subhippo_flag',false,[false true],...
    'erode_flag',false,[false true],...
  ...
    'outdir','analysis',[],...
    'proc_inputlist',[],[],...
    'proc_outputlist',[],[],...
    'proc_indir',[],[],...
    'proc_scalefacts',[],[],...
    'aparc_infix','aparc',[],...
    'fuzzy_dir',[],[],...
    'fuzzy_fstem','fuzzy',[],...
    'fuzzy_order',18,[0,2,4,12,18],...
    'fuzzy_thresh',0,[0,Inf],...
    'fuzzy_smooth',2819,[0,Inf],...
    'projdist_list',[-0.2,0.2],[-5,5],...
    'sphere_flag',false,[false true],...
    'smoothsteps',0,[0,Inf],...
    'sphsmoothsteps',0,[0,Inf],...
    'mask_midbrain_flag',false,[false true],...
    'aseg_roigroups',[],[],...
    'erode_nvoxels',1,[1:100],...
    'check_stale_flag',false,[false true],...
    'check_complete_flag',true,[false true],...
    'FS_version',530,[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % undocumented:
    'icosubj','fsaverage',[],...
    'icosubjdir',[getenv('FREESURFER_HOME') '/subjects'],[],...
    'hemilist',{'lh','rh'},{'lh' 'rh'},...
    'aseg_roilist',[2:5,7,8,10:18,24:26,28,41:44,46,47,49:54,57,58,60,77:79,251:255],[1,Inf],...
    'aseg_aparc_flag',0,[0,1,2],...
    'T1w_scalefact',1.0/256,[1e-10,1e10],...
    'outtype','mgz',{'mgh','mgz'},...
    'proc_intype','mgz',{'mgh','mgz'},...
  ...
    'aseg_tags',{'outdir','outstem','fname_out','csv_flag','fname_aseg',...
                 'aseg_aparc_flag','fname_vals','dispvec',...
                 'disp_roicodes','disp_scalefact','dispfact',...
                 'erode_flag','erode_nvoxels',...
                 'scalefact','minval','M_reg','res_outfix','fname_colorlut',...
                 'verbose','forceflag',...
                 'aseg_roilist','aparc_roilist','exclude_roilist','frames'},[],...
    'cortsurf_tags',{'outdir','outstem','outfix','fnames_aparc',...
                     'fnames_weights','weights_thresh',...
                     'csv_flag','M_reg','resT1flag','res_outfix',...
                     'projdist_list','gwnorm_flag',...
                     'smoothsteps','sphere_flag','sphsmoothsteps',...
                     'mask_midbrain_flag',...
                     'scalefact','minval','fname_colorlut',...
                     'verbose','forceflag'},[],...
    'surf_roi_tags',{'fname_aparc','fname_label','fname_weights','frames',...
                     'minval','scalefact','fname_colorlut','hemi',...
                     'weights_thresh','verbose'},[],...
    'fuzzy_name_tags',{'fuzzy_fstem','fuzzy_order'},[],...
  });

  if parms.run_all_flag
    parms.aseg_flag = 1;
    parms.aparc_flag = 1;
    parms.thick_flag = 1;
    parms.sulc_flag = 1;
    parms.area_flag = 1;
    parms.cortvol_flag = 1;
    parms.T1w_aseg_flag = 1;
    parms.T1w_surf_flag = 1;
    parms.proc_aseg_flag = 1;
    parms.proc_surf_flag = 1;
    parms.aseg_roigroups_flag = 1;
    parms.subhippo_flag = 1;
    parms.erode_flag = 1;
  end;

  % check fspath
  if ~exist(parms.fspath,'dir')
    error('FSContainerPath %s not found',parms.fspath);
  end;
  [parms.subjdir,parms.subj,text] = fileparts(parms.fspath);
  parms.subj = [parms.subj text];

  % check status of FreeSurfer recon
  if parms.check_complete_flag
    [status,message] = MMIL_Get_FSReconStatus(parms.fspath,parms.FS_version);
    if ~ismember(status,[2,5,6])
      fprintf('%s: WARNING: incomplete recon for %s\n',mfilename,parms.subj);
      errcode = 1;
      return;
    end;
    if status==6
      fprintf('%s: WARNING: only volume recon is complete for %s\n',...
        mfilename,parms.subj);
      parms.aparc_flag = 0;
      parms.thick_flag = 0;
      parms.sulc_flag = 0;
      parms.area_flag = 0;
      parms.cortvol_flag = 0;
      parms.T1w_surf_flag = 0;
      parms.proc_surf_flag = 0;
    end;
  end;
  
  % check aseg
  parms.fname_aseg = sprintf('%s/mri/aseg.mgz',parms.fspath);
  if ~exist(parms.fname_aseg)
    fprintf('%s: WARNING: aseg file %s not found\n',...
      mfilename,parms.fname_aseg);
    parms.aseg_flag = 0;
    parms.T1w_aseg_flag = 0;
    parms.proc_aseg_flag = 0;
    parms.aseg_roigroups_flag = 0;
    parms.subhippo_flag = 0;
    parms.erode_flag = 0;
  end;

  % check files necessary for surface analyses
  parms.nhemi = length(parms.hemilist);
  for h=1:parms.nhemi
    if parms.thick_flag || parms.sulc_flag || parms.area_flag ||...
       parms.cortvol_flag || parms.T1w_surf_flag
      spherefile = sprintf('%s/surf/%s.sphere.reg',...
        parms.fspath,parms.hemilist{h});
      if ~exist(spherefile,'file')
        fprintf('%s: WARNING: file %s not found\n',mfilename,spherefile);
        parms.thick_flag = 0;
        parms.sulc_flag = 0;
        parms.area_flag = 0;
        parms.cortvol_flag = 0;
        parms.T1w_surf_flag = 0;
        continue;
      end;
    end;
    if parms.thick_flag
      fname_thick = sprintf('%s/surf/%s.thickness',...
        parms.fspath,parms.hemilist{h});
      if ~exist(fname_thick,'file')
        fprintf('%s: WARNING: file %s not found\n',mfilename,fname_thick);
        parms.thick_flag = 0;
        parms.cortvol_flag = 0;
      end;
    end;
    if parms.sulc_flag
      fname_sulc = sprintf('%s/surf/%s.sulc',...
        parms.fspath,parms.hemilist{h});
      if ~exist(fname_sulc,'file')
        fprintf('%s: WARNING: file %s not found\n',mfilename,fname_sulc);
        parms.sulc_flag = 0;
      end;
    end;
    if parms.area_flag
      fname_area = sprintf('%s/surf/%s.area',...
        parms.fspath,parms.hemilist{h});
      if ~exist(fname_area,'file')
        fprintf('%s: WARNING: file %s not found\n',mfilename,fname_area);
        parms.area_flag = 0;
        parms.cortvol_flag = 0;
      end;
    end;
  end;

  % set fnames_aparc based on parms.aparc_infix
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    parms.fnames_aparc{h} = sprintf('%s/label/%s.%s.annot',...
      parms.fspath,hemi,parms.aparc_infix);
  end;

  % set fnames_weights if fuzzy_flag
  if parms.fuzzy_flag
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
  end;

  % check that T1-weighted volume exists
  if parms.T1w_aseg_flag || parms.T1w_surf_flag
    parms.fname_T1 = sprintf('%s/mri/nu.mgz',parms.fspath);
    if ~exist(parms.fname_T1,'file')
      fprintf('%s: WARNING: T1-weighted file %s not found\n',...
        mfilename,parms.fname_T1);
      parms.T1w_aseg_flag = 0;
      parms.T1w_surf_flag = 0;
    end;
  end;
  
  % check proc_indir and proc_inputlist
  if parms.proc_aseg_flag || parms.proc_surf_flag
    if isempty(parms.proc_inputlist)
      fprintf('%s: WARNING: proc_inputlist is empty\n',mfilename);
      parms.proc_aseg_flag = 0;
      parms.proc_surf_flag = 0;
    elseif isempty(parms.proc_indir)
      fprintf('%s: WARNING: proc_indir is empty\n',mfilename);
      parms.proc_aseg_flag = 0;
      parms.proc_surf_flag = 0;
    elseif ~exist(parms.proc_indir)
      fprintf('%s: WARNING: proc_indir %s does not exist\n',mfilename,parms.proc_indir);
      parms.proc_aseg_flag = 0;
      parms.proc_surf_flag = 0;
    else
      % check proc_inputlist
      ninputs = length(parms.proc_inputlist);
      if ~iscell(parms.proc_inputlist)
        parms.proc_inputlist = {parms.proc_inputlist};
      end;
      % check proc_outputlist
      if isempty(parms.proc_outputlist)
        parms.proc_outputlist = parms.proc_inputlist;
      elseif ~iscell(parms.proc_outputlist)
        parms.proc_outputlist = {parms.proc_outputlist};
      end;
      if length(parms.proc_outputlist) ~= ninputs
        error('number of elements in proc_outputlist does not match proc_inputlist');
      end;
      if isempty(parms.proc_scalefacts)
        parms.proc_scalefacts = ones(ninputs,1);
      elseif length(parms.proc_scalefacts) ~= ninputs
        error('number of elements in proc_scalefacts does not match proc_inputlist');
      end;
      % disable T1w analysis if proc_outputlist contains 'T1w'
      if ismember('T1w',parms.proc_outputlist)
        parms.T1w_aseg_flag = 0;
        parms.T1w_surf_flag = 0;
      end;
      % check existence of each file in proc_inputlist
      keep_flags = zeros(ninputs,1);
      for i=1:ninputs
        fname_in = sprintf('%s/%s.%s',...
          parms.proc_indir,parms.proc_inputlist{i},parms.proc_intype);
        if ~exist(fname_in,'file')
          fprintf('%s: WARNING: file %s not found\n',mfilename,fname_in);
        else
          keep_flags(i) = 1;
          parms.proc_inputlist{i} = fname_in;
        end;
      end;
      ind_keep = find(keep_flags);
      parms.proc_inputlist = parms.proc_inputlist(ind_keep);
      parms.proc_outputlist = parms.proc_outputlist(ind_keep);
      if isempty(ind_keep)
        parms.proc_aseg_flag = 0;
        parms.proc_surf_flag = 0;
      end;
    end;
  end;

  % check whether any thing to do
  if ~parms.aseg_flag && ~parms.aparc_flag &&...
     ~parms.thick_flag && ~parms.sulc_flag &&...
     ~parms.area_flag && ~parms.cortvol_flag &&...
     ~parms.T1w_aseg_flag && ~parms.T1w_surf_flag &&...
     ~parms.proc_aseg_flag && ~parms.proc_surf_flag &&...
     ~parms.aseg_roigroups_flag && ~parms.subhippo_flag && ~parms.erode_flag
     fprintf('%s: WARNING: nothing to do\n',mfilename);
     return;
  end;
  if mmil_isrelative(parms.outdir)
    parms.outdir = [parms.fspath '/' parms.outdir];
  end;
  mmil_mkdir(parms.outdir);

  % check timestamp of final or surf touch file, compare to timestamp of outdir
  if parms.check_stale_flag
    fstem_touch = 'fs.finish.all';
    fname_touch = sprintf('%s/touch/%s.touch',parms.fspath,fstem_touch);
    if ~exist(fname_touch,'file')
      fprintf('%s: WARNING: %s not found for %s (unfinished recon)\n',...
        mfilename,fstem_touch,parms.subj);
    else  
      dir_touch = dir(fname_touch);
      dir_outdir = dir(parms.outdir);
      date_touch = dir_touch.datenum;
      date_outdir = dir_outdir(1).datenum;
      if date_outdir < date_touch
        fprintf('%s: WARNING: fs.finish.all.touch is more recent than %s, removing existing output...\n',...
          mfilename,parms.outdir);
        [s,r] = unix(sprintf('rm -r %s',parms.outdir));
        if s, error('failed to remove outdir %s:\n%s',parms.outdir,r); end;
      end;
    end;
  end;

  % define sets of aseg ROIs to create new ROIs
  if isempty(parms.aseg_roigroups)
    parms.aseg_roigroups = fs_define_aseg_roigroups;
  end;
  
  % set outfix for eroded ROI files
  if parms.erode_flag
    parms.erode_outfix = 'erode';
    if parms.erode_nvoxels>1
      parms.erode_outfix = sprintf('%s%d',...
        parms.erode_outfix,parms.erode_nvoxels);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% read aseg stats, save to mat file
function errcode = read_aseg_stats(parms)
  errcode = 0;
  fname_mat = sprintf('%s/aseg_stats.mat',parms.outdir);
  if ~exist(fname_mat,'file') | parms.forceflag
    if parms.verbose
      fprintf('%s: extracting aseg stats...\n',mfilename);
    end;
    aseg_stats = fs_read_aseg_stats(parms.subj,parms.subjdir,...
      parms.aseg_roigroups);
    if ~isempty(aseg_stats)
      save(fname_mat,'aseg_stats');
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read aparc stats, save to mat file
function errcode = read_aparc_stats(parms)
  errcode = 0;
  fname_mat = sprintf('%s/%s_stats.mat',parms.outdir,parms.aparc_infix);
  if ~exist(fname_mat,'file') | parms.forceflag
    if parms.verbose
      fprintf('%s: extracting aparc stats...\n',mfilename);
    end;
    aparc_stats = ...
      fs_read_aparc_stats(parms.subj,parms.subjdir,parms.aparc_infix);
    if ~isempty(aparc_stats), save(fname_mat,'aparc_stats'); end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resample surface-based measure onto ico sphere surface, smooth on sphere
function errcode = smooth_surfmeas(parms,meas)
  errcode = 0;
  if parms.verbose
    fprintf('%s: resampling %s to ico, applying smoothing kernels...\n',...
      mfilename,meas);
    tic
  end;
  for sphsmoothsteps = parms.sphsmoothsteps
    outstem = sprintf('%s/%s',parms.outdir,meas);
    try
      fs_paint(parms.subj,[],...
        'outstem',outstem,...
        'meas',meas,...
        'outtype',parms.outtype,...
        'subjdir',parms.subjdir,...
        'smoothsteps',parms.smoothsteps,...
        'sphere_flag',1,...
        'sphsmoothsteps',sphsmoothsteps,...
        'mask_midbrain_flag',parms.mask_midbrain_flag,...
        'hemilist',parms.hemilist,...
        'forceflag',parms.forceflag);
    catch
      fprintf('\n%s: WARNING: smoothing %s failed:\n%s\n\n',...
        mfilename,meas,lasterr);
      errcode = 1;
      return;
    end;
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate average surface-based measure for aparc ROIs
function errcode = analyze_aparc_surfmeas(parms,meas)
  errcode = 0;
  if parms.verbose
    fprintf('%s: calculating ROI averages for %s...\n',mfilename,meas);
    tic
  end;
  outstem = sprintf('%s/%s',parms.outdir,meas);
  fname_mat = sprintf('%s_%s_roi_data.mat',outstem,parms.aparc_infix);
  if ~exist(fname_mat,'file') || parms.forceflag
    roi_data = [];
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      parms.fname_aparc = parms.fnames_aparc{h};
      args = mmil_parms2args(parms,parms.surf_roi_tags);
      fname_vals = sprintf('%s-%s.%s',outstem,hemi,parms.outtype);
      try
        tmp_roi_data = mmil_surf_roi(fname_vals,args{:});
      catch
        fprintf('\n%s: WARNING: aparc analysis for %s failed:\n%s\n\n',...
          mfilename,meas,lasterr);
        errcode = 1;
        return;
      end;
      roi_data = [roi_data,tmp_roi_data];
    end;
    save(fname_mat,'roi_data');
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate average surface-based measure for weighted ROIs
function errcode = analyze_fuzzy_surfmeas(parms,meas)
  errcode = 0;
  if parms.verbose
    fprintf('%s: calculating weighted averages for %s...\n',mfilename,meas);
    tic
  end;
  % resample meas to ico, convert to mgh/mgz
  outstem = sprintf('%s/%s',parms.outdir,meas);
  try
    fnames_out = fs_paint(parms.subj,[],...
      'meas',meas,...
      'outstem',outstem,...
      'outtype',parms.outtype,...
      'subjdir',parms.subjdir,...
      'smoothsteps',0,...
      'sphere_flag',1,...
      'sphsmoothsteps',parms.fuzzy_smooth,...
      'mask_midbrain_flag',parms.mask_midbrain_flag,...
      'hemilist',parms.hemilist,...
      'forceflag',parms.forceflag);
  catch
    fprintf('\n%s: WARNING: resampling %s to ico failed:\n%s\n\n',...
      mfilename,meas,lasterr);
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
        fprintf('\n%s: WARNING: fuzzy cluster analysis for %s failed:\n%s\n\n',...
          mfilename,meas,lasterr);
        errcode = 1;
        return;
      end;
      tmp_roi_data = ...
        replace_fuzzy_names(tmp_roi_data,parms.fuzzy_names{h},parms.hemi);
      roi_data = [roi_data,tmp_roi_data];
    end;
    save(fname_mat,'roi_data');
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate cortical volume from thickness and area and smooth
function errcode = smooth_cortvol(parms)
  errcode = 0;
  if parms.verbose
    fprintf('%s: calculating cortical volume...\n',mfilename);
    tic
  end;
  % resample thickness to ico, convert to mgh/mgz, no smoothing
  outstem = [parms.outdir '/thickness'];
  try
    fnames_thick = fs_paint(parms.subj,[],...
      'meas','thickness',...
      'outstem',outstem,...
      'outtype',parms.outtype,...
      'subjdir',parms.subjdir,...
      'smoothsteps',0,...
      'sphere_flag',1,...
      'sphsmoothsteps',0,...
      'mask_midbrain_flag',parms.mask_midbrain_flag,...
      'hemilist',parms.hemilist,...
      'forceflag',parms.forceflag);
  catch
    fprintf('\n%s: WARNING: resampling thickness to ico failed:\n%s\n\n',...
      mfilename,lasterr);
    errcode = 1;
    return;
  end;
  % resample area to ico, convert to mgh/mgz, no smoothing
  outstem = [parms.outdir '/area'];
  try
    fnames_area = fs_paint(parms.subj,[],...
      'meas','area',...
      'outstem',outstem,...
      'outtype',parms.outtype,...
      'subjdir',parms.subjdir,...
      'smoothsteps',0,...
      'sphere_flag',1,...
      'sphsmoothsteps',0,...
      'mask_midbrain_flag',parms.mask_midbrain_flag,...
      'hemilist',parms.hemilist,...
      'forceflag',parms.forceflag);
  catch
    fprintf('\n%s: WARNING: resampling area to ico %s failed:\n%s\n\n',...
      mfilename,parms.outtype,lasterr);
    errcode = 1;
    return;
  end;
  
  % calculate cortical volume from thickness and area
  outstem = sprintf('%s/cortvol-sphere',parms.outdir);
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    fname_thick = fnames_thick{h};
    fname_area = fnames_area{h};
    fname_vol = sprintf('%s-%s.%s',outstem,hemi,parms.outtype);    
    if ~exist(fname_vol,'file')
      vals_thick = fs_load_mgh(fname_thick);
      vals_area = fs_load_mgh(fname_area);
      vals_vol = vals_thick .* vals_area;
      fs_save_mgh(vals_vol,fname_vol);
    end;
  end;
  if parms.verbose
    toc
    fprintf('%s: smoothing cortical volume on ico...\n',mfilename);
    tic
  end;
  if parms.fuzzy_flag && parms.fuzzy_smooth>0
    smoothvec = union(parms.sphsmoothsteps,parms.fuzzy_smooth);
  else
    smoothvec = parms.sphsmoothsteps;
  end;
  smoothvec = setdiff(smoothvec,0);
  for s=1:length(smoothvec)
    smooth = smoothvec(s);
    try
      % smooth cortvol
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        fname_in = sprintf('%s-%s.%s',outstem,hemi,parms.outtype);
        fname_out = sprintf('%s-sm%d-%s.%s',...
          outstem,smooth,hemi,parms.outtype);
        fs_surf2surf(fname_in,parms.icosubj,...
          'fname_out',fname_out,...
          'hemi',hemi,...
          'outtype',parms.outtype,...
          'smooth_in',smooth,...
          'subjdir',parms.icosubjdir,...
          'verbose',0,...
          'forceflag',parms.forceflag);
      end;
    catch
      fprintf('\n%s: WARNING: smoothing cortical volume failed:\n%s\n\n',...
        mfilename,lasterr);
      errcode = 1;
      return;
    end;
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate average cortvol for weighted ROIs
function errcode = analyze_fuzzy_cortvol(parms)
  errcode = 0;
  if parms.verbose
    fprintf('%s: calculating weighted averages for cortical volume...\n',...
      mfilename);
    tic
  end;
  instem = sprintf('%s/cortvol-sphere',parms.outdir);
  outstem = [parms.outdir '/cortvol'];
  % calculate weighted averages
  fname_mat = sprintf('%s_%s_roi_data.mat',outstem,parms.fuzzy_outfix);
  if ~exist(fname_mat,'file') || parms.forceflag
    roi_data = [];
    for h=1:parms.nhemi
      parms.hemi = parms.hemilist{h};
      parms.fname_weights = parms.fnames_weights{h};
      if parms.fuzzy_smooth>0
        fname_out = sprintf('%s-sm%d-%s.%s',...
          instem,parms.fuzzy_smooth,parms.hemi,parms.outtype);
      else
        fname_out = sprintf('%s-%s.%s',instem,parms.hemi,parms.outtype);    
      end;
      args = mmil_parms2args(parms,parms.surf_roi_tags);
      try
        tmp_roi_data = mmil_surf_roi(fname_out,args{:});
      catch
        fprintf('\n%s: WARNING: fuzzy cluster analysis for cortvol failed:\n%s\n\n',...
          mfilename,lasterr);
        errcode = 1;
        return;
      end;
      tmp_roi_data = ...
        replace_fuzzy_names(tmp_roi_data,parms.fuzzy_names{h},parms.hemi);
      roi_data = [roi_data,tmp_roi_data];
    end;
    save(fname_mat,'roi_data');
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample T1 to surface, calculate gray-white contrast, smooth, ROI averages
function errcode = analyze_cortsurf_T1w(parms)
  errcode = 0;
  if parms.verbose
    fprintf('%s: T1w cortsurf analysis for %s...\n',mfilename,parms.fname_T1);
    tic;
  end;
  parms.outstem = [parms.outdir '/T1w'];
  parms.fnames_weights = [];
  parms.scalefact = parms.T1w_scalefact;
  parms.outfix = parms.aparc_infix;
  args = mmil_parms2args(parms,parms.cortsurf_tags);
  try
    mmil_cortsurf_analysis(parms.fname_T1,parms.fspath,args{:});
  catch
    fprintf('\n%s: WARNING: T1w cortsurf analysis failed:\n%s\n\n',...
      mfilename,lasterr);
    errcode = 1;
    return;
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample values from proc files to surface,
%   calculate gray-white contrast, smooth, ROI averages
function errcode = analyze_cortsurf_proc(parms)
  errcode = 0;
  for i=1:length(parms.proc_inputlist)
    fname_in = parms.proc_inputlist{i};
    if parms.verbose
      fprintf('%s: proc cortsurf analysis for %s...\n',mfilename,fname_in);
      tic;
    end;
    parms.outstem = sprintf('%s/%s',parms.outdir,parms.proc_outputlist{i});
    parms.fnames_weights = [];
    parms.scalefact = parms.proc_scalefacts(i);
    parms.outfix = parms.aparc_infix;
    args = mmil_parms2args(parms,parms.cortsurf_tags);
    try
      mmil_cortsurf_analysis(fname_in,parms.fspath,args{:});
    catch
      fprintf('\n%s: WARNING: proc cortsurf analysis failed:\n%s\n\n',...
        mfilename,lasterr);
      errcode = 1;
      return;
    end;
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate average T1w for weighted ROIs
function errcode = analyze_fuzzy_T1w(parms)
  errcode = 0;
  if parms.verbose
    fprintf('%s: T1w cortsurf weighted analysis for %s...\n',...
      mfilename,parms.fname_T1);
    tic;
  end;
  % sample T1w to surface and resample to sphere without smoothing
  parms.scalefact = parms.T1w_scalefact;
  for projdist=parms.projdist_list
    outstem = sprintf('%s/T1w_pdist%0.1f',parms.outdir,projdist);
    try
      fnames_out = fs_paint(parms.subj,parms.fname_T1,...
        'outstem',outstem,...
        'outtype',parms.outtype,...
        'subjdir',parms.subjdir,...
        'projdist',projdist,...
        'smoothsteps',0,...
        'sphere_flag',1,...
        'sphsmoothsteps',parms.fuzzy_smooth,...
        'mask_midbrain_flag',parms.mask_midbrain_flag,...
        'hemilist',parms.hemilist,...
        'forceflag',parms.forceflag);
    catch
      fprintf('\n%s: WARNING: resampling T1w to sphere failed:\n%s\n\n',...
        mfilename,lasterr);
      errcode = 1;
      return;
    end;
    % calculate weighted averages
    fname_mat = sprintf('%s/T1w_%s_pdist%0.1f_roi_data.mat',...
      parms.outdir,parms.fuzzy_outfix,projdist);
    if ~exist(fname_mat,'file') || parms.forceflag
      roi_data = [];
      for h=1:parms.nhemi
        parms.hemi = parms.hemilist{h};
        parms.fname_weights = parms.fnames_weights{h};
        args = mmil_parms2args(parms,parms.surf_roi_tags);
        try
          tmp_roi_data = mmil_surf_roi(fnames_out{h},args{:});
        catch
          fprintf('\n%s: WARNING: fuzzy cluster analysis for T1w failed:\n%s\n\n',...
            mfilename,lasterr);
          errcode = 1;
          return;
        end;
        tmp_roi_data = ...
          replace_fuzzy_names(tmp_roi_data,parms.fuzzy_names{h},parms.hemi);
        roi_data = [roi_data,tmp_roi_data];
      end;
      save(fname_mat,'roi_data');
    end;
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate average proc for weighted ROIs
function errcode = analyze_fuzzy_proc(parms)
  errcode = 0;
  for i=1:length(parms.proc_inputlist)
    fname_in = parms.proc_inputlist{i};
    if parms.verbose
      fprintf('%s: proc fuzzy analysis for %s...\n',mfilename,fname_in);
      tic;
    end;
    % sample values from proc files to surface
    %   and resample to sphere without smoothing
    parms.scalefact = parms.proc_scalefacts(i);
    for projdist=parms.projdist_list
      outstem = sprintf('%s/%s_pdist%0.1f',...
        parms.outdir,parms.proc_outputlist{i},projdist);
      try
        fnames_out = fs_paint(parms.subj,fname_in,...
          'outstem',outstem,...
          'outtype',parms.outtype,...
          'subjdir',parms.subjdir,...
          'projdist',projdist,...
          'smoothsteps',0,...
          'sphere_flag',1,...
          'sphsmoothsteps',parms.fuzzy_smooth,...
          'mask_midbrain_flag',parms.mask_midbrain_flag,...
          'hemilist',parms.hemilist,...
          'forceflag',parms.forceflag);
      catch
        fprintf('\n%s: WARNING: resampling proc values to sphere failed:\n%s\n\n',...
          mfilename,lasterr);
        errcode = 1;
        return;
      end;
      % calculate weighted averages
      fname_mat = sprintf('%s/%s_%s_pdist%0.1f_roi_data.mat',...
        parms.outdir,parms.proc_outputlist{i},parms.fuzzy_outfix,projdist);
      if ~exist(fname_mat,'file') || parms.forceflag
        roi_data = [];
        for h=1:parms.nhemi
          parms.hemi = parms.hemilist{h};
          parms.fname_weights = parms.fnames_weights{h};
          args = mmil_parms2args(parms,parms.surf_roi_tags);
          try
            tmp_roi_data = mmil_surf_roi(fnames_out{h},args{:});
          catch
            fprintf('\n%s: WARNING: fuzzy cluster analysis for proc failed:\n%s\n\n',...
              mfilename,lasterr);
            errcode = 1;
            return;
          end;
          tmp_roi_data = ...
            replace_fuzzy_names(tmp_roi_data,parms.fuzzy_names{h},parms.hemi);
          roi_data = [roi_data,tmp_roi_data];
        end;
        save(fname_mat,'roi_data');
      end;
    end;
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create subdivided hippocampal ROIs
function [errcode,flist_hippo_masks] = subdiv_hippo(parms)
  errcode = 0;
  flist_hippo_masks = [];
  if parms.verbose
    fprintf('%s: creating subdivided hippocampus masks...\n',mfilename);
    tic
  end;
  flist_hippo_masks = fs_create_subhippo_masks(...
    parms.subj,parms.subjdir,parms.outdir,[],parms.forceflag);
  if parms.verbose, toc; end;
  if isempty(flist_hippo_masks)
    fprintf('%s: WARNING: fs_create_subhippo_masks failed\n',mfilename);
    errcode = 1;
    return;
  end;

  % erode hippo ROIs
  if parms.erode_flag
    if parms.verbose
      fprintf('%s: eroding hippo subdivision ROIs...\n',mfilename);
      tic
    end;
    for j=1:length(flist_hippo_masks)
      fname_in = flist_hippo_masks{j};
      [tmp_fpath,tmp_fstem,tmp_fext]=fileparts(flist_hippo_masks{j});
      fname_out = sprintf('%s/%s_%s%s',...
        parms.outdir,tmp_fstem,parms.erode_outfix,tmp_fext);
      fs_erode_mask(fname_in,fname_out,...
        'nvoxels',parms.erode_nvoxels,'forceflag',parms.forceflag);
      flist_hippo_masks{j} = fname_out;
    end;
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create aseg roi group masks (e.g. WholeBrain, LatVentricles, AllVentricles)
function [errcode,aseg_roigroups] = create_roigroups(parms)
  errcode = 0;
  aseg_roigroups = parms.aseg_roigroups;
  if parms.verbose
    fprintf('%s: creating aseg group ROIs...\n',mfilename);
    tic
  end;
  for j=1:length(aseg_roigroups)
    fname_mask = sprintf('%s/aseg_roigroup_%d.mgz',...
      parms.outdir,aseg_roigroups(j).roicode);
    if ~exist(fname_mask) | parms.forceflag
      [vol, M, mr_parms] = fs_load_mgh(parms.fname_aseg);
      % find voxels with matching roi code numbers
      ind = ismember(vol(:),aseg_roigroups(j).roicodes);
      % set those voxels to the new roi code number, others to 0
      vol(ind) = aseg_roigroups(j).roicode;
      vol(~ind) = 0;
      fs_save_mgh(vol,fname_mask,M,mr_parms);
    end;
    aseg_roigroups(j).fname = fname_mask;
  end
  if parms.verbose, toc; end;

  % create eroded WholeBrain, AllVentricles ROIs
  if parms.erode_flag
    if parms.verbose
      fprintf('%s: eroding aseg group ROIs...\n',mfilename);
      tic
    end;
    for j=1:length(aseg_roigroups)
      fname_in = aseg_roigroups(j).fname;
      fname_out = sprintf('%s/aseg_roigroup_%d_%s.mgz',...
        parms.outdir,aseg_roigroups(j).roicode,parms.erode_outfix);
      if ~exist(fname_out,'file') || parms.forceflag
        % load aseg
        [vol, M, mr_parms] = fs_load_mgh(parms.fname_aseg);
        % find ROI's in the ROI group in the aseg
        ind = ismember(vol(:),aseg_roigroups(j).roicodes);
        % set all ROI #'s in this group to group ROI #
        vol(ind) = aseg_roigroups(j).roicode;
        % set remaining to 0 to speed up fs_erode_aseg since they're not needed
        vol(~ind) = 0;
        % save modified aseg w/ only the current ROI group
        fs_save_mgh(vol,fname_in, M, mr_parms);
        % erode the ROI group
        fs_erode_aseg(fname_in,fname_out,...
          'nvoxels',parms.erode_nvoxels,'forceflag',parms.forceflag);
      end;
      aseg_roigroups(j).fname = fname_out;
    end
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% erode ROIs
function [errcode,fname_aseg_erode] = erode_aseg(parms)
  errcode = 0;
  fname_aseg_erode = [];
  if parms.verbose
    fprintf('%s: eroding aseg ROIs...\n',mfilename);
    tic
  end;
  % erode aseg
  fname_aseg_erode = sprintf('%s/aseg_%s.mgz',parms.outdir,parms.erode_outfix);
  try
    fs_erode_aseg(parms.fname_aseg,fname_aseg_erode,...
      'nvoxels',parms.erode_nvoxels,'forceflag',parms.forceflag);
  catch
    fprintf('\n%s: WARNING: aseg erosion failed:\n%s\n\n',...
      mfilename,lasterr);
    errcode = 1;
    return;
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate averages of T1w values for aseg ROIs
function errcode = analyze_aseg_T1w(parms)
  errcode = 0;
  if parms.verbose
    fprintf('%s: T1w aseg analysis for %s...\n',mfilename,parms.fname_T1);
    tic;
  end;
  parms.outstem = [parms.outdir '/T1w'];
  parms.scalefact = parms.T1w_scalefact;
  args = mmil_parms2args(parms,parms.aseg_tags);
  try
    mmil_aseg_analysis(parms.fname_T1,parms.fspath,args{:});
  catch
    fprintf('\n%s: WARNING: T1w aseg analysis failed:\n%s\n\n',...
      mfilename,lasterr);
    errcode = 1;
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate averages of values from proc files for aseg ROIs
function errcode = analyze_aseg_proc(parms)
  errcode = 0;
  for i=1:length(parms.proc_inputlist)
    fname_in = parms.proc_inputlist{i};
    if parms.verbose
      fprintf('%s: proc aseg analysis for %s...\n',mfilename,fname_in);
      tic;
    end;
    parms.outstem = sprintf('%s/%s',parms.outdir,parms.proc_outputlist{i});
    parms.scalefact = parms.proc_scalefacts(i);
    args = mmil_parms2args(parms,parms.aseg_tags);
    try
      mmil_aseg_analysis(fname_in,parms.fspath,args{:});
    catch
      fprintf('\n%s: WARNING: proc aseg analysis failed:\n%s\n\n',...
        mfilename,lasterr);
      errcode = 1;
    end;
    if parms.verbose, toc; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = replace_fuzzy_names(roi_data,fuzzy_names,hemi)
  if isempty(fuzzy_names), return; end;
  for i=1:length(roi_data)
    roi_data(i).roiname = sprintf('ctx-%s-%s',hemi,fuzzy_names{i});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function final_status(errcode,parms)
  if errcode
    if errcode>1
      pstr = 's';
    else
      pstr = [];
    end;
    fprintf('\n%s: WARNING: %d analysis step%s failed\n\n',...
      mfilename,errcode,pstr);
  end;
  if parms.verbose
    fprintf('%s: finished\n',mfilename);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

