function fnames_annot = MMIL_Fibers_to_Surf(varargin)
%function fnames_annot = MMIL_Fibers_to_Surf([options])
%
% Optional Parameters:
%   'subj': FreeSurfer subject name
%     {default = 'fsaverage'}
%   'subjdir': root directory containing subj
%     {default = $FREESURFER_HOME/subjects}
%   'outdir': output directory
%     {default = pwd}
%   'fname_legend': fiber legend csv file
%     if not full path, relative to $MMPS_PARMS/DTI_Fiber directory
%     {default = 'DTI_Fiber_Subdivisions_Legend.csv'}
%   'outfix': output file stem
%     if empty, outfix will be 
%       sprintf('thresh%0.2f_vsm%d_sm',thresh,presmooth,smoothsteps)
%     {default = []}
%   'fibers': vector of fiber numbers to be resampled from atlas
%     {default = [103:110 115:122 133:134]}
%   'subdiv_fibers': vector of fiber subdivision numbers
%     to be resampled from atlas and painted to surface
%     {default = [1031 1041 1051 1061 1071 1081 1091 1101 ...
%                 1151 1152 1161 1162 1173 1175 1183 1185 ...
%                 1193 1195 1203 1205 1213 1215 1223 1225 ...
%                 1331 1332 1333 1341 1342 1343]}
%   'thresh': threshold applied to fiber probability maps
%     {default = 0.1}
%   'presmooth': volume smoothing sigma
%     {default = 2}
%   'projdist': distance (mm) from cortical surface to sample
%     {default = -1}
%   'smoothsteps': smoothing steps on cortical surface
%     {default = 20}
%   'clust_min_val': minimum value for clusters
%     {default = 0.2}
%   'clust_min_area': minimum surface area for clusters
%     {default = 5}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  09/15/12 by Don Hagler
% Last Mod: 02/12/13 Don Hagler
%

%% todo: make get_subdiv_fiber return hemilist,
%%       make necessary changes to loop over hemilist
%%       get hemi/hemilist from DTI_Fiber_Subdivisions_Legend?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames_annot = [];

parms = check_input(varargin);

resample_fibers(parms);

smooth_fibers(parms);

paint_fibers(parms);

find_clusters(parms);

copy_labels(parms);

create_ctab(parms);

fnames_annot = create_annot(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms_filter = {...
    'subj','fsaverage',[],...
    'subjdir',[getenv('FREESURFER_HOME') '/subjects'],[],...
    'outdir',pwd,[],...
    'fname_legend','DTI_Fiber_Subdivisions_Legend.csv',[],...
    'outfix',[],[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % for AtlasTrack
    'fibers',[103:110 115:121 133:134],[],...
    'subdiv_fibers',[1031 1041 1051 1061 1071 1081 1091 1101 ...
                     1151 1152 1161 1162 1173 1175 1183 1185 ...
                     1193 1195 1203 1205 1213 1215 ...
                     1331 1332 1333 1341 1342 1343],[],...
    'locflag',true,[false true],...
    'divide_fibers_flag',true,[false true],...
    'create_paths_flag',false,[false true],...
    'save_mgz_flag',false,[false true],...
  ... % for smooothing fibers
    'fiber_infix','countatlas',[],...
    'thresh',0.1,[0,0.5],...
    'presmooth',2,[0,100],...
  ... % for painting fibers
    'projdist',-1,[-5,5],...
    'smoothsteps',20,[0,100],...
    'outtype','mgz',{'mgh','mgz'},...
  ... % for clustering
    'clust_min_val',0.2,[],...
    'clust_min_area',5,[],...
  ...
    'tags_AtlasTrack',{'fname_T1','fname_FA','fname_V0','M_T1_to_DTI',...
                       'locflag','fibers','subdiv_fibers','divide_fibers_flag',...
                       'combine_fibers_flag','create_paths_flag','save_mgz_flag',...
                       'forceflag','thresh_FA','thresh_prob','min_fiberlen',...
                       'thresh_angle','path_suffix','atlasdir','atlasname',...
                       'fiber_atlasdir','fiber_atlasname','fiber_subdiv_dir',...
                       'outdir','mapsdir','pathsdir'},[],...
    'tags_paint',{'projdist','smoothsteps','subjdir','verbose','forceflag',...
                  'hemi','outstem','outtype'},[],...
    'tags_annot',{'indir','outdir','source_subj','fname_ctab'...
                  'subj','subjdir','annotname','hemilist'...
                  'verbose','forceflag'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);

  % get name of T1 file
  parms.FSPath = [parms.subjdir '/' parms.subj];
  if ~exist(parms.FSPath,'dir')
    error('FreeSurfer recon dir %s not found',parms.FSPath);
  end;  
  parms.fname_T1 = [parms.FSPath '/mri/nu.mgz'];
  if ~exist(parms.fname_T1,'file')
    parms.fname_T1 = [parms.FSPath '/mri/orig.mgz'];
  end;
  if ~exist(parms.fname_T1,'file')
    error('file %s not found',parms.fname_T1);
  end;

  % read fiber legend to get fiber numbers and names
  if mmil_isrelative(parms.fname_legend)
    dir_parms = getenv('MMPS_PARMS');
    parms.fname_legend = [dir_parms '/DTI_Fiber/' parms.fname_legend];
  end;
  FiberLegend = DTI_MMIL_Read_Fiber_Legend(parms.fname_legend);
  parms.fiber_numbers = [FiberLegend.FiberNumber];
  parms.fiber_names = regexprep({FiberLegend.FiberName},' ','_');

  % set output names
  parms.AtlasTrack_dir = [parms.outdir '/resampled_fibers'];
  parms.smoothed_dir = [parms.outdir '/smoothed_fibers'];
  parms.painted_dir = [parms.outdir '/painted_fibers'];
  parms.cluster_dir = [parms.outdir '/cluster_labels'];
  parms.infix = sprintf('thresh%0.3f_vsm%d',parms.thresh,parms.presmooth);
  if parms.smoothsteps
    parms.infix = sprintf('%s-sm%d',parms.infix,parms.smoothsteps);
  end;
  parms.clust_infix = sprintf('%s_cminv%0.2f_cmina%0.1f',...
      parms.infix,parms.clust_min_val,parms.clust_min_area);
  if isempty(parms.outfix)
    parms.outfix = sprintf('thresh%0.3f_vsm%d_sm%d_cminv%0.2f_cmina%0.1f',...
      parms.thresh,parms.presmooth,parms.smoothsteps,...
      parms.clust_min_val,parms.clust_min_area);
  end;
  parms.label_dir = sprintf('%s/labels_%s',parms.outdir,parms.outfix);
  parms.fname_ctab = sprintf('%s/colors_%s.ctab',parms.outdir,parms.outfix);
  parms.annotname = ['fparc_' parms.outfix];

  % check fiber numbers
  bad_fibers = setdiff(10*parms.fibers,parms.fiber_numbers);
  bad_subdiv_fibers = setdiff(parms.subdiv_fibers,parms.fiber_numbers);  
  if ~isempty(bad_fibers)
    error('invalid fiber numbers specified: %s',...
      sprintf('%d ',bad_fibers));
  end;
  if ~isempty(bad_subdiv_fibers)
    error('invalid subdiv fiber numbers specified: %s',...
      sprintf('%d ',bad_subdiv_fibers));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resample_fibers(parms)
  parms.outdir = parms.AtlasTrack_dir;
  args = mmil_parms2args(parms,parms.tags_AtlasTrack);
  dti_AtlasTrack(args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function smooth_fibers(parms)
  mmil_mkdir(parms.smoothed_dir);
  for s=1:length(parms.subdiv_fibers)
    [fiber_name,hemi,fnum] = get_subdiv_fiber(parms,s);
    fname_in = sprintf('%s/fiber_maps/fiber_%d_%s.mat',...
      parms.AtlasTrack_dir,fnum,parms.fiber_infix);
    fname_out = sprintf('%s/%s_thresh%0.3f_vsm%d.mgz',...
      parms.smoothed_dir,fiber_name,parms.thresh,parms.presmooth);    
    if ~exist(fname_out,'file') || parms.forceflag
      [vol,M] = mmil_load_sparse(fname_in);
      vol = 1.0 .* (vol > parms.thresh);
      if parms.presmooth>0
        vol = mmil_smooth3d(vol,parms.presmooth);
        vol = 1.0 * (vol > parms.thresh);
      end;
      fs_save_mgh(vol,fname_out,M);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function paint_fibers(parms)
  mmil_mkdir(parms.painted_dir);
  for s=1:length(parms.subdiv_fibers)
    [fiber_name,parms.hemi] = get_subdiv_fiber(parms,s);
    fname_in = sprintf('%s/%s_thresh%0.3f_vsm%d.mgz',...
      parms.smoothed_dir,fiber_name,parms.thresh,parms.presmooth);
    parms.outstem = sprintf('%s/%s_thresh%0.3f_vsm%d',...
      parms.painted_dir,fiber_name,parms.thresh,parms.presmooth);
    args = mmil_parms2args(parms,parms.tags_paint);
    fs_paint(parms.subj,fname_in,args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function find_clusters(parms)
  mmil_mkdir(parms.cluster_dir);
  for s=1:length(parms.subdiv_fibers)
    [fiber_name,hemi] = get_subdiv_fiber(parms,s);
    fname_in = sprintf('%s/%s_%s-%s.%s',...
      parms.painted_dir,fiber_name,parms.infix,hemi,parms.outtype);
    outstem = sprintf('%s/%s.%s_%s',...
      parms.cluster_dir,hemi,fiber_name,parms.clust_infix);
    fname_sum = sprintf('%s_summary.txt',outstem);
    if ~exist(fname_sum,'file') || parms.forceflag
      % run mri_surfcluster
      cmd = 'mri_surfcluster';
      cmd = sprintf('%s --in %s',cmd,fname_in);
      cmd = sprintf('%s --thmin %0.6f',cmd,parms.clust_min_val);
      cmd = sprintf('%s --subject %s --hemi %s',cmd,parms.subj,hemi);
      cmd = sprintf('%s --sd %s',cmd,parms.subjdir);
      cmd = sprintf('%s --olab %s',cmd,outstem);
      cmd = sprintf('%s --minarea %0.6f',cmd,parms.clust_min_area);
      cmd = sprintf('%s --sum %s',cmd,fname_sum);
      if parms.verbose
        fprintf('%s: cmd = %s\n',mfilename,cmd);
      end;
      [status,result] = unix(cmd);
      if status
        error('cmd %s failed:\n%s',cmd,result);
      end;
      if parms.verbose
        fprintf('%s\n',result);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_labels(parms)
  mmil_mkdir(parms.label_dir);
  for s=1:length(parms.subdiv_fibers)
    [fiber_name,hemi] = get_subdiv_fiber(parms,s);
    % find largest cluster
    flist = dir(sprintf('%s/%s.%s_%s-*.label',...
      parms.cluster_dir,hemi,fiber_name,parms.clust_infix));
    if isempty(flist), continue; end;
    max_f = 1;
    max_nverts = 0;
    for f=1:length(flist)
      v = fs_read_label([parms.cluster_dir '/' flist(f).name]);
      if length(v) > max_nverts
        max_f = f;
        max_nverts = length(v);
      end;
    end;
    fname_in = [parms.cluster_dir '/' flist(max_f).name];
    % get fiber name without 'R_' or 'L_' at beginning
    n = regexp(fiber_name,'^\w_(?<name>\w+)','names');
    if isempty(n)
      name = fiber_name;
    else
      name = n.name;
    end;
    fname_out = sprintf('%s/%s.%s.label',...
      parms.label_dir,hemi,name);
    if exist(fname_in,'file')
      if ~exist(fname_out,'file') || parms.forceflag
        cmd = sprintf('cp %s %s',fname_in,fname_out);
        if parms.verbose
          fprintf('%s: %s\n',mfilename,cmd);
        end;
        [s,r] = unix(cmd);
        if s, error('cmd %s failed:\n%s',cmd,r); end;      
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_ctab(parms)
  if ~exist(parms.fname_ctab,'file') || parms.forceflag
    namelist = [];
    for s=1:length(parms.subdiv_fibers)
      [fiber_name,hemi,fnum] = get_subdiv_fiber(parms,s);
      n = regexp(fiber_name,'\w_(?<name>\w+)','names');
      if isempty(n)
        name = fiber_name;
      else
        name = n.name;
      end;
      fname_in = sprintf('%s/%s.%s.label',parms.label_dir,hemi,name);
      if exist(fname_in,'file'), namelist{end+1} = name; end;
    end;
    if isempty(namelist)
      error('no labels for annotation');
    end;
    namelist = unique(namelist);
    fs_write_ctab(namelist,parms.fname_ctab);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnames_annot = create_annot(parms)
  parms.source_subj = parms.subj;
  args = mmil_parms2args(parms,parms.tags_annot);
  [fnames_annot,fnames_label] = fs_labels2annot(parms.label_dir,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fiber_name,hemi,fnum] = get_subdiv_fiber(parms,s)
  fnum = parms.subdiv_fibers(s);
  ind = find(fnum==parms.fiber_numbers);
  fiber_name = parms.fiber_names{ind};
  if ismember(fnum,[1211,1213,1216,...
                    1221,1223,1226,...
                    1231,1233,1235,...
                    1241,1243,1245,1247])
    hemi = 'rh';
  elseif ismember(fnum,[1212,1215,1217,...
                        1222,1225,1227,...
                        1232,1234,1236,...
                        1242,1244,1246,1248])
    hemi = 'lh';
  elseif mod(floor(fnum/10),2) % is odd, therefore right hemisphere
    hemi = 'rh';
  else
    hemi = 'lh';
  end;
return;


