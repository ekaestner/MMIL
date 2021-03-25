function mmil_group_rois(fname_in,varargin)
%function mmil_group_rois(fname_in,[options])
%
% Purpose: group surface-based ROIs based on similarity of spatial patterns
%   either using kmeans or singular value decomposition (SVD)
%
% Required Input:
%   fname_in: input file name
%     should be multi-frame mgh/mgz file containing surface weighted ROIs
%     expected to be on icosahedral surface (e.g. fsaverage)
%
% Optional Parameters:
%  'kmeans_flag': [0|1] use kmeans, otherwise use svd
%     if kmeans_flag = 1, find ROIs that contribute strongly to
%       only one group identified using kmeans
%     if kmeans_flag = 0, iteratively find the ROIs contributing
%       strongly to the first SVD component, drawing without replacement
%    {default = 0}
%  'kmeans_replicates': number of times to repeat kmeans clustering
%     to avoid local minima
%     {default = 5}
%  'fname_info': file name of csv spreadsheet file containing
%     information for each ROI in fname_in
%     first row should contain column headers
%     must have same number of additional rows as there are frames in fname_in
%     if supplied, info will be passed to output info file
%     {default = []}
%  'fname_out' : output file name
%     if ommitted, fname_out will be constructed from fname_in
%     may be full path or relative to outdir
%     {default = []}
%  'thresh': grouping threshold; ranges from 0 to 1
%     lower values (e.g. 0.5) create larger ROIs
%     higher values (e.g. 0.9) create smaller ROIs
%     {default = 0.8}
%  'ngroups': number of groups of ROIs to be identified
%     depending on input ROIs and threshold, fewer groups may be identified
%     {default = 10}
%  'outdir' : directory for output file
%     if ommitted, will use same directory as input file
%     {default = []}
%  'outfix': suffix appended to output file name
%     ignored if fname_out is supplied
%     {default = 'groups'}
%  'hemi': cortical hemisphere (lh or rh)
%    necessary only if input file name does not have hemi tag at end
%      e.g. 'stem-lh.mgh' or 'stem-rh.mgh'
%    if file name does have hemi tag, this parameter will be ignored
%    {default = []}
%  'ico': [3:6] icosahedral order used for downsampling ROIs
%    will use first 10242 vertices for SVD calculations
%    {default = 5}
%  'corrflag': [0|1] perform SVD on between ROI correlation matrix
%    instead of on ROI surface values
%    ignored if kmeans_flag = 1
%    {default = 0}
%  'aparc_label_flag': [0|1] assign ROI names for each group according
%     to degree of correlation with aparc ROIs
%    {default = 1}
%  'fstem_aparc': file stem for aparc ROIs
%    {default = 'aparc'}
%  'subjname': FreeSurfer subject name used to load aparc ROIs
%    {default = 'fsaverage'}
%  'subjdir': directory containing FreeSurfer subject
%    {default = [getenv('FREESURFER_HOME') '/subjects']}
%  'verbose': [0|1] display mri_surf2surf output
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%
% Created:  08/20/14 by Don Hagler
% Last Mod: 03/05/15 by Don Hagler
%

roi_groups = [];
if ~mmil_check_nargs(nargin,1), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(fname_in,varargin);

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(parms.fname_out,'file') ||...
   ~exist(parms.fname_out_info,'file') ||...
   (parms.aparc_label_flag && ~exist(parms.fname_out_labels,'file')) ||...
    parms.forceflag
  % load ROI file and get downsampled ico vals
  [vals,icovals,parms] = load_rois(parms);

  % load ROI info
  roi_info = load_roi_info(parms);

  % find ROI groups
  if parms.kmeans_flag
    [roi_groups,parms] = find_groups_kmeans(icovals,parms);
  else
    [roi_groups,parms] = find_groups_svd(icovals,parms);
  end;

  % calculate averages for each group
  parms = calculate_averages(vals,roi_groups,parms);

  % label each group according to correlation with aparc ROIs
  if parms.aparc_label_flag
    roi_labels = label_groups(roi_info,parms);
  else
    roi_labels = [];
  end;

  % save info on ROI groups
  save_info(roi_groups,roi_info,roi_labels,parms);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
  ...
    'kmeans_flag',false,[false true],...
    'kmeans_replicates',5,[1,1000],...
    'fname_info',[],[],...
    'fname_out',[],[],...
    'thresh',0.8,[],...
    'ngroups',10,[],...
    'outdir',[],[],...            
    'outfix','groups',[],...
    'hemi',[],[{'lh','rh'}],...
    'ico',5,[3:6],...
    'corrflag',false,[false true],...
    'aparc_label_flag',true,[false true],...
    'fstem_aparc','aparc',[],...
    'subjname','fsaverage',[],...
    'subjdir',[getenv('FREESURFER_HOME') '/subjects'],[],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % undocumented
    'intype','mgh',{'mgh','mgz','w'},...
    'outtype',[],{'mgh','mgz','w'},...
  });

  Ns = [42 162 642 2562 10242 40962 163842];
  parms.nicoverts = Ns(parms.ico);

  % check fname_in
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;

  % check fname_info
  if ~isempty(parms.fname_info) && ~exist(parms.fname_info,'file')
    error('file %s not found',parms.fname_info);
  end;

  % extract the file stem
  [indir,inname,inext] = fileparts(parms.fname_in);
  if strcmp(parms.intype,'mgh') && strcmp(inext,'.mgz')
    parms.intype = 'mgz';
  end;
  if ~strcmp(['.' parms.intype],inext)
    error('fname_in %s missing correct file extension for intype %s',...
      parms.fname_in,parms.intype);
  end;
  if isempty(parms.outtype)
    parms.outtype = parms.intype;
    parms.outext = inext;
  else
    parms.outext = ['.' parms.outtype];
  end;

  pat = sprintf('(?<instem>.+)-(?<hemi>[lr]h)');
  n = regexp(inname,pat,'names');
  if isempty(n)
    if isempty(parms.hemi)
      error('input file name %s missing hemi tag (lh or rh) at end of name and hemi parameter not specified',inname);
    else
      fprintf('%s: WARNING: input file name %s missing hemi tag (lh or rh) at end of name',inname);
    end;
    instem = inname;
  else
    instem = n.instem;
    parms.hemi = n.hemi;
  end;
  if isempty(parms.outdir)
    parms.outdir = indir;
  end;

  % construct the output file name
  if isempty(parms.fname_out)
    outstem = sprintf('%s-%s',instem,parms.outfix);
    parms.fname_out = sprintf('%s-%s%s',outstem,parms.hemi,parms.outext);
  end;
  if mmil_isrelative(parms.fname_out)
    parms.fname_out = fullfile(parms.outdir,parms.fname_out);
  end;
  [parms.outdir,outstem] = fileparts(parms.fname_out);
  parms.fname_out_info = sprintf('%s/%s-info.csv',parms.outdir,outstem);
  parms.fname_out_labels = sprintf('%s/%s-labels.csv',parms.outdir,outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vals,icovals,parms] = load_rois(parms)
  if parms.verbose
    fprintf('%s: loading ROIs from %s...\n',mfilename,parms.fname_in);
  end;
  vals = fs_load_mgh(parms.fname_in);
  vals = squeeze(vals);
  [parms.nverts,parms.nrois] = size(vals);
  icovals = vals(1:parms.nicoverts,:);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_info = load_roi_info(parms)
  roi_info = [];
  if ~isempty(parms.fname_info)
    if parms.verbose
      fprintf('%s: loading ROI info from %s...\n',mfilename,parms.fname_info);
    end;
    roi_info = mmil_csv2struct(parms.fname_info);
    if length(roi_info)~=parms.nrois
      error('number of ROIs (%d) do not match number of entries in ROI info (%d)',...
        parms.nrois,length(roi_info));
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [roi_groups,parms] = find_groups_kmeans(icovals,parms)
  roi_groups = cell(parms.ngroups,1);
  if parms.verbose
    fprintf('%s: finding groups using kmeans...\n',mfilename);
    tic;
  end;
  X = icovals;
  idx = kmeans(X,parms.ngroups+1,...
    'EmptyAction','drop','Replicates',parms.kmeans_replicates);
  % exclude most frequenly represented group (no ROI)
  ind_groups = setdiff(1:parms.ngroups+1,mode(idx));
  % create mask of each group
  masks = zeros(parms.nicoverts,parms.ngroups);
  for r=1:parms.ngroups
    ind_mask = find(idx == ind_groups(r));
    masks(ind_mask,r) = 1;
  end;
  % identify seed points that strongly contribute to each group
  group_vals = icovals'*masks;
  maxvals = max(group_vals);
  norm_group_vals = bsxfun(@rdivide,group_vals,maxvals);
  group_mask = 1.0*(norm_group_vals > parms.thresh);
  % exclude ROIs if they contribute to more than one group
  ngroups_per_seed = sum(group_mask,2);
  ind_zero = find(ngroups_per_seed>1);
  group_mask(ind_zero,:) = 0;
  % compile ROI indices for each group
  for r=1:parms.ngroups
    ind_rois = find(group_mask(:,r))';
    if isempty(ind_rois), continue; end;
    roi_groups{r} = ind_rois;
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [roi_groups,parms] = find_groups_svd(icovals,parms)
  roi_groups = cell(parms.ngroups,1);
  if parms.verbose
    fprintf('%s: finding groups using SVD...\n',mfilename);
    tic;
  end;
  if parms.corrflag
    X = corr(icovals);
  else
    X = icovals;
  end;
  ind_rem = [1:parms.nrois];
  for r=1:parms.ngroups
    % find ROIs for first component of all remaining ROIs
    [X0,ind_rois0] = find_first_comp(X,ind_rem,parms);
    % find ROIs for first component of data for just those ROIs
    [X1,ind_rois1] = find_first_comp(X0,ind_rois0,parms);
    % select ROIs for this group
    ind_rois = ind_rem(ind_rois0(ind_rois1));
    roi_groups{r} = ind_rois;
    % exclude selected ROIs from remaining
    ind_rem = setdiff(ind_rem,ind_rois);
  end;
  if parms.verbose, toc; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X1,ind_rois] = find_first_comp(X0,ind_rem,parms)
  if parms.corrflag
    X1 = X0(ind_rem,ind_rem);
  else
    X1 = X0(:,ind_rem);
  end;
  [U,S,V] = svd(X1,0);
  % find ROIs for first component
  Sx = zeros(size(S));
  Sx(1,1) = S(1,1);
  Y = U*Sx*V';
  maxy = max(Y,[],1);
  maxy = maxy / max(maxy);
  ind_rois = find(maxy > parms.thresh);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = calculate_averages(vals,roi_groups,parms)
  if parms.verbose
    fprintf('%s: calculating ROI group averages and saving to %s...\n',...
      mfilename,parms.fname_out);
  end;
  ind_keep = find(~cellfun(@isempty,roi_groups));
  roi_groups = roi_groups(ind_keep);
  parms.ngroups = length(roi_groups);
  vals_groups = zeros(parms.nverts,parms.ngroups);
  for r=1:parms.ngroups
    ind_rois = roi_groups{r};
    % calculate mean across selected ROIs
    mean_vals = mean(vals(:,ind_rois),2);
    % normalize so maximum is 1
    vals_groups(:,r) = mean_vals / max(mean_vals);
  end;
  vals_groups = reshape(vals_groups,[parms.nverts,1,1,parms.ngroups]);
  fs_save_mgh(vals_groups,parms.fname_out);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_labels = label_groups(roi_info,parms)
  if parms.verbose
    fprintf('%s: labeling ROI groups...\n',mfilename);
  end;
  % load averaged group ROIs
  roivals = squeeze(fs_load_mgh(parms.fname_out));
  [num_verts,num_rois] = size(roivals);
  % load aparc ROIs
  fname_aparc = sprintf('%s/%s/label/%s.aparc.annot',...
    parms.subjdir,parms.subjname,parms.hemi);
  [aparc_roinums,aparc_roilabels] = fs_read_annotation(fname_aparc);
  num_aparc_rois = length(aparc_roilabels);
  if length(aparc_roinums) ~= num_verts
    error('number of vertices in aparc does not match\n',mfilename);
  end;
  % create matrix of aparc ROIs
  aparc_roivals = zeros(num_verts,num_aparc_rois);
  for i=1:num_aparc_rois
    ind_roi = find(aparc_roinums==i);
    aparc_roivals(ind_roi,i) = 1;
  end;
  % calculate correlation between roivals and aparc_roivals
  R = corr(roivals,aparc_roivals);
  % find maximum correlation values for each ROI
  [maxR,ind_max] = max(R,[],2);
  % create labels for each ROI
  roi_labels = cell(num_rois,1);
  for i=1:num_rois
    roi_labels{i} = aparc_roilabels{ind_max(i)};
  end;
  % check for duplicates, add numbers to distinguish them
  uind = unique(ind_max);
  if length(uind)<length(ind_max)
    for i=1:length(uind)
      ind = uind(i);
      ind_rep = find(ind_max == ind);
      nreps = length(ind_rep);
      if nreps>1
        repR = maxR(ind_rep);
        [repR_sorted,ind_sort] = sort(repR,'descend');
        for j=1:nreps
          tmp_label = roi_labels{ind_rep(ind_sort(j))};
          roi_labels{ind_rep(ind_sort(j))} = sprintf('%s-%d',tmp_label,j);
        end;  
      end;
    end;
  end;
  % write to csv file
  label_data = cat(2,num2cell(1:num_rois)',roi_labels);
  label_data = cat(1,{'group','label'},label_data);
  mmil_write_csv(parms.fname_out_labels,label_data);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_info(roi_groups,roi_info,roi_labels,parms)
  if parms.verbose
    fprintf('%s: saving ROI group information to %s...\n',...
      mfilename,parms.fname_out_info);
  end;
  % get info for ROIs in groups and add group numbers
  ind_all = [roi_groups{:}];
  roi_out_info = roi_info(ind_all);
  for r=1:length(roi_groups)
    ind_rois = roi_groups{r};
    [tmp,ind] = intersect(ind_all,ind_rois);
    for j=1:length(ind)
      roi_out_info(ind(j)).group = r;
      if ~isempty(roi_labels)
        roi_out_info(ind(j)).label = roi_labels{r};
      end;
    end;
  end;
  % reorder fields to put group first
  nfields = length(fieldnames(roi_out_info));
  if isempty(roi_labels)
    roi_out_info = orderfields(roi_out_info,[nfields,1:nfields-1]);
  else
    roi_out_info = orderfields(roi_out_info,[nfields-1,nfields,1:nfields-2]);
  end;
  % write to csv file
  mmil_struct2csv(roi_out_info,parms.fname_out_info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


