function mmil_relabel_rois(fname_in,fname_info,varargin)
%function mmil_relabel_rois(fname_in,fname_info,[options])
%
% Purpose: relabel surface-based ROIs according
%  to degree of correlation with aparc ROIs
%
% Required Input:
%   fname_in: input file name
%     should be multi-frame mgh/mgz file containing surface weighted ROIs
%     expected to be on icosahedral surface (e.g. fsaverage)
%  fname_info: file name of csv spreadsheet file containing
%     information for each ROI in fname_in
%     first row should contain column headers
%     must have same number of additional rows as there are frames in fname_in
%     {default = []}
%
% Optional Parmaeters:
%  'fname_out_info' : output file name for ROI info
%     if ommitted, fname_out_info will be constructed from fname_info
%     may be full path or relative to outdir
%     {default = []}
%  'fname_out_labels' : output file name for ROI labels
%     if ommitted, fname_out_labels will be constructed from fname_info
%     may be full path or relative to outdir
%     {default = []}
%  'outdir' : directory for output file
%     if ommitted, will use same directory as input file
%     {default = []}
%  'outfix': suffix appended to output file name
%     ignored if fname_out is supplied
%     {default = 'relabeled'}
%  'hemi': cortical hemisphere (lh or rh)
%    necessary only if input file name does not have hemi tag at end
%      e.g. 'stem-lh.mgh' or 'stem-rh.mgh'
%    if file name does have hemi tag, this parameter will be ignored
%    {default = []}
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
% Created:  03/05/15 by Don Hagler
% Last Mod: 03/05/15 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = check_input(fname_in,fname_info,varargin);

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(parms.fname_out_info,'file') ||...
   ~exist(parms.fname_out_labels,'file') || parms.forceflag
  % load ROI info
  roi_info = load_roi_info(parms);

  % label each ROI according to correlation with aparc ROIs
  roi_labels = label_rois(roi_info,parms);

  % save info on ROIs
  save_info(roi_info,roi_labels,parms);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,fname_info,options)
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
    'fname_info',fname_info,[],...
  ...
    'fname_out_info',[],[],...
    'fname_out_labels',[],[],...
    'outdir',[],[],...            
    'outfix','relabeled',[],...
    'hemi',[],[{'lh','rh'}],...
    'fstem_aparc','aparc',[],...
    'subjname','fsaverage',[],...
    'subjdir',[getenv('FREESURFER_HOME') '/subjects'],[],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  });

  % check fname_in
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;

  % check fname_info
  if ~exist(parms.fname_info,'file')
    error('file %s not found',parms.fname_info);
  end;

  % extract the file stem for input info file
  [indir,inname,inext] = fileparts(parms.fname_info);

  % get hemisphere
  pat = sprintf('(?<instem>.+)-(?<hemi>[lr]h)-info');
  n = regexp(inname,pat,'names');
  if isempty(n)
    if isempty(parms.hemi)
      error('input file name %s does not match expected pattern and hemi parameter not specified',inname);
    else
      fprintf('%s: WARNING: input file name %s does not match expected pattern',inname);
    end;
    instem = inname;
  else
    instem = n.instem;
    parms.hemi = n.hemi;
  end;

  % set output dir
  if isempty(parms.outdir)
    parms.outdir = indir;
  end;

  % construct the output file names
  if isempty(parms.fname_out_info)
    parms.fname_out_info = sprintf('%s-%s-%s-info.csv',...
      instem,parms.outfix,parms.hemi);
  end;
  if mmil_isrelative(parms.fname_out_info)
    parms.fname_out_info = fullfile(parms.outdir,parms.fname_out_info);
  end;
  [parms.outdir,outstem] = fileparts(parms.fname_out_info);
  
  if isempty(parms.fname_out_labels)
    parms.fname_out_labels = regexprep(parms.fname_out_info,...
      '-info\.csv','-labels\.csv');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_info = load_roi_info(parms)
  roi_info = [];
  if parms.verbose
    fprintf('%s: loading ROI info from %s...\n',mfilename,parms.fname_info);
  end;
  roi_info = mmil_csv2struct(parms.fname_info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_labels = label_rois(roi_info,parms)
  if parms.verbose
    fprintf('%s: labeling ROIs...\n',mfilename);
  end;
  % load ROIs
  roivals = squeeze(fs_load_mgh(parms.fname_in));
  [num_verts,num_rois] = size(roivals);
  if length(roi_info)~=num_rois
    error('number of ROIs (%d) do not match number of entries in ROI info (%d)',...
      num_rois,length(roi_info));
  end;
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
  label_data = cat(1,{'roi','label'},label_data);
  mmil_write_csv(parms.fname_out_labels,label_data);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_info(roi_info,roi_labels,parms)
  if parms.verbose
    fprintf('%s: saving ROI information to %s...\n',...
      mfilename,parms.fname_out_info);
  end;
  % get info for ROIs and add labels
  roi_out_info = roi_info;
  for r=1:length(roi_info)
    roi_out_info(r).label = roi_labels{r};
  end;
  % reorder fields to put label first
  nfields = length(fieldnames(roi_out_info));
  roi_out_info = orderfields(roi_out_info,[nfields,1:nfields-1]);
  % write to csv file
  mmil_struct2csv(roi_out_info,parms.fname_out_info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


