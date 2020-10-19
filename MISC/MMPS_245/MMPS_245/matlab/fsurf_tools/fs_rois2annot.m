function fs_rois2annot(fnames_roi,fnames_info,varargin)
%function fs_rois2annot(fnames_roi,fnames_info,[options])
%
% Purpose; convert multi-frame ROI files to annot parcellation file
%
% Required parameters:
%   fnames_roi: full path of multi-frame ROI file in mgh/mgz format
%     may be cell array containing files for each hemisphere
%   fnames_info: full path of spreadsheet with ROI names in csv format
%               should contain "label" column unless specified by roi_label
%     may be cell array containing files for each hemisphere
%       must match fnames_roi
%
% Optional paramers:
%   'hemilist': cell array of hemispheres
%     must match fnames_roi
%     {default = {'lh','rh'}}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'cparc'}
%   'roi_label': column label used for ROI names
%     {default = 'label'}
%   'subj': subject name used in creating label and annotation files
%     {default = 'fsaveage'}
%   'subjdir': full path of directory containing FreeSurfer recons
%     if empty, will use $FREESURFER_HOME/subjects
%     {default = []}
%   'thresh': threshold applied to ROI values
%     {default = 0.5}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  09/11/15 by Don Hagler
% Last Mod: 04/01/16 by Don Hagler
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin,{...
  'hemilist',{'lh','rh'},[],...
  'outdir',pwd,[],...
  'outstem','cparc',[],...
  'roi_label','label',[],...
  'subj','fsaverage',[],...
  'subjdir',[getenv('FREESURFER_HOME') '/subjects'],[],...
  'thresh',0.5,[],...
  'forceflag',false,[false true],...
});

%% todo: check that fnames_roi, fnames_info, hemilist are cell arrays
%%       with matching lengths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labeldir = [parms.outdir '/labels'];
mmil_mkdir(labeldir);

% compile ROI names
roi_names = [];
nhemi = length(parms.hemilist);
for h=1:nhemi
  hemi = parms.hemilist{h};
  roi_info = mmil_csv2struct(fnames_info{h});
  roi_names = cat(2,roi_names,{roi_info.(parms.roi_label)});
end;
roi_names = unique(roi_names);
roi_names = cat(2,'unlabeled',roi_names);

% create color table file
fname_ctab = sprintf('%s/%s.ctab',parms.outdir,parms.outstem);
if ~exist(fname_ctab,'file') || parms.forceflag
  fs_write_ctab(roi_names,fname_ctab);
end;

nhemi = length(parms.hemilist);
for h=1:nhemi
  hemi = parms.hemilist{h};
  roi_vals = squeeze(fs_load_mgh(fnames_roi{h}));
  roi_info = mmil_csv2struct(fnames_info{h});
  nroi = length(roi_info);
  if size(roi_vals,2) ~= nroi
    error('mismatch in number of ROIs in %s and %s',...
      fname_roi,fname_info);
  end;

  % find ROI with greatest value for each vertex
  [max_vals,i_max] = max(roi_vals,[],2);

  % find vertices with max values greater than threshold
  i_thresh = find(max_vals>parms.thresh);
  i_unlabeled = find(max_vals<=parms.thresh);

  % for each ROI, identify the vertices
  for r=1:nroi
    roi_name = roi_info(r).label;
    % write label file
    fname_label = sprintf('%s/%s.%s.label',...
      labeldir,hemi,roi_name);
    if ~exist(fname_label,'file') || parms.forceflag
      % find vertices in this ROI
      v = find(i_max == r);
      % exclude those below threshold
      v = intersect(v,i_thresh);
      fs_write_label(v,fname_label,parms.subj);
    end;
  end;
  
  % create an extra ROI for unlabeled vertices
  if ~isempty(i_unlabeled)
    roi_name = 'unlabeled';
    fname_label = sprintf('%s/%s.%s.label',...
      labeldir,hemi,roi_name);
    if ~exist(fname_label,'file') || parms.forceflag
      fs_write_label(i_unlabeled,fname_label,parms.subj);
    end;
  end;
    
  % convert labels to annotation
  fs_labels2annot(labeldir,...
    'fname_ctab',fname_ctab,'outdir',parms.outdir,...
    'source_subj',parms.subj,'subj',parms.subj,...
    'subjdir',parms.subjdir,'annotname',parms.outstem,...
    'hemilist',{hemi},...
    'forceflag',parms.forceflag);
end;

