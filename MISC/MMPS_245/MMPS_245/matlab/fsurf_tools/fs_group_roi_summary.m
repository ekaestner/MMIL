function fs_group_roi_summary(fname_out,varargin);
%function fs_group_roi_summary(fname_out,[options]);
%
% Usage:
%  fs_group_roi_summary(fname_out,'key1', value1,...);
%
% Required input:
%  fname_out: output file name
%    This function will generate a csv (comma separated value) format
%    spreadsheet containing morphological data from Freesurfer aseg
%    and aparc ROIs for a given subjects directory
%
% Optional parameters:
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    {default = $SUBJECTS_DIR}
%  'subjlist': cell array of subject names
%    output will restricted only to those subjects (must exist in subjdir)
%    if empty or ommitted, all subjects in subjdir will be included
%    {default: []}
%  'grouplist': cell array of group label for each subject
%    this list must be the same length as subjlist
%    this information will simply be included as an extra column
%    {default: []}
%  'a2005_flag': [0|1] toggle use of a2005s.stats files
%    (otherwise use aparc.stats)
%    {default: 0}
%  'checkstatus_flag': [0|1] toggle check that recon is complete
%    {default: 0}
%  'FS_version': FreeSurfer version (important for checkstatus)
%    {default: 450}
%  'forceflag': [0|1] overwrite existing output file
%    {default: 1}
%
% Parameters controlling what information is included in output file
%  'aseg': [0|1] toggle whether to include aseg volumes
%    {default: 1}
%  'aparc': [0|1] toggle whether to include aparc measures
%    if this value is 0, no aparc measures will be included and the
%      following 8 parameters will be ignored
%    {default: 1}
%  'grayvol': [0|1] toggle whether to include gray matter volume
%    {default: 1}
%  'surfarea': [0|1] toggle whether to include surface area
%    {default: 1}
%  'thickavg': [0|1] toggle whether to include average thickness
%    {default: 1}
%  'thickstd': [0|1] toggle whether to include standard deviation of thickness
%    {default: 1}
%  'meancurv': [0|1] toggle whether to include mean curvature
%    {default: 1}
%  'gausscurv': [0|1] toggle whether to include gaussian curvature
%    {default: 1}
%  'foldind': [0|1] toggle whether to include folding index
%    {default: 1}
%  'curvind': [0|1] toggle whether to include curvature index
%    {default: 1}
%
% created:  07/13/07 by Don Hagler
% last mod: 11/02/11 by Don Hagler
%

hemilist = {'Left','Right'};
aseg_roilist = [1:28,40:60,77:79,20001:20003];
exclude_aseg_roilist = [1,6,9,19:23,27,40,45,48,55,56,59];
aseg_roilist = setdiff(aseg_roilist,exclude_aseg_roilist);

num_aseg_rois = length(aseg_roilist);
num_aparc_rois = 36;
num_aparc_rois_2005 = 79;

%% todo: accept cell matrix of extra info

%% todo: keep example aseg and aparc stats in example matfile for writing header
%% todo: accept spreadsheet with subject info
%% todo: accept list of aseg roi numbers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'subjdir',[],[],...
  'subjlist',[],[],...
  'grouplist',[],[],... 
  'a2005_flag',false,[false true],...
  'checkstatus_flag',false,[false true],...
  'FS_version',450,[],...
  'forceflag',false,[false true],...
  'aseg',true,[false true],...
  'aparc',true,[false true],...
  'grayvol',true,[false true],...
  'surfarea',true,[false true],... 
  'thickavg',true,[false true],...
  'thickstd',true,[false true],...
  'meancurv',true,[false true],...
  'gausscurv',true,[false true],...
  'foldind',true,[false true],...
  'curvind',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

if isempty(parms.subjdir)
  parms.subjdir = getenv('SUBJECTS_DIR');
  if isempty(parms.subjdir)
    fprintf('%s: ERROR: SUBJECTS_DIR not defined as an environment variable\n',mfilename);
    return;
  end;
end;
if ~exist(parms.subjdir,'dir')
  fprintf('%s: ERROR: subjects dir %s not found\n',mfilename,parms.subjdir);
  return;
end;

if exist(fname_out,'file') & ~parms.forceflag
  fprintf('%s: output file %s already exists\n',...
    mfilename,fname_out);
  fprintf('      delete first or set forceflag = 1\n');
  return;
end;

aparc_meas_list = {};
if parms.aparc
  if parms.grayvol, aparc_meas_list{end+1} = 'grayvol'; end;
  if parms.surfarea, aparc_meas_list{end+1} = 'surfarea'; end;
  if parms.thickavg, aparc_meas_list{end+1} = 'thickavg'; end;
  if parms.thickstd, aparc_meas_list{end+1} = 'thickstd'; end;
  if parms.meancurv, aparc_meas_list{end+1} = 'meancurv'; end;
  if parms.gausscurv, aparc_meas_list{end+1} = 'gausscurv'; end;
  if parms.foldind, aparc_meas_list{end+1} = 'foldind'; end;
  if parms.curvind, aparc_meas_list{end+1} = 'curvind'; end;
  if isempty(aparc_meas_list), parms.aparc = 0; end;
end;

if isempty(parms.aseg) & isempty(parms.aparc)
  error('nothing selected for output');
end;

if ~isempty(parms.grouplist) & length(parms.grouplist) ~= length(parms.subjlist)
  error('length of grouplist must match subjlist');
end;

if parms.a2005_flag
  num_aparc_rois = num_aparc_rois_2005;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(fname_out,'wt');
if fid==-1
  fprintf('%s: ERROR: failed to open %s for writing\n',mfilename,fname_out);
  return;
end;
wrote_header_flag = 0;
dirlist = dir(sprintf('%s/*',parms.subjdir));

tmp_subjlist = {dirlist.name};
if isempty(parms.subjlist)
  parms.subjlist = tmp_subjlist;
else
  [parms.subjlist,i_subs] = intersect(parms.subjlist,tmp_subjlist);
  if ~isempty(parms.grouplist)
    parms.grouplist = {parms.grouplist{i_subs}};
  end;
end;
if isempty(parms.subjlist)
  error('subjlist contains no valid subjects');
end;

for i=1:length(parms.subjlist)
  subjname = parms.subjlist{i};
  if ismember(subjname,{'.','..'}), continue; end;
  subjpath = sprintf('%s/%s',parms.subjdir,subjname);
  if ~isdir(subjpath), continue; end;
  if parms.checkstatus_flag
    [allexist,volexist] = fs_check_status(subjname,parms.subjdir,parms.FS_version);
  else
    allexist=1;
    volexist=1;
  end;
  aseg_stats = []; lh_aparc_stats = []; rh_aparc_stats = [];
  if allexist | volexist
    % get cortical thickness, volume for aseg and aparc ROIs
    [aseg_stats,lh_aparc_stats,rh_aparc_stats] = ...
      fs_read_seg_stats(subjname,parms.subjdir,parms.a2005_flag);
    % reduce to subset of ROIs
    if ~isempty(aseg_stats)
      roicodes = cell2mat({aseg_stats.roicode});
      i_roicodes = find(ismember(roicodes,aseg_roilist));
      aseg_stats = aseg_stats(i_roicodes);    
    end;
  else
    fprintf('%s: WARNING: recon incomplete for %s\n',...
      mfilename,subjpath);
    continue;
  end;
  if isempty(aseg_stats) & isempty(lh_aparc_stats) & isempty(rh_aparc_stats)
    fprintf('%s: WARNING: unable to get stats for %s\n',...
      mfilename,subjpath);
    %% todo: output NaN's?
    continue;
  end;
  % check that extracted stats have correct number of ROIs, otherwise exclude from output
  %% todo: could output zeros, but then what if very first subject has wrong number of ROIs
  %%         (must be correct for header)
  if parms.aseg
    tmp_num_aseg_rois = length(aseg_stats);
    if tmp_num_aseg_rois ~= num_aseg_rois
      fprintf('%s: WARNING: number of aseg ROIs (%d) for %s does not match expected (%d)\n',...
        mfilename,tmp_num_aseg_rois,subjname,num_aseg_rois);
      continue;
    end;
  end;
  if parms.aparc
    tmp_num_aparc_rois = length(lh_aparc_stats);
    if tmp_num_aparc_rois ~= num_aparc_rois
      fprintf('%s: WARNING: number of lh aparc ROIs (%d) for %s does not match expected (%d)\n',...
        mfilename,tmp_num_aparc_rois,subjname,num_aparc_rois);
      continue;
    end;
    tmp_num_aparc_rois = length(rh_aparc_stats);
    if tmp_num_aparc_rois ~= num_aparc_rois
      fprintf('%s: WARNING: number of rh aparc ROIs (%d) for %s does not match expected (%d)\n',...
        mfilename,tmp_num_aparc_rois,subjname,num_aparc_rois);
      continue;
    end;
  end;
  if ~wrote_header_flag
    % write header row
    fprintf(fid,'"SubjectID"');
    if ~isempty(parms.grouplist)
      fprintf(fid,',"GroupID"');
    end;
    if parms.aseg
      for k=1:num_aseg_rois
        fprintf(fid,',"%s"',aseg_stats(k).roiname);
      end;
    end;
    if parms.aparc
      for k=1:num_aparc_rois
        for h=1:length(hemilist)
          hemi = hemilist{h};
          switch hemi
            case 'Left'
              tmp_roiname = lh_aparc_stats(k).roiname;
            case 'Right'
              tmp_roiname = rh_aparc_stats(k).roiname;
          end;
          for m=1:length(aparc_meas_list)
            meas = aparc_meas_list{m};
            fprintf(fid,',"%s-ctx-%s-%s"',hemi,tmp_roiname,meas);
          end;
        end;
      end;
    end;
    fprintf(fid,'\n');
    wrote_header_flag=1;
  end;
  % write one row of values
  fprintf(fid,'"%s"',subjname);
  if ~isempty(parms.grouplist)
    fprintf(fid,',"%s"',parms.grouplist{i});
  end;
  if parms.aseg
    for k=1:num_aseg_rois
      fprintf(fid,',%0.6f',aseg_stats(k).volume);
    end;
  end;
  if parms.aparc
    for k=1:num_aparc_rois
      for h=1:length(hemilist)
        hemi = hemilist{h};
        switch hemi
          case 'Left'
            aparc_stats = lh_aparc_stats(k);
          case 'Right'
            aparc_stats = rh_aparc_stats(k);
        end;
        for m=1:length(aparc_meas_list)
          meas = aparc_meas_list{m};
          val = getfield(aparc_stats,meas);
          fprintf(fid,',%0.6f',val);
        end;
      end;
    end;
  end;
  fprintf(fid,'\n');
end;
fclose(fid);

