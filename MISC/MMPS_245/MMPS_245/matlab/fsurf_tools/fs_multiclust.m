function fs_multiclust(fname,varargin)
%function fs_multiclust(fname,varargin)
%
% Purpose: find cortical surface clusters using sliding
%   threshold, for multiple frames
%
% Required Input:
%   fname: full path file name of mgh format surface data
%     may have multiple frames (e.g. time points or conditions)
%
% Optional Parameters:
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'mclust'}
%   'frames': index numbers for frames of fname
%     if empty, use all frames
%     {default = []}
%   'subj': freesurfer subject
%     {default = 'fsaverage'}
%   'subjdir': freesurfer SUBJECTS_DIR
%     if empty, use FREESURFER_HOME/subjects
%     {default = []}
%   'hemi': cortical hemisphere
%     If empty will attempt to infer from fname
%     {default = []}
%   'val_thresh': vector of thresholds for which to find clusters
%     {default = [2:15]}
%   'clust_thresh': minimum cluster size (mm^2)
%     {default = 100}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  11/15/11 by Don Hagler
% Last Mod: 11/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms_filter = {...
  'fname',fname,[],...
  'outdir',pwd,[],...
  'outstem','mclust',[],...
  'frames',[],[],...
  'subj','fsaverage',[],...
  'subjdir',[],[],...
  'hemi',[],[],...
  'val_thresh',[2:15],[],...
  'clust_thresh',100,[],...
  'forceflag',false,[false true],...
...
  'surfname','white',[],...
  'growflag',false,[false true],...
  'growthresh',2,[],...
  'trimflag',false,[false true],...% 1: trim larger clusters, 0: remove larger clusters
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = check_input(varargin,parms_filter);

run_multiclust(parms);

sort_labels(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options,parms_filter)
  parms = mmil_args2parms(options,parms_filter);
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;
  if isempty(parms.subjdir)
    parms.subjdir = sprintf('%s/subjects',getenv('FREESURFER_HOME'));
  end;
  if isempty(parms.hemi)
    n = regexp(parms.fname,'\w+-(?<hemi>[lr]h)','names');
    if isempty(n)
      error('unspecified hemi for %s',parms.fname);
    end;
    parms.hemi = n.hemi;
  end;
  [parms.indir,parms.fstem] = fileparts(parms.fname);
  parms.nthresh = length(parms.val_thresh);
  if length(parms.clust_thresh)==1 & parms.nthresh~=1
    parms.clust_thresh = parms.clust_thresh * ones(size(parms.val_thresh));
  end;
  tmp = fs_load_subj(parms.subj,parms.hemi,parms.surfname,1,parms.subjdir);
  parms.nverts = tmp.nverts;
  [tmp,volsz] = fs_read_header(parms.fname);
  if volsz(1) ~= parms.nverts
    error('mismatch between number of vertices in fname %s (%d) and surface for subj %s',...
      parms.fname,volsz(1),parms.subj);
  end;
  if isempty(parms.frames)
    parms.frames = [1:volsz(4)];
  end;
  parms.nframes = numel(parms.frames);
  parms.extradir = [parms.outdir '/extra'];
  mmil_mkdir(parms.extradir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_multiclust(parms)
  for i=1:parms.nframes
    f = parms.frames(i);
    fname_out = sprintf('%s/%s-clust-f%d-%s.mgh',parms.extradir,parms.fstem,f,parms.hemi);
    if ~exist(fname_out,'file') || parms.forceflag
      fname_log = sprintf('%s/%s-clust-f%d-%s.log',parms.extradir,parms.fstem,f,parms.hemi);
      if ~exist(fname_log,'file') || parms.forceflag
        % convert to w format
        instem = sprintf('%s-f%d',parms.fstem,f);
        outstem = sprintf('%s-multiclust-f%d',parms.fstem,f);
        tmp_fstem = sprintf('%s/%s',parms.extradir,instem);
        tmp_fname = sprintf('%s-%s.w',tmp_fstem,parms.hemi);
        if ~exist(tmp_fname,'file') || parms.forceflag
          fs_mgh2w(parms.fname,tmp_fstem,parms.hemi,f,f,1);
        end;
        % run multiclust
        cmd = sprintf('setenv SUBJECTS_DIR %s',parms.subjdir);
        cmd = sprintf('%s\nmulticlust',cmd);
        cmd = sprintf('%s -instem %s',cmd,instem);
        cmd = sprintf('%s -subj %s',cmd,parms.subj);
        cmd = sprintf('%s -outstem %s',cmd,outstem);
        cmd = sprintf('%s -outdir %s',cmd,parms.extradir);
        cmd = sprintf('%s -indir %s',cmd,parms.extradir);
        cmd = sprintf('%s -hemi %s',cmd,parms.hemi);
        cmd = sprintf('%s -nthresh %d',cmd,parms.nthresh);
        cmd = sprintf('%s -thresh %s',cmd,sprintf(' %0.6f',parms.val_thresh));
        cmd = sprintf('%s -minarea %s',cmd,sprintf(' %0.6f',parms.clust_thresh));
        if parms.growflag
          cmd = sprintf('%s -growthresh %0.6f',cmd,parms.growthresh);
        else
          cmd = sprintf('%s -nogrow',cmd);
        end;      
        cmd = sprintf('%s -surf %s',cmd,parms.surfname);
        cmd = sprintf('%s -summary',cmd);
        cmd = sprintf('%s -outputmasks',cmd);

        fprintf('%s: cmd = %s\n',mfilename,cmd);
        [s,r] = mmil_unix(cmd);
        if s, error('failed:\n%s',r); end;
        fprintf('%s\n',r);

        fid = fopen(fname_log,'wt');
        if fid<0, error('failed to open %s for writing',fname_log); end;
        fprintf(fid,'%s\n\n%s\n',cmd,r);
        fclose(fid);        
      end;
      % combine into a single mgh file for each time point
      flist = dir(sprintf('%s/%s-*-%s.w',parms.extradir,outstem,parms.hemi));
      nclust = length(flist);
      vals = zeros(parms.nverts,1);
      for c=1:nclust
        fname = char(flist(c).name);
        regpat = sprintf('%s-(?<clustnum>\\d+)-%s.w',outstem,parms.hemi);
        n = regexp(fname,regpat,'names');
        clustnum = str2double(n.clustnum) + 1;
        fname = sprintf('%s/%s',parms.extradir,fname);
        [w,v] = fs_read_wfile(fname);
        vals(v) = clustnum;
      end;
      fs_save_mgh(vals,fname_out);
    end;
  end;

  % combine into time course of masks
  fname_out = sprintf('%s/%s-clust_tc-%s.mgh',parms.extradir,parms.fstem,parms.hemi);
  if ~exist(fname_out,'file') || parms.forceflag
    vals = zeros(parms.nverts,1,1,parms.nframes);
    for i=1:parms.nframes
      f=parms.frames(i);
      parms.fname = sprintf('%s/%s-clust-f%d-%s.mgh',...
        parms.extradir,parms.fstem,f,parms.hemi);
      tmp_vals = fs_load_mgh(parms.fname);
      vals(:,1,1,i) = tmp_vals;
    end;
    fs_save_mgh(vals,fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = sort_labels(parms)
  clust_info = [];
  fname_in = sprintf('%s/%s-clust_tc-%s.mgh',...
    parms.extradir,parms.fstem,parms.hemi);
  fname_out = sprintf('%s/%s-%s.mgh',...
    parms.outdir,parms.outstem,parms.hemi);
  fname_mat = sprintf('%s/%s-%s.mat',...
    parms.outdir,parms.outstem,parms.hemi);
  fname_annot = sprintf('%s/%s.%s_fparc.annot',...
    parms.outdir,parms.hemi,parms.outstem);
  if ~exist(fname_out,'file') ||...
     ~exist(fname_mat,'file') ||...
     ~exist(fname_annot,'file') || parms.forceflag
    vals = squeeze(fs_load_mgh(fname_in));
    ind_labeled = find(sum(vals,2)>0);
    nlabeled = length(ind_labeled);
    labeled_vals = vals(ind_labeled,:);

    % gather information about all clusters
    clust_info = [];
    c = 1;
    for t=1:parms.nframes
      tmp_vals = labeled_vals(:,t);
      uniq_vals = unique(tmp_vals(tmp_vals>0));
      nclust = length(uniq_vals);
      for i=1:nclust
        ind = find(tmp_vals==uniq_vals(i));
        v = ind_labeled(ind);
        clust_info(c).nverts = numel(v);
        clust_info(c).verts = v;
        clust_info(c).t = t;
        c = c + 1;
      end;
    end;
    nclust = length(clust_info);

    % sort ascending order by size
    [sorted_sizes,ind_sort] = sort([clust_info.nverts]);
    clust_info = clust_info(ind_sort);

    % check for overlap with larger clusters
    for i=1:nclust
      for j=i+1:nclust  
        [overlap,ii,ij] = intersect(clust_info(i).verts,clust_info(j).verts);
        if ~isempty(overlap)
          if parms.trimflag  % remove verts from larger cluster
            [v_keep,i_keep] = setdiff(clust_info(j).verts,overlap);
            clust_info(j).verts = v_keep;
            clust_info(j).nverts = numel(v_keep);
            fprintf('%s: removing %d verts from cluster %d...\n',...
              mfilename,length(overlap),j);
          else          % remove larger cluster
            clust_info(j).verts = [];
            clust_info(j).nverts = 0;
            fprintf('%s: removing cluster %d...\n',...
              mfilename,j);
          end;
        end;
      end;
    end;

    % exclude clusters with 0 vertices (removed overlap)
    ind_nonzero = find([clust_info.nverts]);
    clust_info = clust_info(ind_nonzero);
    nclust = length(clust_info);

    % save mgh file with all clusters
    new_vals = zeros(parms.nverts,1);
    for i=1:nclust
      v = clust_info(i).verts;
      new_vals(v) = i;
    end;
    fs_save_mgh(new_vals,fname_out);

    % save label files
    new_vals = zeros(parms.nverts,1);
    roinames = cell(nclust,1);
    for i=1:nclust
      roinames{i} = sprintf('%s_roi%02d',parms.outstem,i);
      fname_label = sprintf('%s/%s.%s.label',...
        parms.outdir,parms.hemi,roinames{i});
      v = clust_info(i).verts;
      fs_write_label(v,fname_label)
    end;

    % create color table file
    fname_ctab = sprintf('%s/%s-%s.ctab',...
      parms.outdir,parms.outstem,parms.hemi);
    if ~exist(fname_ctab,'file') || parms.forceflag
      fs_write_ctab(roinames,fname_ctab);
    end;

    % convert labels to annotation
    fs_labels2annot(parms.outdir,...
      'fname_ctab',fname_ctab,'outdir',parms.outdir,...
      'source_subj',parms.subj,'subj',parms.subj,...
      'subjdir',parms.subjdir,'annotname',[parms.outstem '_fparc'],...
      'hemilist',{parms.hemi},...
      'forceflag',parms.forceflag);

    % print summary
    fprintf('Clust #     nverts\n');
    for i=1:nclust
      fprintf('%3d       %10d\n',i,clust_info(i).nverts);
    end;
    save(fname_mat,'clust_info','parms');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

