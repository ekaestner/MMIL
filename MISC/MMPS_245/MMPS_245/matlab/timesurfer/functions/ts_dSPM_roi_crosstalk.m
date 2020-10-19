function ts_dSPM_roi_crosstalk(varargin)
%function ts_dSPM_roi_crosstalk([options])
%
% Purpose: calculate dSPM crosstalk maps for regions of interest
%
% Optional Parameters:
%  'rootdir': root directory containing stcfiles directory
%    {default: pwd} (current working directory)
%  'indir': input directory, relative to rootdir unless full path is given
%    {default: 'matfiles'}
%  'prefix': prefix of stc files
%    If full path is given, rootdir is used only for output
%    If cell array, will loop over multiple prefixes like different conditions
%    {default: 'dSPM'}
%  'outdir': output directory, relative to rootdir unless full path is given
%    {default: 'roi_crosstalk'}
%  'outstem': output file stem
%    {default: 'CT'}
%  'roidir': full path of directory containing label files
%    If empty, looks in subjdir/subjname/label
%    {default = []}
%  'subjname': name of FreeSurfer subject with label files (ROIs)
%    {default = fsaverage}
%  'subjdir': FreeSurfer root directory containing subjname
%    {default = $SUBJECTS_DIR}
%  'ico': icosahedron order number:
%      Order     Number of Vertices
%        1              42
%        2             162
%        3             642
%        4            2562
%        5           10242
%        6           40962
%        7          163842
%     If 0, assumes stcfiles are in native subject vertices
%    {default: 4}
%  'roinames': cell array of ROI names
%      Label files have this format: hemi.roiname.label or 
%                                    hemi.roiname.ico4.label
%      If empty, will search for label files in subjdir/subjname/label
%      If none found, will quit with error
%    {default: []}
%  'hemilist': cell array of cortical hemispheres
%    {default: {'lh','rh'}}
%  'forceflag': [0|1] overwrite existing output
%    {default: 0}
%
% Created:   11/29/10 by Don Hagler
% Last Mod:  11/30/10 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'rootdir',pwd,[],...
  'indir','matfiles',[],...
  'prefix','dSPM',[],...
  'conditions',[],[],...
  'outdir','roi_crosstalk',[],...
  'outstem','CT',[],...
  'roidir',[],[],...
  'subjname','fsaverage',[],...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'ico',4,[0:7],...
  'roinames',[],[],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check indir
if isempty(regexp(parms.indir,'^/')) % relative
  parms.indir = [parms.rootdir '/' parms.indir];
end;
if ~exist(parms.indir,'dir')
  error('input dir %s not found',parms.indir);
end;

% check roi dir
if isempty(parms.roidir)
  parms.roidir = [parms.subjdir '/' parms.subjname '/label'];
end;
if ~exist(parms.roidir,'dir')
  error('dir %s not found',parms.roidir);
end;

% find roi names
if isempty(parms.roinames)
  flist = dir(sprintf('%s/*.label',parms.roidir));
  if isempty(flist)
    error('no label files found in roidir %s',parms.roidir);
  end;
  tmp_roinames = {flist.name};
  if parms.ico>0
    pat = sprintf('^[lr]h.(?<roiname>.+).ico%d.label$',parms.ico);
  else
    pat = '^[lr]h.(?<roiname>.+).label$';
  end;
  n = regexp(tmp_roinames,pat,'names');
  keep_flags = ~cellfun(@isempty,n);
  ind_keep = find(keep_flags);
  for k=1:length(ind_keep)
    parms.roinames{end+1} = n{ind_keep(k)}.roiname;
  end;
  parms.roinames = unique(parms.roinames);
end;

% check roi files
if isempty(parms.roinames)
  error('roinames is empty');
else
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    for r=1:length(parms.roinames)
      roiname = parms.roinames{r};
      if parms.ico
        fname_roi = sprintf('%s/%s.%s.ico%d.label',...
          parms.roidir,hemi,roiname,parms.ico);
      else
        fname_roi = sprintf('%s/%s.%s.label',...
          parms.roidir,hemi,roiname);
      end;
      if ~exist(fname_roi,'file')
        error('file %s not found',fname_roi);
      end;
    end;
  end;
end;

% create outdir
if isempty(regexp(parms.outdir,'^/')) % relative
  parms.outdir = [parms.rootdir '/' parms.outdir];
end;
if ~exist(parms.outdir,'dir')
  [success,msg] = mkdir(parms.outdir);
  if ~success
    error('failed to created %s:\n%s',parms.outdir,msg);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load parms and inverse
fname = sprintf('%s/%s_parms.mat',parms.indir,parms.prefix);
if ~exist(fname,'file'), error('file %s not found',fname); end;
dSPM_parms = load(fname);
dSPM_parms = dSPM_parms.parms;

fname = sprintf('%s/%s_inverse.mat',parms.indir,parms.prefix);
if ~exist(fname,'file'), error('file %s not found',fname); end;
load(fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_dips_lh = find(dSPM_parms.lh_dec_dips);
v_dips_rh = find(dSPM_parms.rh_dec_dips);
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  for r=1:length(parms.roinames)
    roiname = parms.roinames{r};
    if parms.ico
      fname_roi = sprintf('%s/%s.%s.ico%d.label',...
        parms.roidir,hemi,roiname,parms.ico);
    else
      fname_roi = sprintf('%s/%s.%s.label',...
        parms.roidir,hemi,roiname);
    end;
    runflag = 0;
    for h2=1:length(parms.hemilist)
      hemi2 = parms.hemilist{h2};
      fname_out = sprintf('%s/%s-%s_%s-%s.mgh',...
          parms.outdir,roiname,hemi,parms.outstem,hemi2);
      if ~exist(fname_out,'file') || parms.forceflag
        runflag=1;
        break;
      end;
    end;
    if runflag
      fprintf('%s: calculating crosstalk for ROI %s-%s...\n',...
        mfilename,roiname,hemi);
      roi = fs_read_label(fname_roi);
      avg_CT = single(0);
      for k=1:length(roi)
        v = roi(k);
        switch hemi
          case 'lh'
            p = 0;
            q = find(v==v_dips_lh);
          case 'rh'
            p = dSPM_parms.num_dips_lh;
            q = find(v==v_dips_rh);
        end;
        if isempty(q), error('bad vertex number %d',v); end;
        p = 3*(p + q - 1) + 1;
        q = p + 2;      
        tmp_M = M(p:q,:);
        tmp_G = G_xyz(:,p:q);

        % multiply forward matrix by inverse for current vertex
        tmp_numer = (tmp_M*G_xyz).^2;
        tmp_denom = (tmp_M*tmp_G).^2;
        % sum over dipole xyz
        tmp_numer = sum(tmp_numer,1);
        tmp_numer = tmp_numer(1:3:end) + tmp_numer(2:3:end) + tmp_numer(3:3:end);
        tmp_denom = sum(tmp_denom(:));
        % calculate normalized crosstalk
        tmp_CT = tmp_numer/tmp_denom;
        avg_CT = avg_CT + tmp_CT;
      end;
      avg_CT = avg_CT/length(roi);
      for h2=1:length(parms.hemilist)
        hemi2 = parms.hemilist{h2};
        fname_out = sprintf('%s/%s-%s_%s-%s.mgh',...
            parms.outdir,roiname,hemi,parms.outstem,hemi2);
        switch hemi2
          case 'lh'
            vals = avg_CT(1:dSPM_parms.num_dec_dips_lh)';
          case 'rh'
            vals = avg_CT(dSPM_parms.num_dec_dips_lh+1:end)';
        end;
        fprintf('%s: writing %s...\n',mfilename,fname_out);
        fs_save_mgh(vals,fname_out);
      end;
    end;
  end;
end;

