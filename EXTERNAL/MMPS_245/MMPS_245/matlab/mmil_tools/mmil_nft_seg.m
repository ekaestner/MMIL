function mmil_nft_seg(subj,varargin)
%function mmil_nft_seg(subj,[options])
%
% Purpose: run NFT's Segmentation GUI
%
% Required Input:
%   subj: freesurfer subject name
%
% Optional Parameters:
%   'root_indir': root input directory containing subj
%     {default = pwd}
%   'root_outdir': root output directory
%     if empty, use root_indir
%     {default = []}
%   'outdir': output directory
%     either full path or relative to root_outdir/subj
%     {default = 'bem_nft/seg'}
%   'aseg_flag': [0|1] use aseg to replace nft brain mask
%     with iterative dilation to force separation between layers
%     {default = 0}
%   'fstem_aseg': file stem of aseg
%     {default = 'aparc+aseg'}
%   'aseg_outfix': suffix added to output file name for segments file with aseg
%     {default = 'aseg'}
%   'plot_flag': [0|1] save images of resulting mask volumes
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  02/17/14 by Don Hagler
% Last Mod: 06/10/14 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
parms = check_input(subj,varargin);

% convert data to necessary orientation and format
if parms.verbose
  fprintf('%s: preparing data...\n',mfilename);
end;
parms = prep_data(parms);

% run NFT segmentation
parms = nft_segmentation(parms);

% combine with aseg
if parms.aseg_flag
  parms = combine_aseg(parms);
end;

% save mask images
if parms.plot_flag
  parms = create_images(parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(subj,options)
  parms = mmil_args2parms(options,{...
    'subj',subj,[],...
  ...
    'root_indir',pwd,[],...
    'root_outdir',[],[],...
    'outdir','bem_nft/seg',[],...
    'aseg_flag',false,[false true],...
    'dilate_flag',false,[false true],...
    'fstem_aseg','aparc+aseg',[],...
    'aseg_outfix','aseg',[],...
    'plot_flag',false,[false true],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % undocumented
    'mask_smooth1',40,[],...
    'mask_thresh1',0.5,[],...
    'mask_smooth2',12,[],...
    'mask_thresh2',0.2,[],...
    'mask_smooth3',0,[],...
    'mask_names',{'brainmask','innerskullmask','outerskullmask','scalpmask'},[],...
    'dilate_smooth',25,[],...
    'dilate_thresh',0.2,[],...
    'image_slice',140,[],...
    'visible_flag',false,[false true],...
    'seg_fstem','nu_PSR',[],...
  });

  parms.fname_in = sprintf('%s/%s/mri/nu.mgz',parms.root_indir,parms.subj);
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;

  if parms.aseg_flag
    parms.fname_aseg = sprintf('%s/%s/mri/%s.mgz',...
      parms.root_indir,parms.subj,parms.fstem_aseg);
    if ~exist(parms.fname_aseg,'file')
      error('file %s not found',parms.fname_aseg);
    end;
  end;

  % create output dir
  if isempty(parms.root_outdir)
    parms.root_outdir = parms.root_indir;
  end;
  if mmil_isrelative(parms.outdir)
    parms.outdir = sprintf('%s/%s/%s',...
      parms.root_outdir,parms.subj,parms.outdir);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_data(parms)
  mmil_mkdir(parms.outdir);
  % copy data file
  fname_out = sprintf('%s/nu.mgz',parms.outdir);
  if ~exist(fname_out,'file') || parms.forceflag
    fs_copy_mgh(parms.fname_in,fname_out);
  end;
  fname_in = fname_out;
  % reorient to orientation expected by NFT
  %   use fs_reorient so it is easily reversible (unlike mri_convert)
  fname_out = sprintf('%s/%s.mgz',parms.outdir,parms.seg_fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    [vol,M] = fs_load_mgh(fname_in);
    [vol,M] = fs_reorient(vol,M,'PSR');
    fs_save_mgh(vol,fname_out,M);
  end;
  fname_in = fname_out;
  % convert data
  fname_out = sprintf('%s/%s.img',parms.outdir,parms.seg_fstem);
  fs_mri_convert(fname_in,fname_out,...
    'options','-ot spm','forceflag',parms.forceflag);
  parms.fname_in = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = nft_segmentation(parms)
  fname_out = sprintf('%s/%s_segments.mat',...
    parms.outdir,parms.seg_fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    fprintf('%s: ready to run NFT Segmentation\n',mfilename)
    fprintf('   This requires you to load as input %s.hdr (File->Open menu button)\n',...
      parms.seg_fstem);
    fprintf('     and click "Run" for multiple steps\n');
    fprintf('       setting seed points before each as needed\n');
    fprintf('     and waiting several minutes for the step to complete before clicking "Next"\n');
    fprintf('     For Outer Skull Segmentation, a window will popup\n');
    fprintf('       and you should click on each eye ball\n');
    fprintf('   Finally, you must click "Segmentation" button to save results\n');
    fprintf('     you may click "Output Folder" to set the output directory\n');
    fprintf('   press a key when ready to begin\n',mfilename)
    pause;
    Segmentation('subjectdir',parms.outdir,'subject',parms.seg_fstem);
  end;
  if ~exist(fname_out,'file')
    error('output file %s not found',fname_out);
  end;
  parms.fname_seg = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = combine_aseg(parms)
  % copy aseg
  if parms.verbose
    fprintf('%s: creating mask from aseg...\n',mfilename);
  end;
  fname_out = sprintf('%s/%s.mgz',parms.outdir,parms.fstem_aseg);
  if ~exist(fname_out,'file') || parms.forceflag
    fs_copy_mgh(parms.fname_aseg,fname_out);
  end;
  parms.fname_aseg = fname_out;
  % create brain mask from aseg
  parms = create_aseg_brain_mask(parms);
  % iteratively dilate mask
  if parms.verbose
    fprintf('%s: dilating aseg brain mask...\n',mfilename);
  end;
  parms = dilate_aseg_brain_mask(parms);  
  % combine with masks in Segm structure
  if parms.verbose
    fprintf('%s: combining masks with segmentation...\n',mfilename);
  end;
  parms = combine_segments(parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_aseg_brain_mask(parms)
  % create tight brain mask from aseg
  tmp_parms = [];
  tmp_parms.smooth1 = parms.mask_smooth1;
  tmp_parms.thresh1 = parms.mask_thresh1;
  tmp_parms.smooth2 = parms.mask_smooth2;
  tmp_parms.thresh2 = parms.mask_thresh2;
  tmp_parms.smooth3 = parms.mask_smooth3;
  tmp_parms.forceflag = parms.forceflag;
  tmp_parms.fname_in = parms.fname_aseg;
  tmp_parms.fname_out = sprintf('%s/%s_mask0.mgz',...
    parms.outdir,parms.fstem_aseg);
  args = mmil_parms2args(tmp_parms);
  mmil_dilate_mask([],args{:});
  parms.fname_brain_mask = tmp_parms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = dilate_aseg_brain_mask(parms)
  tmp_parms = [];
  tmp_parms.smooth1 = parms.dilate_smooth;
  tmp_parms.thresh1 = parms.dilate_thresh;
  tmp_parms.smooth2 = 0;
  tmp_parms.smooth3 = 0;
  tmp_parms.forceflag = parms.forceflag;
  if parms.dilate_flag
    nmask = 3;
  else
    nmask = 2;
  end;
  for j=0:nmask
    k = j + 1;
    tmp_parms.fname_in = sprintf('%s/%s_mask%d.mgz',...
      parms.outdir,parms.fstem_aseg,j);
    tmp_parms.fname_out = sprintf('%s/%s_mask%d.mgz',...
      parms.outdir,parms.fstem_aseg,k);
    args = mmil_parms2args(tmp_parms);
    mmil_dilate_mask([],args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = combine_segments(parms)
  parms.seg_fstem = [parms.seg_fstem '_' parms.aseg_outfix];
  fname_out = sprintf('%s/%s_segments.mat',...
    parms.outdir,parms.seg_fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    % load existing seg struct
    load(parms.fname_seg);
    if parms.dilate_flag
      m = 1;
    else
      m = 0;
    end;
    nmasks = length(parms.mask_names);
    % combine existing masks with aseg-derived masks
    for i=1:nmasks
      mask_name = parms.mask_names{i};
      fname_mask = sprintf('%s/%s_mask%d.mgz',...
        parms.outdir,parms.fstem_aseg,m);
      m = m + 1;
      vol_mask = mmil_reorient_for_nft(fname_mask,...
        parms.outdir,parms.forceflag);
      if i==1
        Segm.(mask_name) = (vol_mask>0);
      else
        Segm.(mask_name) = (Segm.(mask_name) | (vol_mask>0));
      end;
    end;
    save(fname_out,'Segm');
  end;
  parms.fname_seg = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_images(parms)
  fname_out = sprintf('%s/%s_sag%d.tif',...
    parms.outdir,parms.seg_fstem,parms.image_slice);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: saving mask image...\n',mfilename);
    end;
    load(parms.fname_seg);
    im = 0;
    for i=1:length(parms.mask_names)
      im = im + squeeze(Segm.(parms.mask_names{i})(parms.image_slice,:,:));
    end;
    figure(1); clf;
    imagesc(im',[0,4]);
    axis xy;
    if ~parms.visible_flag
      set(gcf,'visible','off');
    end;
    print(gcf,'-dtiff',fname_out);
    if ~parms.visible_flag
      close(gcf);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


