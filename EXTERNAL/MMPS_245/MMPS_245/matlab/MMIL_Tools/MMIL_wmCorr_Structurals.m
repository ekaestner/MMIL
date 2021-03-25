function errcode = MMIL_wmCorr_Structurals(ContainerPath,varargin)
%function errcode = MMIL_wmCorr_Structurals(ContainerPath,[options])
%
% Purpose: correct structural images for bias-field intensity inhomogeneities
%   using nu correction (N3) and white matter segmentation-based bias correction
%
% Required Input:
%   ContainerPath: full path of MRIPROC container
%
% Optional Input:
%  'scantypes': cell array of scan types to correct
%    {default = {'MPR','XetaT2','FLASHhi','MEDIChi',}}
%  'abcd_T1_flag': [0|1] use customized ABCD options for T1 correction
%    if true, some options below are ignored
%    {default = 1}
%  'n3flag': [0|1] apply N3 bias correction as pre-normalization step
%    {default = 1}
%  'vol_target': target volume for white matter intensity
%    may be a scalar value
%    if empty, will use median white matter intensity value
%    {default = 110}
%  'ratio_range': range of values to constrain ratio
%    if empty, no clipping
%    {default = []}
%  'image_range': range of values to clip output image
%    {default = [0,255]}
%  'T2w_image_range': range of values to clip output for T2-weighted image
%    {default = [0,2500]}
%  'lambda2': weighting factor for Laplacian of smooth volume
%    higher values mean stronger smoothness constraint
%    {default = 5}
%  'niters_spsm': number of iterations of sparse smoothing
%    {default = 3}
%  'cleanup_flag': [0|1] remove temporary files before quitting
%    {default = 1}
%  'verbose': [0|1] display status messages
%    {default = 0}
%  'ext': file extension (i.e. '.mgh' or '.mgz')
%    {default = '.mgz'}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Created:  11/30/15 by Don Hagler
% Prev Mod: 08/02/17 by Don Hagler
% Last Mod: 11/06/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'scantypes',{'MPR','XetaT2','FLASHhi','MEDIChi'},[],...
  'abcd_T1_flag',true,[false true],...
  'n3flag',true,[false true],...
  'vol_target',110,[],...
  'ratio_range',[],[],...
  'image_range',[0,255],[],...
  'T2w_image_range',[0,2500],[],...
  'lambda2',5,[0,1e10],...
  'niters_spsm',3,[1,10],...
  'tmpdir','tmp_wm_corr',[],...
  'cleanup_flag',true,[false true],...
  'verbose',false,[false true],...
  'ext','.mgz',{'.mgh','.mgz'},...
  'forceflag',false,[false true],...
...
  'outfix','_wmbc',[],...
  'wm_outfix','_wmseg',[],...
  'bm_outfix','_brain',[],...
  'T2w_scantypes',{'XetaT2'},[],...
...
  'interpm',2,[0:5],...
  'smoothmask_flag',true,[false,true],...
...
  'corr_tags',{'fname_out','fname_wm','fname_bm','n3flag',...
               'fname_n3','vol_target','ratio_range','image_range',...
               'lambda0','lambda1','lambda2','niters_spsm','tmpdir',...
               'cleanup_flag','verbose','forceflag','estop','maxiter',...
               'smooth1','thresh1','smooth2','thresh2','smooth3','thresh3',...
               'ext'},[],...
  'abcd_T1_corr_tags',{'fname_out','fname_wm','fname_bm','image_range','ext'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load container info
[ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
if errcode ~= 0, return; end;

if isempty(ContainerInfo)
  fprintf('ERROR: ContainerInfo is empty\n');
  errcode = 1;
  return;
end

if ~isfield(ContainerInfo,'ScanInfo')
  fprintf('ERROR: ContainerInfo.ScanInfo missing\n');
  errcode = 1;
  return;
end
ScanInfo = ContainerInfo.ScanInfo;

% suffixlist allows to check whether input file exist
% checks for the input file in the order specified
suffixlist = {'_B1_uw','B1','_uw',''};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply wm bias corrections
fname_ref = []; fname_ref_wm = []; fname_ref_bm = [];
for f=1:length(parms.scantypes)
  fnamestem = parms.scantypes{f};
  scaninfo = mmil_getfield(ScanInfo,fnamestem,[]);
  T2flag = ismember(fnamestem,parms.T2w_scantypes);
  for i = 1:length(scaninfo)
    for s = 1:length(suffixlist)
      fname_in = sprintf('%s/%s%d%s%s',...
        ContainerPath,fnamestem,i,suffixlist{s},parms.ext);
      fname_out = sprintf('%s/%s%d%s%s%s',...
        ContainerPath,fnamestem,i,suffixlist{s},parms.outfix,parms.ext);
      fname_wm = sprintf('%s/%s%d%s%s%s',...
        ContainerPath,fnamestem,i,suffixlist{s},parms.wm_outfix,parms.ext);
      fname_bm = sprintf('%s/%s%d%s%s%s',...
        ContainerPath,fnamestem,i,suffixlist{s},parms.bm_outfix,parms.ext);
      if exist(fname_in, 'file')
        if isempty(fname_ref) && ~T2flag
          refflag = 1;
          fname_ref = fname_in;
          fname_ref_wm = fname_wm;
          fname_ref_bm = fname_bm;
        end;
        tparms = parms;
        tparms.fname_out = fname_out;
        tparms.fname_wm = fname_wm;
        tparms.fname_bm = fname_bm;
        tparms.tmpdir = sprintf('%s/tmp_wmbc_%s%d',...
          ContainerPath,fnamestem,i);
        if ~exist(fname_out,'file') || parms.forceflag
          % if scan is T2w image, register to reference T1w image
          if T2flag
            tparms.image_range = parms.T2w_image_range;
            if isempty(fname_ref_wm) || isempty(fname_ref_bm) ||...
               ~exist(fname_ref_wm,'file') || ~exist(fname_ref_bm,'file')
               fprintf('%s: missing reference masks for %s, skipping...\n',...
                 mfilename,fname_in);
              errcode = 1;
              continue;
            end;
            fprintf('%s: registering %s to %s...\n',...
              mfilename,fname_in,fname_ref);
            M_ref_to_orig = mmil_jpdfreg_T1T2(fname_ref,fname_in,...
              'outdir',tparms.tmpdir,'fname_T1_mask',fname_ref_bm,...
              'smoothmask_flag',parms.smoothmask_flag,...
              'cleanup_flag',parms.cleanup_flag,'forceflag',parms.forceflag);
            resample_mask(fname_ref_wm,tparms.fname_wm,...
              fname_in,M_ref_to_orig,parms,1);
            resample_mask(fname_ref_bm,tparms.fname_bm,...
              fname_in,M_ref_to_orig,parms,0);
          end;
          fprintf('%s: wm bias correction %s...\n',mfilename,fname_in);
          if parms.abcd_T1_flag && ~T2flag
            args = mmil_parms2args(tparms,parms.abcd_T1_corr_tags);
            abcd_T1_wm_corr(fname_in,args{:});
          else
            args = mmil_parms2args(tparms,parms.corr_tags);
            mmil_wm_corr(fname_in,args{:});
          end;
        end;
        % cleanup
        if parms.cleanup_flag && exist(tparms.tmpdir,'dir')
          cmd = sprintf('rm -r %s',tparms.tmpdir);
          [status,msg] = unix(cmd);
          if status
            error('tmpdir cleanup failed:\n%s\n%s\n',cmd,msg);
          end;
        end;
        break; % breaks from suffixlist loop after first fname_in occurrence
      end;
    end;  
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resample_mask(fname_in,fname_out,fname_targ,M_reg,parms,wm_flag)
  if ~exist(fname_out,'file') || parms.forceflag
    fprintf('%s: resampling %s into alignment with %s...\n',...
      mfilename,fname_in,fname_targ);
    [vol,M] = fs_load_mgh(fname_in);
    indir = fileparts(fname_targ);
    [M_targ,volsz_targ] = mmil_load_mgh_info(fname_targ,...
      parms.forceflag,indir);
    [vol,M] = mmil_resample_vol(vol,M,...
      'M_ref',M_targ,'nvox_ref',volsz_targ(1:3),...
      'interpm',parms.interpm,'bclamp',1,...
      'M_reg',inv(M_reg));
    vol = 1.0*(vol>0.5);
    if parms.abcd_T1_flag
      % if abcd_T1_flag, no erosion or dilation was done
      %   so do it now for use with T2
      vol_ctx = ctx_mgh2ctx(vol,M);
      if wm_flag
        if parms.verbose
          fprintf('%s: eroding white matter mask...\n',mfilename);
        end;
        vol_ctx = abcd_erode_wm(vol_ctx);
      else
        if parms.verbose
          fprintf('%s: dilating brain mask...\n',mfilename);
        end;
        vol_ctx = abcd_dilate_bm(vol_ctx);
      end;
      vol = vol_ctx.imgs;
    end;  
    fs_save_mgh(vol,fname_out,M);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

