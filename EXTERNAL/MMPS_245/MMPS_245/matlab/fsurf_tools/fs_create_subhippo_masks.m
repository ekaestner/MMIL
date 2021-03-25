function flist=fs_create_subhippo_masks(subjname,subjdir,outdir,listflag,forceflag)
%function flist=fs_create_subhippo_masks(subjname,subjdir,outdir,listflag,forceflag)
%
% Required:
%   subjname: freesurfer recon subject name
%
% Optional:
%   subjdir: root freesurfer subjects dir
%       if empty, will use getenv('SUBJECTS_DIR')
%     {default = []}
%   outdir: output directory
%      if empty, will put output in mri subdirectory in subjdir
%      if not full path, will interpret as subdirecory in subjdir
%     {default = []}
%   listflag: [0|1] whether to only list roi names rather
%      than try to create them
%     {default = 0}
%   forceflag: [0|1] whether to overwrite existing output
%     {default = 0}
%
% created:  04/01/09 by Don Hagler
% last mod: 04/01/09 by Don Hagler
%

flist = [];
if (~mmil_check_nargs(nargin,2)) return; end;
if ~exist('outdir','var') | isempty(outdir), outdir = 'mri'; end
if ~exist('listflag','var') | isempty(listflag), listflag = 0; end
if ~exist('forceflag','var') | isempty(forceflag), forceflag = 0; end

hemilist = {'lh','rh'};
roicodes = [17,53];
roiname = 'hippomask';
num_div_z = 3; % x = L-R, y = S-I, z = P-A

subjpath = sprintf('%s/%s',subjdir,subjname);
if ~exist(subjpath,'file')
  error('freesurfer recon dir %s not found',subjpath);
end;

if ~strncmp(outdir,'/',1)
  outdir = [subjpath '/' outdir];
end;
if ~exist(outdir,'dir')
  [success,msg]=mkdir(outdir);
  if ~success
    error('failed to created output directory %s:\n%s',...
      outdir,msg);
  end;
end;

if listflag
  f=1;
  for h=1:length(hemilist)
    hemi = hemilist{h};
    for r=1:num_div_z
      fname_mask = sprintf('%s/%s_%s_%dof%d.mgz',...
        outdir,hemi,roiname,r,num_div_z);
      flist{f} = fname_mask;
      f = f+1;
    end;
  end;
else
  % affine transform of aseg to atlas space
  fname_in = sprintf('%s/mri/aseg.mgz',subjpath);
  fname_out = sprintf('%s/aseg_lta.mgz',outdir);
  fname_trans = sprintf('%s/mri/transforms/talairach.lta',subjpath);
  status=fs_warp_vol2atlas(subjname,fname_in,fname_out,...
    'transform',fname_trans,'resamp_type','nearest',...
    'subjdir',subjdir,'overwrite_flag',forceflag);
  if status
    fprintf('%s: ERROR: Failed to transforma aseg to atlas space\n',mfilename);
    return;
  end;

  % load aseg volume
  [aseg_vol,M_atlas,mrparms,volsz] = fs_load_mgh(fname_out);

  f=1;
  for h=1:length(hemilist)
    hemi = hemilist{h};
    roicode = roicodes(h);

    % extract hippocampi
    roi = find(aseg_vol==roicode);
    if isempty(roi)
      fprintf('%s: ERROR: missing roicode %d from aseg',mfilename,roicode);
      flist = [];
      return;
    end;  

    % find A-P boundaries
    [i,j,k] = ind2sub(volsz(1:3),roi);

    % calculate cumulative volume as a function of z
    tmp_vol = zeros(volsz);
    tmp_vol(roi) = 1;
    vol_vec = squeeze(sum(sum(tmp_vol,1),2));
    cum_vol_vec = cumsum(vol_vec)/sum(vol_vec);

    % create masks for hippocampus subdivided by fractional volume
    for r=1:num_div_z
      fname_mask = sprintf('%s/%s_%s_%dof%d.mgz',...
        outdir,hemi,roiname,r,num_div_z);
      if ~exist(fname_mask,'file') || forceflag
        if r==1
          tmp = find(cum_vol_vec>0);
          lb = tmp(1);
        else
          tmp = (r-1)/num_div_z;
          tmp = cum_vol_vec-tmp;
          [tmp,lb] = min(abs(tmp));
        end;
        lb = max(1,lb);

        if r==num_div_z
          tmp = find(cum_vol_vec==1);
          ub = tmp(1);
        else
          tmp = r/num_div_z;
          tmp = cum_vol_vec-tmp;
          [tmp,ub] = min(abs(tmp));
        end;
        ub = min(256,ub);

        p = find(k>=lb & k<=ub);
        tmp_roi = roi(p);
        tmp_vol = zeros(volsz);
        tmp_vol(tmp_roi) = 1;

        % save as individual masks
        fname_tmp = sprintf('%s/%s_%s_%dof%d_atlas.mgh',...
          outdir,hemi,roiname,r,num_div_z);
        fs_save_mgh(tmp_vol,fname_tmp,M_atlas);

        % transform back to single subject space
        fname_mask = sprintf('%s/%s_%s_%dof%d.mgz',...
          outdir,hemi,roiname,r,num_div_z);
        status=fs_warp_vol2atlas(subjname,fname_tmp,fname_mask,...
          'transform',fname_trans,'resamp_type','nearest',...
          'inverse_flag',1,...
          'subjdir',subjdir,'overwrite_flag',forceflag);
        if status
          fprintf('%s: ERROR: Failed to transform mask to subject space\n',mfilename);
          flist = [];
          delete(fname_tmp);
          return;
        end;
        delete(fname_tmp);
      end;
      flist{f} = fname_mask;
      f = f+1;
    end;
  end;
end;
