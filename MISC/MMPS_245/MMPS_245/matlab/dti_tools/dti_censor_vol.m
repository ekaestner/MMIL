function vol = dti_censor_vol(vol,vol_synth,censor_mat,varargin)
%function vol = dti_censor_vol(vol,vol_synth,censor_mat,[options])
%
% Purpose: Replace values in diffusion images with tensor synthesized values
%
% Required Parameters:
%   vol: 4D diffusion data                            size = [nx,ny,nz,nf]
%   vol_synth: synthesized 4D diffusion data          size = [nx,ny,nz,nf]
%   censor_mat: matrix of censored slices by frames   size = [nz,nf]
%
% Optional Parameters:
%  'volmask': 3D brain mask volume                     size = [nx,ny,nz]
%    restrict replacement to voxels inside brainmask
%    {default=[]}
%  'censor_min_ndirs': minimum number of diffusion directions (not including
%     b=0 images) required for tensor fit after censoring
%     will not do censoring if it means reducing number of directions below min
%     {default=6}
%  'nb0': number of b=0 images in vol
%     {default=1}
%
%
% Created:  02/22/10 Don Hagler
% Last Mod: 10/27/12 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse arguments
if ~mmil_check_nargs(nargin, 3), return; end;
parms = mmil_args2parms( varargin, { ...
  'nb0',1,[0 Inf],...
  'volmask',[],[],...
  'censor_min_ndirs',6,[],...
});

volsz = size(vol);

if any(volsz~=size(vol_synth))
  error('dimensions of vol and vol_synth do not match');
end;

nz = size(vol,3);
nf = size(vol,4);

if any([nz,nf]~=size(censor_mat))
  error('size of censor_mat must be [%d,%d]',nz,nf);
end;

% how many directions in tensor fit?
nframes = size(vol,4) - sum(censor_mat,2);
censor_mat(nframes<parms.censor_min_ndirs + parms.nb0,:) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for z=1:nz
  for f=1:nf
    if censor_mat(z,f)
      if isempty(parms.volmask)
        vol(:,:,z,f) = vol_synth(:,:,z,f);
      else
        im_orig = squeeze(vol(:,:,z,f));
        im_synth = squeeze(vol_synth(:,:,z,f));
        im_mask = squeeze(parms.volmask(:,:,z));
        im_orig(im_mask>0) = im_synth(im_mask>0);
        vol(:,:,z,f) = im_orig;
      end;      
    end;
  end;
end;

