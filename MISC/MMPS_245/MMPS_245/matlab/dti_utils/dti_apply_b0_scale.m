function vol = dti_apply_b0_scale(vol,i_b0,scalefacts)
%function vol = dti_apply_b0_scale(vol,i_b0,scalefacts)
%
% Purpose: scale diffusion weighted images according to scale factors
%   calculated from preceding b=0 images
%
% Required Parameters:
%   vol: 4D volume containing diffusion MRI data
%   i_b0: vector of frame numbers corresponding to b=0 images
%   scalefacts: vector of scale factors, one for each i_b0
%
% Output:
%    vol : 4D volume with scaled intensities
%
% Created:  05/08/12 Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

nf = size(vol,4);
nb0 = length(i_b0);
if nb0<=1, return; end;
if max(i_b0)>nf
  error('i_b0 contains value greater than number of frames (%d)',nf);
end;
if min(i_b0)<1
  error('i_b0 contains value less than one');
end;

scalefacts_all = ones(nf,1);
b = 1;
f_bnext = i_b0(b+1);
for f=1:nf
  if f>=f_bnext
    if b+1>=length(i_b0)
      f_bnext = nf;
      b = length(i_b0);
    else
      b = b+1;
      f_bnext = i_b0(b+1);
    end;
  end;
  scalefacts_all(f) = scalefacts(b);
end;
% apply scalefacts
for f=1:nf
  vol(:,:,:,f) = vol(:,:,:,f)*scalefacts_all(f);
end;


return;

