function vol = dti_synth_vol(DTfit)
%function vol = dti_synth_vol(DTfit)
%
% Required Parameters:
%   DTfit: output structure from dti_fit_tensor with fields:
%     volDT : tensor fit parameters    size = [nx,ny,nz,7]
%     volsz : size of input 4D vol          = [nx,ny,nz,nf]
%     Q : tensor fit forward matrix    size = [nf,7]
%
% Created:  02/18/10 Don Hagler
% Last Mod: 10/27/12 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse arguments
if ~mmil_check_nargs(nargin, 1), return; end;

nx = DTfit.volsz(1);
ny = DTfit.volsz(2);
nz = DTfit.volsz(3);
nf = DTfit.volsz(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fprintf('%s: synthesizing volume from tensor fit...\n',mfilename);
vol = zeros(DTfit.volsz,'single');
for z = 1:nz
  T = reshape(DTfit.volDT(:,:,z,:),[nx*ny,7])';
  tmp_mat = exp(DTfit.Q*T);
  vol(:,:,z,:) = single(reshape(tmp_mat',[nx ny nf]));
end;

