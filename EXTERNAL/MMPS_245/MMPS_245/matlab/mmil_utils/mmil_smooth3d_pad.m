function vol_out = mmil_smooth3d_pad(vol_in,sig)
%function vol_out = mmil_smooth3d_pad(vol_in,sig)
%
% Purpose: smooth 3D volume with zero-padding on edges
%   to prevent wrap-around effect
%
% Created:   01/20/15 by Don Hagler
% Last Mod:  01/20/15 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

if numel(sig)==1
  sig = sig*ones(3,1);
elseif numel(sig)~=3
  error('sig must have 1 or 3 elements');
end;
np = max(sig);

[nx,ny,nz] = size(vol_in);
nxp = nx + 2*np;
nyp = ny + 2*np;
nzp = nz + 2*np;
% pad vol with zeros on all sides
vol_pad = zeros(nxp,nyp,nzp);
vol_pad(np:np+nx-1,np:np+ny-1,np:np+nz-1) = vol_in;
% smooth padded image
vol_out = mmil_smooth3d(vol_pad,sig(1),sig(2),sig(3));
% remove padding
vol_out = vol_out(np:np+nx-1,np:np+ny-1,np:np+nz-1);


