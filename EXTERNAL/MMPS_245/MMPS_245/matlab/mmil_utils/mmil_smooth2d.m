function im_out = mmil_smooth2d(im_in,sig1,sig2)
%function im_out = mmil_smooth2d(im_in,sig1,sig2)
%
% Created:   03/09/10 by Don Hagler based on code from Anders Dale
% Last Mod:  05/21/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

orig_type = class(im_in);
if ~strcmp(orig_type,'single')
  im_in = single(im_in);
end;

idim = size(im_in,1);
jdim = size(im_in,2);
filtim = exp(-1/2*(((([1:idim]-(idim/2+1))'*ones(1,jdim))/idim*sig1).^2+((ones(idim,1)*([1:jdim]-(jdim/2+1)))/jdim*sig2).^2));
kim = fftshift(ifftn(im_in));
im_out = fftn(ifftshift(kim.*filtim));

im_out = real(im_out);
if strcmp(orig_type,'double')
  im_out = double(im_out);
end;


