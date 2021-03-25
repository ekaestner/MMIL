function vec_out = mmil_smooth1d(vec_in,sig1)
%function vec_out = mmil_smooth1d(vec_in,sig1)
%
% Created:   03/09/10 by Don Hagler based on code from Anders Dale
% Last Mod:  05/21/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

orig_type = class(vec_in);
if ~strcmp(orig_type,'single')
  vec_in = single(vec_in);
end;

idim = length(vec_in);
filtvec = exp(-1/2*(((([1:idim]-(idim/2+1)))/idim*sig1).^2));
if (size(filtvec,1) ~= size(vec_in,1))
  filtvec = filtvec.';
end
kvec = fftshift(fft(vec_in));
vec_out = ifft(ifftshift(kvec.*filtvec));

vec_out = real(vec_out);
if strcmp(orig_type,'double')
  vec_out = double(vec_out);
end;


