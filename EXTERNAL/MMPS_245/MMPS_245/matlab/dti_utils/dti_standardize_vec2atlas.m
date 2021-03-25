function vec = dti_standardize_vec2atlas(vec,M,varargin)
%function vec = dti_standardize_vec2atlas(vec,M,[options])
%
% Purpose: standardize matrix of diffusion vectors
%   to match atlas (LPI), according to M matrix
%
%
% Required Input:
%   vec: matrix of [x y z] vectors (should be n x 3)
%   M: vox2ras matrix
%     will flip x, y, and z of vecs to match LPI, 
%     depending on orientation of M
%
% Optional Input:
%  'orient_ref': reference orientation
%    {default = 'LPI'}
%
% Created:  03/10/10 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'orient_ref','LPI',[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vec_orig = vec;
if all(size(vec_orig)==[3 1])
  vec = reshape(vec_orig,[1 3]);
elseif size(vec_orig,2)~=3
  error('vec should be n x 3');
end;

% flip or reorder vec to match orient_ref based on M
orient = fs_read_orient([],M);
[permvec,flipvec] = fs_compare_orient(orient,parms.orient_ref);
vec = vec(:,permvec);
if any(flipvec>0)
  for i=1:3
    if flipvec(i)<0
      vec(:,i) = -vec(:,i);
    end;
  end;
end;

% in case input was 3 x 1
vec = reshape(vec,size(vec_orig));

