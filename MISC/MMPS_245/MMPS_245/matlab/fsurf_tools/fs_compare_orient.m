function [permvec,flipvec] = fs_compare_orient(orient,orient_ref)
%function [permvec,flipvec] = fs_compare_orient(orient,orient_ref)
%
% Purpose: compare orientation of a volume to a reference
%   to determine permutation and dimension flipping required
%   to reorient one volume into same orientation as a reference
%
% Required Input:
%   orient: three character string specifying volume orientation
%     e.g. 'RAS', 'LPI', 'PRI', etc.
%   orient_ref: reference orientation
%
% Output:
%   permvec: permutation vector needed to reorient volume to be
%     same as reference
%   flipvec: vector of -1 or 1 for each dimension (after permutation)
%
% Created:  03/10/10 by Don Hagler
% Last Mod: 03/10/10 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

permvec = [];
flipvec = [];
if (~mmil_check_nargs(nargin,2)) return; end;

% check for valid orient
orient = upper(orient);
if any(~ismember(orient(:),{'R','L','A','P','S','I'})) ||...
   length(unique(orient(:)))~=3 ||...
   isempty(regexp(orient,'[LR]')) ||...
   isempty(regexp(orient,'[AP]')) ||...
   isempty(regexp(orient,'[SI]'))
  error('invalid orient %s',orient);
end;

% check for valid orient_ref
orient_ref = upper(orient_ref);
if any(~ismember(orient_ref(:),{'R','L','A','P','S','I'})) ||...
   length(unique(orient_ref(:)))~=3 ||...
   isempty(regexp(orient_ref,'[LR]')) ||...
   isempty(regexp(orient_ref,'[AP]')) ||...
   isempty(regexp(orient_ref,'[SI]'))
  error('invalid orient_ref %s',orient_ref);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine how to permute and flip dims
in_order = zeros(1,3);
in_dirs = zeros(1,3);
out_order = zeros(1,3);
out_dirs = zeros(1,3);
permvec = zeros(1,3);
flipvec = zeros(1,3);

j = regexp(orient,'R');
k = regexp(orient,'L');
if ~isempty(j)
  in_order(j) = 1;
  in_dirs(j) = 1;
else
  in_order(k) = 1;
  in_dirs(k) = -1;
end;

j = regexp(orient,'A');
k = regexp(orient,'P');
if ~isempty(j)
  in_order(j) = 2;
  in_dirs(j) = 1;
else
  in_order(k) = 2;
  in_dirs(k) = -1;
end;

j = regexp(orient,'S');
k = regexp(orient,'I');
if ~isempty(j)
  in_order(j) = 3;
  in_dirs(j) = 1;
else
  in_order(k) = 3;
  in_dirs(k) = -1;
end;

j = regexp(orient_ref,'R');
k = regexp(orient_ref,'L');
if ~isempty(j)
  out_order(j) = 1;
  out_dirs(j) = 1;
else
  out_order(k) = 1;
  out_dirs(k) = -1;
end;

j = regexp(orient_ref,'A');
k = regexp(orient_ref,'P');
if ~isempty(j)
  out_order(j) = 2;
  out_dirs(j) = 1;
else
  out_order(k) = 2;
  out_dirs(k) = -1;
end;

j = regexp(orient_ref,'S');
k = regexp(orient_ref,'I');
if ~isempty(j)
  out_order(j) = 3;
  out_dirs(j) = 1;
else
  out_order(k) = 3;
  out_dirs(k) = -1;
end;

for i=1:3
  j = out_order(i);
  k = find(in_order == j);
  permvec(i) = k;
  flipvec(i) = in_dirs(k)*out_dirs(i);
end;

