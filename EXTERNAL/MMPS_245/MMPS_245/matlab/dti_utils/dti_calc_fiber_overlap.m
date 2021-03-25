function overlap = dti_calc_fiber_overlap(fname_mask,fname_prob,...
  thresh_list,fname_FA,thresh_FA)
%function overlap = dti_calc_fiber_overlap(fname_mask,fname_prob,...
%  [thresh_list],[fname_FA],[thresh_FA])
%
% Created:  10/18/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

overlap = [];

if ~mmil_check_nargs(nargin,2), return; end;

if ~exist('thresh_list','var') | isempty(thresh_list)
  thresh_list = 0.5;
end;
if ~exist('fname_FA','var')
  fname_FA = [];
end;
if ~exist('thresh_FA','var') | isempty(thresh_FA)
  thresh_FA = 0.2;
end;

if ~exist(fname_mask,'file'), error('file %s not found',fname_mask); end;
if ~exist(fname_prob,'file'), error('file %s not found',fname_prob); end;

vol_mask = fs_load_mgh(fname_mask);
vol_prob = fs_load_mgh(fname_prob);

if length(size(vol_mask))~=length(size(vol_prob)) |...
   any(size(vol_mask)~=size(vol_prob))
  error('mismatch in input file dimensions');
end;

if ~isempty(fname_FA)
  if ~exist(fname_FA,'file'), error('file %s not found',fname_FA); end;
  vol_FA = fs_load_mgh(fname_FA);
  if length(size(vol_mask))~=length(size(vol_FA)) |...
     any(size(vol_mask)~=size(vol_FA))
    error('mismatch in input FA file dimensions');
  end;
  vol_FA(vol_FA<thresh_FA) = 0;
  vol_FA(vol_FA>=thresh_FA) = 1;
  vol_prob = vol_prob.*vol_FA;
end;

vol_mask(vol_mask>0)=1;

overlap = zeros(size(thresh_list));
for i=1:length(thresh_list)
  vol_prob_mask = zeros(size(vol_prob));
  vol_prob_mask(vol_prob>thresh_list(i))=1;
  overlap(i) = calc_volume_overlap(vol_mask,vol_prob_mask);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function overlap = calc_volume_overlap(volmask1,volmask2)

volintersect = volmask1+volmask2;

volintersect(volintersect<2)=0;
n1 = length(find(volmask1));
n2 = length(find(volmask2));
ni = length(find(volintersect));
if (n1+n2)>0
  overlap = 2*ni/(n1+n2);
end;

return;
