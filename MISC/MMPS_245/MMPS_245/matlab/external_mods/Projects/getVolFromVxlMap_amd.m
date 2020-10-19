function volm = getVolFromVxlMap_amd(vol, vxlmap, btype, varargin)
%
% Morph Volume to and from subject space to atlas space
%
%   volm = getVolFromVxlMap(vol, vxlmap, btype, [interpm], [padding], [bclamp])
%
%   btype: 0: mapping from subject space to rgeistration atlas space
%             vol:  volume in subject space 
%             volm: volume in regsitration atlas space
%
%   btype: 1: mapping from regsitration atlas space to subject space
%             vol:  volume in regstration atlas space 
%             volm: volume in subject space space
%
%   interpm :     0: Nearest Neighbor 1:Linear  2:cubic (default)
%                 3: Key's spline 4: Cubic spline. 5: Hamming_Sinc
%
%   padding :     half width of number of points used in the interpolation,
%                 only for Cubic spline and Hamming_Sinc. 
%                 For example, 3 means it would use
%                 6*6*6 points for interpolation.
%   bclamp : set negative value to zero default = true
%
%
% Created:  ?        by Anders Dale
% Last Mod: 07/09/21 by Don Hagler
%


interpm = 2;
if nargin >= 4
  interpm = varargin{1};
end


padding=3;
if nargin >= 5
  padding = varargin{2};
end

if (interpm==0) % Nearest Neighbor padding=0
    padding=0;
elseif (interpm==1) % Lineae  padding=1
    padding=1;
elseif (interpm==2) % cubic and key's padding =2
    padding=2;
elseif (interpm==3)
    padding=2;
end

bclamp=true;
if nargin >= 6
  bclamp = varargin{3};
end

if (btype==0) 
  volm = vxlmap.atl2subj.indVol;
  volm.imgs = zeros(size(volm.imgs));
  ind = find(vxlmap.atl2subj.indVol.imgs>0);
  indx = double(vxlmap.atl2subj.indVol.imgs(ind)); 
  vxl =double([vxlmap.atl2subj.I(indx) vxlmap.atl2subj.J(indx) vxlmap.atl2subj.K(indx) ones(size(vxlmap.atl2subj.I(indx)))]);
else
  volm = vxlmap.subj2atl.indVol;
  volm.imgs = zeros(size(volm.imgs));
  ind = find(vxlmap.subj2atl.indVol.imgs>0);
  indx = double(vxlmap.subj2atl.indVol.imgs(ind)); 
  vxl =double([vxlmap.subj2atl.I(indx) vxlmap.subj2atl.J(indx) vxlmap.subj2atl.K(indx) ones(size(vxlmap.subj2atl.I(indx)))]);
end

clear indx vxlmap;

lph = (vol.Mvxl2lph)*vxl';
vxlval = vol_getvxlsval(lph', vol, eye(4,4), interpm, padding);
volm.imgs(ind) = vxlval;
if (bclamp)
  ind0=find(volm.imgs<0);
  volm.imgs(ind0)=0;
end
 
volm.maxI = vol.maxI;
volm.minI = vol.minI;
