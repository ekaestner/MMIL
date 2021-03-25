function [sum_log10_val, val, bins1, bins2, jentropy]=vols_jhist_mask_amd(vol1, vol2, vol1mask, varargin)
% Calculate Joint histogram for two image volume based on vol1's mask
%   
%   [sum_log10_val, val, bins1, bins2, jentropy]=vols_jhist(vol1, vol2, vol1mask, [mstep], [Mreg1], [Mreg2], [numbins1], [numbins2], [thmin1], [thmax1], [thmin2], [thmax2], [interpm], [padding])
%  
% Input:
%       vol*:       Input Volume
%       vol1mask:   Mask Volume;
%       mstep:      sampling rate in volmask (default =1)     
%       Mreg*:       Linear Transformation matrix 4X4 (defualt=I) LPH based
%       numbins*:   Number of bins
%       thmin*:      Threshold for min val 
%       thmax*:      Threshold for max val
%       interpm :   0: : Nearest 1:Linear (default) 2:cubic 
%                   3: Key's spline 4: Cubic spline. 5: Hamming_Sinc
%       padding :   half width of number of points used in the interpolation,
%                   only for 4 and 5. For example, 3 means it would use
%                    6*6*6 points for interpolation.
%  Output:
%       sum_log10_val: sum of log10 of (val+0.5)
%       val:        counts
%       bins*:      bins
%       entropy:    joint entropy

mstep =1;
if nargin>=4
    mstep=varargin{1};
end

Mreg1= eye(4,4);
if nargin >= 5
  Mreg1 = varargin{2};
end

Mreg2= eye(4,4);
if nargin >= 6
  Mreg2 = varargin{3};
end

if isempty(vol1mask)
  vol1mask = vol1;
  vol1mask.imgs(:) = 1;
end

ind = find(vol1mask.imgs>0);
vsize = length(ind);
numbins1 = ceil(vsize^(1/3));
if nargin>=7
    numbins1=varargin{4};
end

numbins2 = ceil(vsize^(1/3));
if nargin>=8
    numbins2=varargin{5};
end

bmin1=false;
thmin1=0;
if nargin>=9
    bmin1=true;
    thmin1=varargin{6};
end

bmax1=false;
thmax1=0;
if nargin>=10
    bmax1=true;
    thmax1=varargin{7};
end


bmin2=false;
thmin2=0;
if nargin>=11
    bmin2=true;
    thmin2=varargin{8};
end

bmax2=false;
thmax2=0;
if nargin>=12
    bmax2=true;
    thmax2=varargin{9};
end

interpm=1;
if nargin>=13
    interpm=varargin{10};
end

padding=0;
if nargin >= 14
  padding = varargin{11};
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


Mlpho2vxli1=inv(vol1.Mvxl2lph)*Mreg1;
Mlpho2vxli2=inv(vol2.Mvxl2lph)*Mreg2;

[sum_log10_val, val, bins1, bins2, jentropy]=volsjhistMEX(size(vol1.imgs), vol1.imgs, size(vol2.imgs), vol2.imgs, ...
                                        numbins1, numbins2, [], Mlpho2vxli1, Mlpho2vxli2, ...
                                        bmin1, thmin1, bmax1, thmax1, bmin2, thmin2, bmax2, thmax2, interpm, padding, vol1mask.imgs, vsize, vol1mask.Mvxl2lph, size(vol1mask.imgs), mstep);
