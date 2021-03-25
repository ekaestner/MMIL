function mi=vols_getMI_mask_djh(vol1, vol2, volmask,varargin)
% Calculate Mutual Infomation for two image volume based on masked volume
%   
%   mi=vols_getMI_mask_djh(vol1, vol2, volmask, [Mreg1], [Mreg2], [mstep], [thmin1], [thmax1], [thmin2], [thmax2], [interpm], [padding])
%  
% Input:
%       vol*:        Input Volume
%       volmask:     Mask Volume
%       Mreg*:       Linear Transformation matrix 4X4 (defualt=I) LPH based
%       mstep:       Samping rate in masked volume (default =1)
%       thmin*:      Threshold for min val 
%       thmax*:      Threshold for max val
%       interpm :   0: : Nearest 1:Linear(default) 2:cubic 
%                   3: Key's spline 4: Cubic spline. 5: Hamming_Sinc
%       padding :   half width of number of points used in the interpolation,
%                   only for 4 and 5. For example, 3 means it would use
%                    6*6*6 points for interpolation.
%  Output:
%      mi:      Mutual Information


Mreg1= eye(4,4);
if nargin >= 4
  Mreg1 = varargin{1};
end

Mreg2= eye(4,4);
if nargin >= 5
  Mreg2 = varargin{2};
end

mstep =1;
if nargin>=6
    mstep=varargin{3};
end

bmin1=false;
thmin1=0;
if nargin>=7
    bmin1=true;
    thmin1=varargin{4};
end

bmax1=false;
thmax1=0;
if nargin>=8
    bmax1=true;
    thmax1=varargin{5};
end


bmin2=false;
thmin2=0;
if nargin>=9
    bmin2=true;
    thmin2=varargin{6};
end

bmax2=false;
thmax2=0;
if nargin>=10
    bmax2=true;
    thmax2=varargin{7};
end

interpm=1;
if nargin>=11
    interpm=varargin{8};
end

padding=0;
if nargin >= 12
  padding = varargin{9};
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

ind = find(volmask.imgs>0);
vsize = length(ind);
numbins1 = ceil(vsize^(1/3));
numbins2 = numbins1;
Mlpho2vxli1=inv(vol1.Mvxl2lph)*Mreg1;
Mlpho2vxli2=inv(vol2.Mvxl2lph)*Mreg2;

[val1, bins1, entropy1]=volhistMEX(size(vol1.imgs), vol1.imgs, numbins1, [], Mlpho2vxli1, bmin1, thmin1, bmax1, thmax1,interpm, padding, volmask.imgs, vsize, volmask.Mvxl2lph, size(volmask.imgs), mstep);
[val2, bins2, entropy2]=volhistMEX(size(vol2.imgs), vol2.imgs, numbins2, [], Mlpho2vxli2, bmin2, thmin2, bmax2, thmax2,interpm, padding, volmask.imgs, vsize, volmask.Mvxl2lph, size(volmask.imgs), mstep);
[sum_log10_val, jval, jbins1, jbins2, jentropy]=volsjhistMEX(size(vol1.imgs), vol1.imgs, size(vol2.imgs), vol2.imgs, ...
                                        numbins1, numbins2, [], Mlpho2vxli1, Mlpho2vxli2, ...
                                        bmin1, thmin1, bmax1, thmax1, bmin2, thmin2, bmax2, thmax2,interpm, padding, volmask.imgs, vsize, volmask.Mvxl2lph, size(volmask.imgs), mstep);
mi=entropy1+entropy2-jentropy;


