function [VL VP VH Tr Var] = dctMorph_vol2hatl_amd(vol, volm, volstd, ...
                                               M_atl_to_vol_af, nK, ...
                                               sampling, debug, ...
                                               titer, bcombined, ...
                                               regBend, regMem, ...
                                               stablize, bsmooth, boutput)
%
% dctM = dctMorph_vol2hatl(vol, volm, volstd, M_volt_to_vol_af, nK, sampling, ...
%                          [debug, titer,bcombined, regBend, regMem, stablize, bsmooth])
%
% Input:
%       vol: vol structure
%       volm: atlas mean volume
%       volstd: atlas std volume
%       M_atl_to_vol_af: Affine transfrom from atlas to vol in lph space
%       nK (1X3): nubmber of DCT bases functions in each directions
%       sampling(1X3): sample size for masking voxles
%       bcombined: Combined vector field with Affine transfrom (default=1) 
%       titer: totol iteration (defulat = 16);
%       regBend: regularisation for bending Energy (default = 1)
%       regMem:  regularisation for membrance Energy (default = 0)
%       debug: show morph volume (default = 0)
%       stablize: stablize factor for Inverse (default = 1)
%       bsmooth: smooth volume (default = true)
%
%
% Output:
%   VL, VP, VH: Vector field in LPH space
%   Tr: DCT coeffients.
%   min_cost: MI
%
% Created:  06/23/07 by Anders Dale
% Last Mod: 06/23/07 by Anders Dale
%


% Refs:
% 
% 1. J. Ashburner and K. J. Friston
% Nonlinear Spatial Normalization using Basis Functions.
% Human Brain Mapping 7(4):in press (1999)
%
% Modified by AMD to use less memory
%

if (~exist('bcombined'))
  bcombined=1;
end

if (~exist('regBend'))
  regBend=1;
end

if (~exist('regMem'))
  regMem=0;
end

if (~exist('debug'))
  debug=0;
end

if (~exist('stablize'))
  stablize=1;
end

if (~exist('titer'))
  titer=16;
end

if (~exist('bsmooth'))
  bsmooth=true;
end

if (~exist('boutput'))
  boutput=true;
end

if boutput, fprintf('\n     Resampling vol to atlas and Smoothing...'); end
volr = vol_resample(vol, volm, M_atl_to_vol_af);
clear vol
if (bsmooth)
  volr = vol_filter(volr,1);
end
volr.imgs = single(volr.imgs);

if boutput, fprintf('\n     Calculate Gradient...'); end
volgi = vol_filter(volr,5);
volgj = vol_filter(volr,4);
volgk = vol_filter(volr,6);

volr.imgs = double(volr.imgs);
volgi.imgs = double(volgi.imgs);
volgj.imgs = double(volgj.imgs);
volgk.imgs = double(volgk.imgs);

if boutput, fprintf('\n     Precalculate DCT functions and its derivatives...'); end
dim = size(volm.imgs);
bi = stablize*dctb(dim(1), nK(1), (0:dim(1)-1)',0);
bj = stablize*dctb(dim(2), nK(2), (0:dim(2)-1)',0);
bk = stablize*dctb(dim(3), nK(3), (0:dim(3)-1)',0);
bi1 = stablize*dctb(dim(1), nK(1), (0:dim(1)-1)',1);
bj1 = stablize*dctb(dim(2), nK(2), (0:dim(2)-1)',1);
bk1 = stablize*dctb(dim(3), nK(3), (0:dim(3)-1)',1);
bi2 = stablize*dctb(dim(1), nK(1), (0:dim(1)-1)',2);
bj2 = stablize*dctb(dim(2), nK(2), (0:dim(2)-1)',2);
bk2 = stablize*dctb(dim(3), nK(3), (0:dim(3)-1)',2);

if boutput, fprintf('\n     Calculate Regulization Term...'); end
R = getRegularization(bi, bj, bk, bj1, bj1, bk1, bj2, bj2, bk2, regBend, regMem);

if boutput, fprintf('\n     Levenberg-Marguardt Minimization...'); end
edgeskip=sampling;

pVar = Inf;

s1 = 3*prod(nK);
T  = zeros(s1,1);

for iter=1:titer
  if boutput, fprintf('\n         iteration %2d: ', iter); end
  [Alpha,Beta, Var] = getABMEX(dim, volr.imgs, volgi.imgs, volgj.imgs, ...
                                    volgk.imgs, volm.imgs, volstd.imgs, ...
                                    bi, bj, bk, bi1, bj1, bk1, T, ...
                                    sampling, edgeskip);
  if boutput, fprintf(' Var %f: ', Var); end
  if Var>pVar
    scal = pVar/Var ; 
    Var = pVar; 
  else
    scal = 1;
  end;
  pVar = Var;
  T = (Alpha + R*scal)\(Alpha*T + Beta);
end;

Tr = reshape(T,[nK 3]);
Tr = Tr*stablize.^3;
bi = bi./stablize;
bj = bj./stablize;
bk = bk./stablize;

clear volstd volgi volgj volgk

if boutput, fprintf('\n     Calculate Vector field in LPH coordinates...'); end

if (bcombined)
  M = (M_atl_to_vol_af)*volm.Mvxl2lph;
else
  M = volm.Mvxl2lph;
end
i=1:dim(1);
j=1:dim(2);
k=1:dim(3);
[I,J] = ndgrid(i, j);
VL=volm;
VL.imgs=zeros(size(volm.imgs));
VP=VL;
VH=VL;
lph =ones(dim(1)*dim(2),4);

if (debug)
  volrm = volm;
  volrm.imgs=zeros(dim);
end

for l=k,   % Cycle over planes
           % Nonlinear deformations
           %----------------------------------------------------------------------------
           ti = get_2Dtrans(Tr(:,:,:,1),bk,l);
           tj = get_2Dtrans(Tr(:,:,:,2),bk,l);
           tk = get_2Dtrans(Tr(:,:,:,3),bk,l);
           vj = J    + bi*tj*bj';
           vi = I    + bi*ti*bj';
           vk = k(l) + bi*tk*bj';

           [vl vp vh] = mmult(vi, vj, vk, M);
           [vlo vpo vho] = mmult(I, J, k(l), volm.Mvxl2lph);
           
           
           VL.imgs(:,:,l)=vl-vlo;
           VP.imgs(:,:,l)=vp-vpo;
           VH.imgs(:,:,l)=vh-vho;
end      

[VL.maxI VL.minI]=maxmin(VL.imgs);
[VP.maxI VP.minI]=maxmin(VP.imgs);
[VH.maxI VH.minI]=maxmin(VH.imgs);

