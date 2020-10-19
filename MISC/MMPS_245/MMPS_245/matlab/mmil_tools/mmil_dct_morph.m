function regStruct = mmil_dct_morph(vol,volm,varargin)
%function regStruct = mmil_dct_morph(vol,volm,[options])
%
%  Required Input:
%   vol: input brain volume (ctx format)
%   volm: target brain volume (ctx format)
%
%  Optional Input:
%  'volstd': target brain std volume (ctx format)
%     {default = []}
%  'volmask': target brain mask (ctx format)
%     {default  = []}
%  'stdbgval': value of of std outside brain mask
%     {default = 75}
%  'smoothflag': [0|1] apply Gaussian smoothing to input volume
%     {default = 1}
%  'sampling': sampling rate (num voxels) in each dimension
%     {default = [4 4 4]}
%  'nK': number of DCT basis functions in each dimension
%     {default = [5 5 5]}
%  'tstep': size of translation step (mm)
%     {default = 0.5}
%  'tstep2': size of translation step (mm) for second iteration
%     {default = 0.25}
%  'astep': size of angle step (degrees)
%     {default = 0.25}
%  'astep2': size of angle step (degrees) for second iteration
%     {default = 0.5}
%  'scales': vector of scales for multi-scale search
%     {default = [0 83 49 27 16 9 5 3 2 1]}
%  'ns': number of samples
%     {default = 64}
%  'ns2': number of samples for second iteration
%     {default = 32}
%  'sf': scaling factor
%     {default = 1}
%  'thresh': threshold applied to vol
%     to prevent zero values affecting registration
%     {default = 0}
%  'verbose': [0|1] whether to display progress messages
%     {default = 1}
%
% Output:
%   regStruct: structure containing fields:
%     M_atl_to_vol_rb: rigid body registration matrix
%     M_atl_to_vol_af: affine registration matrix
%     VL: volume containing displacements in L-R axis
%     VP: volume containing displacements in A-P axis
%     VH: volume containing displacements in I-S axis
%     Tr: output from nonlinear registration
%     sf_rb: scaling factor for rigid body registration
%     bconverged
%     bcombined
%     min_cost_rb
%     min_cost_af
%     min_cost_m
%
% Created:  04/02/12 by Don Hagler
% Last Mod: 05/13/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'volstd',[],[],...
  'volmask',[],[],...
  'stdbgval',75,[],...
  'smoothflag',true,[false true],...
  'sampling',[4 4 4],[],...
  'nK',[5 5 5],[],...
  'tstep',0.5,[],...
  'astep',0.25,[],...
  'scales',[0 83 49 27 16 9 5 3 2 1],[],...
  'ns',64,[],...
  'sf',1,[],...
  'thresh',0,[0,Inf],...
  'verbose',true,[false true],...
  'ns2',32,[],...
  'tstep2',0.25,[],...
  'astep2',0.5,[],...
});
regStruct = [];

% convert astep to radians
parms.astep = parms.astep*(pi/180);

volm.imgs = parms.sf * volm.imgs;

if isempty(parms.volstd)
  parms.volstd = volm;
  parms.volstd.imgs = ones(size(volm.imgs));
end

if isempty(parms.volmask)
  parms.volmask = volm;
else
  ind_nonbrain = find(parms.volmask.imgs<=eps);
  parms.volmask = volm;
  parms.volmask.imgs(ind_nonbrain) = 0;
  parms.volstd.imgs(ind_nonbrain) = parms.stdbgval;
end;

parms.volstd.imgs = parms.sf*parms.volstd.imgs;

% resample to orientation and resolution of atlas
%% todo: more specific test for whether resampling is needed
if vol.dimr < (0.75*volm.dimr)
  fprintf('%s: resampling input volume...\n',mfilename);
  vol = vol_resample_pad(vol,volm,eye(4));
end

% Gaussian smoothing
if parms.smoothflag
  fprintf('%s: smoothing input volume...\n',mfilename);
  vol = vol_filter(vol,1);
end;
bsmooth = false;

if ~isempty(parms.scales) % skip linear reg, if "scales" empty
  if parms.thresh>0
    if parms.verbose
     fprintf('%s: rigid body registration by multiscale search with thresh = %0.2g...\n',...
      mfilename,parms.thresh);
     tic;
    end
    [regStruct.M_atl_to_vol_rb, regStruct.min_cost_rb, regStruct.sf_rb] = ...
        reg_vol2hatl_rb_mmil(vol,parms.volmask,parms.volstd,parms.ns,false,...
           bsmooth,parms.scales,parms.tstep,parms.astep,2,eye(4,4),parms.thresh);
    % second iteration in case first registration failed
    % VJ NOTE: does not provide major improvement
    %%         but included for backward compatibility
    if (regStruct.min_cost_rb>8)
      if (parms.verbose)
        fprintf('%s:   alignment may have failed (min cost = %0.1f)\n',...
          mfilename,regStruct.min_cost_rb);
        fprintf('%s:   realigning with smaller steps, larger angular steps and more points...\n',...
          mfilename);
      end
      % convert astep2 to radians
      parms.astep2 = parms.astep2*(pi/180);
      [M_atl_to_vol_rb, min_cost_rb, sf_rb] = ...
        reg_vol2hatl_rb_mmil(vol,parms.volmask,parms.volstd, parms.ns2, false, bsmooth, ...
        parms.scales, parms.tstep2, parms.astep2, 2, regStruct.M_atl_to_vol_rb,parms.thresh);
      if min_cost_rb < regStruct.min_cost_rb
        regStruct.M_atl_to_vol_rb = M_atl_to_vol_rb;
        regStruct.min_cost_rb = min_cost_rb;
        regStruct.sf_rb =sf_rb;
      end
    end
  else
    if parms.verbose
     fprintf('%s: rigid body registration by multiscale search...\n',mfilename);
     tic;
    end
    [regStruct.M_atl_to_vol_rb, regStruct.min_cost_rb, regStruct.sf_rb] = ...
        reg_vol2hatl_rb(vol,parms.volmask,parms.volstd,parms.ns,false,...
           bsmooth,parms.scales,parms.tstep,parms.astep,2,eye(4,4));
  end;
  if parms.verbose
    toc;
    fprintf('%s:   min cost = %0.1f\n',...
      mfilename,regStruct.min_cost_rb);
    fprintf('%s:   scaling factor = %0.2f\n',...
      mfilename,regStruct.sf_rb);
  end;
  vol.imgs = vol.imgs.*regStruct.sf_rb;

  if parms.verbose
    fprintf('%s: affine registration by simplex...\n',mfilename);
    tic;
  end
  errfn = 'calmi_hatl';
  maxiter= 2000;
  ftol=0.001;
  mindp=0.01;
  initp=[0 0 0 1 1 1 0 0 0 0 0 0 ];
  wei=[0.01 0.01 0.01 0.2 0.2 0.2 0.01 0.01 0.01 0.01 0.01 0.01];
  [regStruct.M_atl_to_vol_af, regStruct.min_cost_af, regStruct.bconverged] = ...
      reg_vol2hatl_af(vol, volm, parms.volstd, parms.sampling,...
         regStruct.M_atl_to_vol_rb, bsmooth, ...
         errfn, initp, wei, maxiter, ftol, mindp);
  if parms.verbose
    toc;
    fprintf('%s:   min cost = %0.1f\n',...
      mfilename,regStruct.min_cost_af);
  end;
else
  regStruct.M_atl_to_vol_rb = eye(4);
  regStruct.M_atl_to_vol_af = eye(4);
end

if parms.verbose
  fprintf('%s: nonlinear registration by DCT...\n',mfilename);
  tic;
end
titer=16;
regStruct.bcombined = true;
regBend=1.0;
regMem=0;
stablize=1.;
[regStruct.VL regStruct.VP regStruct.VH regStruct.Tr, regStruct.min_cost_m] = ...
    dctMorph_vol2hatl_amd(vol, volm, parms.volstd, regStruct.M_atl_to_vol_af,...
      parms.nK, parms.sampling, false, titer, regStruct.bcombined,...
      regBend, regMem, stablize, bsmooth, parms.verbose);
if parms.verbose
  fprintf('\n');
  toc;
end;
