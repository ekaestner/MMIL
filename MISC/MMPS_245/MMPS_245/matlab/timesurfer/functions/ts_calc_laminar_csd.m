function [data_csd,data_ecd,data_prep] = ts_calc_laminar_csd(data,varargin)
%function [data_csd,data_ecd,data_prep] = ts_calc_laminar_csd(data,[options])
%
% Purpose: calculate current source density from laminar data
%
% Required Input:
%   data: 2D data matrix, with size = [nchans,ntpoints]
%     units of data should be mV
%     assumed to be potential gradient
%
% Optional Parameters for data preparation:
%   'badchans': vector of bad channels
%     will be replaced with weighted average of surrounding good channels
%     {default = []}
%   'chan_wf': exponential decay factor (0 to 1)
%     used to weight neighbors of bad channels based on relative distance
%     {default = 0.1}
%   'smoothing': smoothing sigma (# of channels)
%     {default = 0}
%
% Optional Parameters that control method of CSD calculation:
%   'hamming_flag': [0|1] use 5-point hamming filter
%     otherwise perform calculation by inverse
%     {default = 0}
%   'edge_corr_flag': [0|1] correct for surface potentials
%     applies to both inverse and hamming methods (but done differently)
%     {default = 1}
%
% Optional Parameters for CSD calculation by inverse:
%   'SNR': assumed signal to noise ratio along the shaft of the cylinder
%     {default = 20}
%   'dz': intrachannel spacing (mm)
%     {default = 0.15}
%   'd1': inner diameter of cylinder (mm)
%     {default = 0}
%   'd2': outer diameter of cylinder (mm)
%     {default = 10}
%   'sigma': conductivity (S/m or mS/mm) of tissue
%     {default = 0.3}
%
% Optional Parameters for CSD calculation by hamming filter:
%   'r': resistance in ohm/mm
%     {default = 20}
%
% Output:
%   data_csd: current source density
%     matrix with size = [nchans+1,ntpoints]
%     units = uA/mm^3
%   data_ecd: equivalent current dipole
%     vector with size = [1,ntpoints]
%     units = uA*mm/mm^2
%   data_prep: smoothed data with badchans censored
%     matrix with size = [nchans,ntpoints]
%     units = mV
%
% Created:  08/12/14 by Don Hagler
% Last Mod: 04/02/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_csd = []; data_ecd = []; data_prep = [];
if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters
[parms,data] = check_input(data,varargin);

% replace bad chans, apply smoothing
data_prep = prep_data(data,parms);

% calculate csd
data_csd = calc_csd(data_prep,parms);

% calculate ecd
data_ecd = calc_ecd(data_csd,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,data] = check_input(data,options)
  parms = mmil_args2parms(options,{...
  ... % data preparation
    'badchans',[],[],...
    'chan_wf',0.1,[0,1],...
    'smoothing',0,[],...
  ... % csd calculation method
    'hamming_flag',false,[false true],...
    'edge_corr_flag',true,[false true],...
  ... % csd inverse
    'SNR',20,[1e-10,1e10],...
    'dz',0.15,[],...
    'd1',0,[],...
    'd2',10,[],...
    'sigma',0.3,[],...
  ... % csd hamming
    'r',20,[],... 
...
    'prep_tags',{'badchans','chan_wf','smoothing'},[],...
  });

  % check data matrix
  if ndims(data)~=2
    error('data matrix must be 2D (nchans x ntpoints)');
  end;
  [parms.nchans,parms.ntpoints] = size(data);
  
  % channel indices for electrodes and local field potential
  parms.pg_chans = [1:parms.nchans];
  parms.nchans_lfp = parms.nchans + 1;
  parms.lfp_chans = [1:parms.nchans_lfp];

  % distance along the electrode in mm
  parms.z = (parms.lfp_chans-1)*parms.dz; 
  ind_mid = round(parms.nchans_lfp/2);
  parms.z_cent = parms.z - parms.z(ind_mid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = prep_data(data,parms)
  if isempty(parms.badchans) && parms.smoothing==0, return; end;
  args = mmil_parms2args(parms,parms.prep_tags);
  data = ts_prep_laminar(data,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data_csd = calc_csd(data,parms)
  if parms.hamming_flag
    data_csd = calc_csd_ham(data,parms);
  else
    data_csd = calc_csd_inv(data,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data_csd,data_lfp] = calc_csd_ham(data,parms)
  data_csd = zeros(parms.nchans_lfp,parms.ntpoints);
  data = -data;
  if parms.edge_corr_flag
    data_pad = zeros(parms.nchans+4,parms.ntpoints);
    data_pad(3:parms.nchans+2,:) = data;
    for k=1:parms.nchans+1
      data_csd(k,:)=.23*(data_pad(k,:)-data_pad(k+1,:))+.54*(data_pad(k+1,:)-data_pad((k+2),:))+.23*(data_pad((k+2),:)-data_pad((k+3),:));
    end
  else
    for k=2:(parms.nchans-2)
      data_csd(k+1,:)=.23*(data((k-1),:)-data(k,:))+.54*(data(k,:)-data((k+1),:))+.23*(data((k+1),:)-data((k+2),:));
    end
  end;
  data_csd=data_csd./(parms.r*parms.dz^2);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data_csd,data_lfp] = calc_csd_inv(data,parms)
  % calculate forward and inverse operators
  [VtoL,CtoL,LtoC] = calc_operators(parms);

  % initialize output
  data_csd = zeros(parms.nchans_lfp,parms.ntpoints);
  data_lfp = zeros(parms.nchans_lfp,parms.ntpoints);

  % Mfirst and Mlast used to account for the surface terms when estimating CSD
  M0 = 0.5*(1-(parms.z-parms.z(1))./sqrt((parms.z-parms.z(1)).^2+parms.d2^2/4));
  M1 = 0.5*(1-abs(parms.z-parms.z(end))./sqrt((parms.z-parms.z(end)).^2+parms.d2^2/4)); 

  % calculate local field potential
  data_lfp=VtoL*data;

  % remove surface contribution to the field
  if parms.edge_corr_flag
    data_lfp = correct_potential(data_lfp,M0,M1,CtoL,parms.dz);
  end;

  % calculate CSD from corrected LFP
  data_csd = LtoC*data_lfp;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [VtoL,CtoL,LtoC] = calc_operators(parms)
  VtoL = []; CtoL = []; LtoC = [];

  % construct matrix to calculate LFP from the potential differences (V)
  D = spdiags(ones(parms.nchans,1),0,parms.nchans,parms.nchans);
  D = D+spdiags(-ones(parms.nchans,1),-1,parms.nchans,parms.nchans);
  VtoL = zeros(parms.nchans_lfp,parms.nchans);
  VtoL(2:parms.nchans_lfp,:) = -full(inv(D));  % use minus sign because V_{i}=LFP_{i}-LFP{+1}

  % construct CSD operator (CSD to LFP)
  CtoL = zeros(parms.nchans_lfp,parms.nchans_lfp);
  for i_z = 1:parms.nchans_lfp
    zi=parms.z(i_z);
    % Green's function for potential
    CtoL(i_z,:) = (sqrt((zi-parms.z).^2+parms.d2^2/4) - sqrt((zi-parms.z).^2+parms.d1^2/4));
  end

  % subtract the mean 
  CtoL = CtoL - min(abs(CtoL(:)))*ones(size(CtoL));

  % dimensional forward operator
  A = 1/2/parms.sigma*CtoL*(parms.dz);

  % calcualte the inverse operator (W)
  LtoC = inverse_tikh(A,parms.SNR);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data_ecd = calc_ecd(data_csd,parms)
  % initialize output
  data_ecd = zeros(1,parms.ntpoints);
  % scale CSD by distance from center
  z_data_csd = bsxfun(@times,parms.z_cent',data_csd);
  % evaluate ECD using trapezoidal integration
  data_ecd = trapz(parms.z_cent,z_data_csd,1);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W,R] = inverse_tikh(A,SNR)
  AAt = A*A';
  W = A'/(AAt+1/SNR^2*mean(abs(diag(AAt)))*eye(size(A,1)));
  % R = model resolution matrix = W*A;
  % D = data resolution matrix = A*W;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LFP_view_cor = correct_potential(LFP,M0,M1,CtoL,dz)
  mLFP = repmat(mean(LFP,1),size(LFP,1),1);
  LFP_view = LFP - mLFP;

  dLFPtop = (LFP_view(1,:) - LFP_view(2,:))/dz;
  dLFPbot = (LFP_view(end,:) - LFP_view(end-1,:))/dz;

  LFP_cor1 = 0.5*(CtoL(:,1)*dLFPtop + CtoL(:,end)*dLFPbot);
  LFP_cor2 = M0'*LFP_view(1,:) + M1'*LFP_view(end,:);

  LFP_view_cor = LFP_view-LFP_cor1-LFP_cor2;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
