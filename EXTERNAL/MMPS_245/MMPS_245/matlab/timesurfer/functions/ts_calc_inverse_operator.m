function [M,nnf]=ts_calc_inverse_operator(G,varargin)
%function [M,nnf]=ts_calc_inverse_operator(G,[options])
%
% ts_calc_inverse_operator: calculate dSPM inverse operator
%   using cortically constrained minimum norm solution
%
% Required Input:
%   G: forward gain matrix of sensor amplitudes for each dipole
%      (or dipole orientation)   (num_sensors x num_sources)
%
% Optional Input:
%   'prewhiten_flag': [0|1] prewhiten gain matrix based on sensor covariance
%     before inverse calculation and use SVD method
%     {default = 0}
%   'SNR': estimated signal to noise ratio
%     {default = 10}
%   'C': sensor covariance matrix used for inverse calculation
%     {default = identity matrix}
%   'C_nn': sensor covariance matrix used for noise normalization
%     ignored if prewhiten_flag = 1
%     if empty, will use C
%     {default = []}
%   'R': source covariance matrix
%     {default = scaled identity matrix}
%   
% Output:
%   M: inverse operator
%   nnf: vector of source scaling factors (noise normalization factor)
%
% Acknowledgements:
%   method from: Dale et al Neuron, Vol. 26, 55-67, April, 2000
%     with prewhiten_flag = 1:
%       notation based on MNE manual by M. Hamalainen Chap 4
%       code modified from mn_inverse_operator.m by M. Huang, Ph.D. 05/19/2005
%     with prewhiten_flag = 0:
%      modified by: Anders Dale, 01/20/2008 and 02/24/2011
%
% Created:  08/02/05 by Don Hagler
% Last Mod: 10/29/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=[]; nnf=[];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(G,varargin);

% calculate inverse operator and noise normalization factors
if ~parms.prewhiten_flag
  [M,nnf] = calc_inverse_amd(G,parms);
else
  [M,nnf] = calc_inverse_prewhiten(G,parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(G,options)
  if isempty(G)
    error('gain matrix G is empty');
  end
  parms = mmil_args2parms(options,{...
    'prewhiten_flag',false,[false true],...
    'SNR',10,[eps Inf],...
    'C',[],[],...
    'C_nn',[],[],...
    'R',[],[],...    
  });
  [num_sensors,num_sources] = size(G);
  if isempty(parms.C), parms.C = eye(num_sensors); end;
  if isempty(parms.C_nn), parms.C_nn = parms.C; end;
  if isempty(parms.R)
    parms.R = speye(num_sources,num_sources);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M,nnf] = calc_inverse_amd(G,parms)
  % scale source covariance matrix so that trace(GRG')/trace(I)=1
  R = scale_sourcecov(parms.R,G);
  
  % calculate inverse operator
  lamda_sq_SNR = 1.0/parms.SNR^2; % power SNR is needed for regularization
  M = R*G'*inv(G*R*G' + ...
               lamda_sq_SNR*mean(diag(G*R*G'))/mean(diag(parms.C))*parms.C);

  % noise sensitivity normalization
  num_sources = size(G,2);
  nnf = zeros(num_sources,1);
  for k=1:num_sources
    nnf(k) = sqrt(M(k,:)*parms.C_nn*M(k,:)');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M,nnf] = calc_inverse_prewhiten(G,parms)
  % calculate C^-1/2 for prewhitening
  [Uc,gammac] = eig(parms.C);
  rank_C = rank(parms.C);
  sqrt_inv_gammac = zeros(size(gammac));
  num_sensors = size(G,1);
  for i=1:num_sensors
    if (num_sensors-i+1)<=rank_C
      sqrt_inv_gammac(i,i) = sqrt(1.0/gammac(i,i));
    end;
  end;
  Cinv_sqrt = sqrt_inv_gammac*Uc';

  % prewhiten gain matrix
  G = Cinv_sqrt*G;

  % scale source covariance matrix so that trace(GRG')/trace(I)=1
  R = scale_sourcecov(parms.R,G);

  % Cholesky factorization for computational efficiency
  Rc = chol(R);             

  % calculate inverse operator with singular value decomposition
  A = G*Rc;
  [V,lamda,U] = svd(A',0);  % economy SVD preventing memory overflow
  lamda = diag(lamda);
  lamda_sq_SNR = 1.0/parms.SNR^2; % power SNR is needed for regularization
  gamma = diag(lamda./(lamda.^2+lamda_sq_SNR));
  gamma = sparse(gamma);
  Vbar_gamma = Rc*V*gamma;
  M = Vbar_gamma*U'*Cinv_sqrt; % this is the inverse operator

  % noise sensitivity normalization
  num_sources = size(G,2);
  nnf = zeros(num_sources,1);
  for k=1:num_sources
    nnf(k) = sqrt(sum(Vbar_gamma(k,:).^2));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = scale_sourcecov(R,G)
  sf = trace(G*R*G');
  if isnan(sf)
    error('trace of GRG'' is NaN');
  end;
  num_sensors = size(G,1);
  R = num_sensors*R/sf;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

