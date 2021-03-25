function [E,norm_var_E,...
Yeeg,Yfiteeg,Eeeg,var_Eeeg,var_Yeeg,norm_var_Eeeg,...
Ygrad,Yfitgrad,Egrad,var_Egrad,var_Ygrad,norm_var_Egrad,...
Ymag,Yfitmag,Emag,var_Emag,var_Ymag,norm_var_Emag]...
= rc_RCSE_calc_fit_err(Y,Yfit,parms)
%
% Purpose: dipole optimization using specified vectors of
%   r and th offsets (exhaustive search)
%
% Required Input:
%   Y: data matrix
%   Yfit: fitted data matrix
%   parms: RCSE parms struct
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 05/16/11 by Don Hagler
%

% return struct?

if ~mmil_check_nargs(nargin,3), return; end;

fprintf('%s: calculating fit and residual error...\n',mfilename);
% calculate separate error for each channel type
Y2 = reshape(Y,[size(Y,1),length(parms.goodchans),parms.nconds]);
E = Y-Yfit;
Yfit2 = reshape(Yfit,[size(Yfit,1),length(parms.goodchans),parms.nconds]);
norm_var_E = ones(size(Y,1),1);
ntypes = 0;
if parms.useEEG_flag
  [good_chans,ind,ind_good_eeg] = intersect(parms.EEG_chans,parms.goodchans);
  Yeeg = Y2(:,ind_good_eeg,:);
  Yfiteeg = Yfit2(:,ind_good_eeg,:);
  Eeeg = Yeeg - Yfiteeg;
  Yeeg = reshape(Yeeg,[size(Yeeg,1),size(Yeeg,2)*size(Yeeg,3)]);
  Yfiteeg = reshape(Yfiteeg,[size(Yfiteeg,1),size(Yfiteeg,2)*size(Yfiteeg,3)]);
  Eeeg = reshape(Eeeg,[size(Eeeg,1),size(Eeeg,2)*size(Eeeg,3)]);
  var_Eeeg = var(Eeeg,0,2);
  var_Yeeg = var(Yeeg,0,2);
  var_Yfiteeg = var(Yfiteeg,0,2);
  tmp_var_Yeeg = var_Yeeg;
  tmp_var_Yeeg(find(~var_Yeeg))=1;
  norm_var_Eeeg = var_Eeeg/max(var_Yeeg);
  norm_var_E = norm_var_E.*norm_var_Eeeg;
  ntypes = ntypes + 1;
else
  Yeeg = [];
  Yfiteeg = [];
  Eeeg = [];
  var_Eeeg = [];
  var_Yeeg = [];
  norm_var_Eeeg = [];
end;
if parms.usegrad_flag
  [good_chans,ind,ind_good_grad] = intersect(parms.grad_chans,parms.goodchans);
  Ygrad = Y2(:,ind_good_grad,:);
  Yfitgrad = Yfit2(:,ind_good_grad,:);
  Egrad = Ygrad - Yfitgrad;
  Ygrad = reshape(Ygrad,[size(Ygrad,1),size(Ygrad,2)*size(Ygrad,3)]);
  Yfitgrad = reshape(Yfitgrad,[size(Yfitgrad,1),size(Yfitgrad,2)*size(Yfitgrad,3)]);
  Egrad = reshape(Egrad,[size(Egrad,1),size(Egrad,2)*size(Egrad,3)]);
  var_Egrad = var(Egrad,0,2);
  var_Ygrad = var(Ygrad,0,2);
  tmp_var_Ygrad = var_Ygrad;
  tmp_var_Ygrad(find(~var_Ygrad))=1;
  norm_var_Egrad = var_Egrad/max(var_Ygrad);
  norm_var_E = norm_var_E.*norm_var_Egrad;
  ntypes = ntypes + 1;
else
  Ygrad = [];
  Yfitgrad = [];
  Egrad = [];
  var_Egrad = [];
  var_Ygrad = [];
  norm_var_Egrad = [];
end;
if parms.usemag_flag
  [good_chans,ind,ind_good_mag] = intersect(parms.mag_chans,parms.goodchans);
  Ymag = Y2(:,ind_good_mag,:);
  Yfitmag = Yfit2(:,ind_good_mag,:);
  Emag = Ymag - Yfitmag;
  Ymag = reshape(Ymag,[size(Ymag,1),size(Ymag,2)*size(Ymag,3)]);
  Yfitmag = reshape(Yfitmag,[size(Yfitmag,1),size(Yfitmag,2)*size(Yfitmag,3)]);
  Emag = reshape(Emag,[size(Emag,1),size(Emag,2)*size(Emag,3)]);
  var_Emag = var(Emag,0,2);
  var_Ymag = var(Ymag,0,2);
  tmp_var_Ymag = var_Ymag;
  tmp_var_Ymag(find(~var_Ymag))=1;
  norm_var_Emag = var_Emag/max(var_Ymag);
  norm_var_E = norm_var_E.*norm_var_Emag;
  ntypes = ntypes + 1;
else
  Ymag = [];
  Yfitmag = [];
  Emag = [];
  var_Emag = [];
  var_Ymag = [];
  norm_var_Emag = [];
end;
norm_var_E = norm_var_E.^(1/ntypes);

