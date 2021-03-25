function scaled_data = ts_ICAfit_auto (data_in,ica_data,varargin)
%function scaled_data = ts_ICAfit_auto (data_in,ica_data,varargin)
%
% Usage:
%  scaled_data = ts_ICAfit (data_in,ica_data)
%
% Required Input:
%   NOTE: These inputs are only the data field from an
%         epoch_data structure
% 
%  data_in   - original data prior to ICA rejection
%  ica_data  - data after ICA rejection
%
% Output:
%
%   scaled_data - data after fitting
%
% Created:   03/18/08 by Andrei Irimia
% Early Mod: 04/14/08 by Rajan Patel
% Last Mod:  11/17/11 by Don Hagler
%

% 04/14/08 - Work on data only and not files - more universal

%% Check Inputs

if nargin < 2, help(mfilename); end

parms = mmil_args2parms(varargin,...
                        { ...
                         'verbose',true,sort([false true]),...
                         'logfile',[],[],...
                         'logfid',1,[],...
                        },...
                        false);
                    
if size(data_in) ~= size(ica_data)
    mmil_error(parms,'The size of the data sets do not match.');
end

%% rescale the projected data    

[channels samples trials] = size(data_in);
for c = 1:channels
    orig_data  = reshape(data_in (c,:,:),1,samples*trials);
    post_ica   = reshape(ica_data(c,:,:),1,samples*trials);
    mean_orig  = mean(orig_data);
    mean_ica   = mean(post_ica);
    warning off matlab:dividebyzero;
    beta       = sum((post_ica - mean_ica).*(orig_data - mean_orig))/sum((post_ica - mean_ica).^2);
    warning on  matlab:dividebyzero
    if isnan(beta), beta = 0; end
    alpha = mean_orig - mean_ica*beta;
    scaled_data(c,:,:) = reshape(alpha + post_ica.*beta,1,samples,trials);
end
