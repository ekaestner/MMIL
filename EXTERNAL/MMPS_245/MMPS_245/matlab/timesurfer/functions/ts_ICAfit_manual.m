function epoch_data = ts_ICAfit_manual(rawfile,icafile,varargin)
%function ts_ICAfit (rawfile,icafile,outfile,varargin)
%
% Usage:
%  ts_ICAfit (rawfile,icafile,outfile,varargin)
%
% Required Input:
%  rawfile  - full path to epoch_data structure prior to ICA
%                  -or-
%             timesurfer epoch_data structure prior to ICA
%  icafile  - full path to epoch_data structure after ICA
%                  -or-
%             timesurfer epoch_data structure after ICA
%
% Optional Input:
%   outfile  - full path of desired output ICA epoch data structure after 
%              rescaling
%
% Output:
%
%   epoch_data - ICA epoch data scaled to original means and variances
%
% created:         03/18/08    by Andrei Irimia
% last modified:   04/11/08    by Matt Leonard
% last modified:   09/07/09    by Matt Leonard

%% Check Inputs

if nargin < 1, help(mfilename); end

parms = mmil_args2parms(varargin,...
                        { ...
                         'event_codes',[],[],...
                         'conditions',[],[],...
                         'outfile',[],[],...
                         'verbose',true,sort([false true]),...
                         'logfile',[],[],...
                         'logfid',1,[],...
                         'rootoutdir',pwd,[],...
                         'prefix','proc',[],...
                        },...
                        false);
                    

% errors = ts_checkdata(epoch_data);
% if ~isempty(errors)
%     mmil_error(parms,'Errors in provided data structure: %s.',errors);
% end                   
%% rescale the projected data


if ischar(rawfile)
    load(rawfile); raw_epoch_data = epoch_data; clear epoch_data;
elseif isstruct(rawfile)
    raw_epoch_data = rawfile;
else
    error('rawfile neither string nor struct')
end

if ischar(icafile)
    load(icafile); ica_epoch_data = epoch_data; clear epoch_data;
elseif isstruct(icafile)
    ica_epoch_data = icafile;
else
    error('icafile neither string nor struct')
end





% fin_epoch_data = ica_epoch_data;
fprintf('%s: rescaling ICA data to original...\n',mfilename);
for k = 1:length(raw_epoch_data.epochs)
    sz  = size(raw_epoch_data.epochs(1,k).data);
    if numel(sz) == 2
      sz(3) = 1;
    end
%     fin_epoch_data.epochs(1,k).data = [];
    for c = 1:sz(1)
        raw = reshape(raw_epoch_data.epochs(1,k).data(c,:,:),1,sz(2)*sz(3));
        ica = reshape(ica_epoch_data.epochs(1,k).data(c,:,:),1,sz(2)*sz(3));
        mr  = mean(raw); mi = mean(ica);
        warning off matlab:dividebyzero;
        beta  = sum((ica - mi).*(raw - mr))/sum((ica - mi).^2);
        warning on  matlab:dividebyzero
        if isnan(beta), beta = 0; end
        alpha = mr - mi*beta;
        ica_epoch_data.epochs(1,k).data(c,:,:) = reshape(alpha + ica.*beta,1,sz(2),sz(3));
    end
    fprintf('%s: rescaled data in epoch %d.\n',mfilename,k);
end
epoch_data = ica_epoch_data;
if ~isempty(parms.outfile)
save (parms.outfile, '-v7.3', 'epoch_data');
end

return