function dat = notch_filter(dat, nskip, notch_thresh)
% function dat_filtered = notch_filter(dat, nskip, thresh)
%
% Simple notch filter for a timeseries of k-space data to help reduce the
% effects of white-pixel noise.
%
% Inputs
%   dat    - Input data. Dim: [FE, PE, Echo, Slice, Coil, Time].
%   nskip  - Number of time frames in the beginning to ignore. (Typically
%   the number of calibration scans, plus 2 if external cal was used.)
%   thresh - Notch threshold (z-score units). voxels with intensities below
%   this value will be replaced by the corresponding value from a
%   neighboring time point.
%
% Output
%   dat_filt - Filtered data (same dims as dat).
%
%
% Notes:
% We don't try to filter the cal data. That would require a different
% (more complicated) algorithm. And it *might* not be necessary?
%
% (c) 2014 Bob Dougherty, Stanford University

% Operate on coil-combined (using SS) magnitude data:

dat_start = nskip + 1;
if size(dat,6) > dat_start+1
    datc = abs(sqrt(sum(dat.^2, 5)));
    datc(:,:,:,:,:,dat_start:end) = NaN;

    mn = nanmean(datc,6);
    sd = nanstd(datc,0,6);
    z = NaN(size(datc));
    for t=dat_start:size(datc,6)
        z(:,:,:,:,:,t) = (datc(:,:,:,:,:,t)-mn)./sd;
    end

    % Do it again after the really bad outliers are removed:
    bad_vals = abs(z)>abs(notch_thresh);
    datc(bad_vals) = NaN;
    mn = nanmean(datc,6);
    sd = nanstd(datc,0,6);
    z = NaN(size(datc));
    for t=dat_start:size(datc,6)
        z(:,:,:,:,:,t) = (datc(:,:,:,:,:,t)-mn)./sd;
    end

    %% try to avoid using frames with lots of motion, or the corrupt
    %% frames that we sometimes get at the end of an aborted scan
    %bad_frames = squeeze(sum(sum(z>3)))./(size(z,1)*size(z,2)) > 0.01;
    %if any(bad_frames)
    %    bad_frames(1:dat_start) = 1;
    %    mn = mean(datc(:,:,:,:,:,~bad_frames),6);
    %    sd = std(datc(:,:,:,:,:,~bad_frames),0,6);
    %    z = zeros(size(datc));
    %    for t=1:size(datc,6)
    %        if ~bad_frames(t)
    %            z(:,:,:,:,:,t) = (datc(:,:,:,:,:,t)-mn)./sd;
    %        end
    %    end
    %end

    [x,y,e,s,c,t] = ind2sub(size(datc), find(z<notch_thresh));
    fprintf('mux_epi_main: Fixing %d notches\n', numel(x));
    %if p.debug
    %    disp(unique(t)')
    %end
    % Apply the filter.
    for ii=1:numel(x)
        % Replace this pixel in each of the coil images.
        for jj=1:size(dat,5)
            % For all but the first time point, replace with the
            % corresponding value from the previous time point. (Note that
            % if there are notches in the same voxel in two adjact
            % timepoints t and t+1, this will have the effect of replacing
            % both t and t+1 with the value in t-1.)
            if(t(ii)==dat_start)
                dat(x(ii),y(ii),e(ii),s(ii),jj,t(ii)) = dat(x(ii),y(ii),e(ii),s(ii),jj,t(ii)+1);
            else
                dat(x(ii),y(ii),e(ii),s(ii),jj,t(ii)) = dat(x(ii),y(ii),e(ii),s(ii),jj,t(ii)-1);
            end
        end
    end
end

return

function out = nanmean(d, dim)
if nargin < 2
    if size(d,1) ~= 1
        dim = 1;
    elseif size(d,2) ~= 1
        dim = 2;
    else
        dim = 3;
    end
end
nans = isnan(d);
d(nans) = 0;
out = sum(d, dim) ./ sum(~nans,dim);
return

function out = nanstd(d, flag, dim)
if nargin==1
    flag = 0;
elseif isempty(flag)
    flag = 0;
end
if nargin< 3,
    if size(d,1) ~= 1
        dim = 1;
    elseif size(d,2) ~= 1
        dim = 2;
    else
        dim = 3;
    end
end
tmp = isnan(d);
d(tmp) = 0;
tmp = ~tmp;
tmp = sum(tmp, dim);
tmp(tmp==0) = 1;

out = sqrt((sum(d.^2, dim) - sum(d, dim).^2./tmp)./(tmp-abs(flag-1)));
%out(nononnans) = NaN;
return

