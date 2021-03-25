function Y = rc_load_ret_avgt(avg_data,cond_order,channels,scale_matrix)
%function Y = rc_load_ret_avgt(avg_data,cond_order,channels,[scale_matrix])
%
% Purpose: load values from average structure
%
% Required Input:
%   avg_data: average structure created by ts_average_fif_data
%   cond_order: vector of condition numbers indicating how to order the conditions
%   channels: vector of channel/sensor indices -- if empty, get all sensors
%
% Optional Input:
%   scale_matrix: matrix of scale factors for each channel and time point
%               if empty or omitted, no scaling will be done
%
% Output:
%   Y: matrix of values for each timepoint for each sensor and each condition
%
% created:  12/01/05 by Don Hagler
% last mod: 02/19/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('channels','var'), channels = []; end;
if ~exist('scale_matrix','var'), scale_matrix = []; end;

Y=[];

if(isempty(avg_data))
  error('avg_data structure is empty');
end

if(length(avg_data.averages)<max(cond_order) | ...
   length(avg_data.averages)<length(cond_order) | min(cond_order) < 0)
  error('something is wrong with cond_order');
end

[num_sensors,num_tpoints] = size(avg_data.averages(1).data);

% if not specified, use all channels
if isempty(channels)
  channels = 1:num_sensors;
end
if(length(channels)>num_sensors | max(channels)>num_sensors | min(channels)<0)
  error('something is wrong with channels');
end

if size(scale_matrix,1)~=num_sensors | size(scale_matrix,2)~=num_tpoints
  error('scale_matrix is wrong size');
end

k1=1;
ntheta = length(cond_order);
nchannels = length(channels);
Y=zeros(num_tpoints,nchannels*ntheta,1);
for j=1:ntheta
  k2=k1+nchannels-1;
  theta=cond_order(j);
  if(theta==0)
    Y(:,k1:k2) = 0;
  else
    if isempty(scale_matrix)
      Y(:,k1:k2) = avg_data.averages(theta).data(channels,:)';
    else
      data = avg_data.averages(theta).data;
      data = data.*scale_matrix;
      Y(:,k1:k2) = data(channels,:)';
    end;
  end
  k1=k2+1;
end

return;

