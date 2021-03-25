% function [parms crows] = load_csv_parms(fname,varargin)
function [parms crows] = ts_csv2funparms(fname,varargin)
% Purpose: read function parameters from csv file
% Output: parms.(functions).(params) = values
% Modified by JSS on 08-Apr-2009: replaced by ts_csv2parms()
[parms crows,rawdata] = ts_csv2parms(fname,varargin{:});
% 
% CONTROL_OPTIONS = [0 1 2 3];
% % located in second column to right of function name in first column
% % 0 - skip function
% % 1 - run on local machine
% % 2 - run on cluster
% % 3 - auto script only
% 
% % fname     = 'csv_process_fif_data.csv';
% delimiter = ',';
% comment   = '';
% quotes    = '';
% options   = 'empty2NaN';
% [data, result] = mmil_readtext(fname, delimiter, comment, quotes, options);
% 
% % only keep data satisfying the following criteria
% % - col1 & col2 are strings w/o whitespace
% % - col3 is not empty or nan
% % exception: col2 is a CONTROL_OPTIONS
% keep = [];
% cidx = []; % function control
% for row = 1:size(data,1)
%   dat = data(row,1:3);
%   if isstr(dat{1}) && ~any(regexp(dat{1},'\W')) &&...
%      isstr(dat{2}) && ~any(regexp(dat{2},'\W')) &&...
%      ~isempty(dat{3})
%     if isnan(dat{3}), continue; end
%     keep = [keep row];
%   elseif isstr(dat{1}) && ~any(regexp(dat{1},'\W')) &&...
%       ismember(dat{2},CONTROL_OPTIONS)
%     cidx = [cidx row];
%   end
% end
% crows = data(cidx,:);
% data  = data(keep,1:3);
% 
% % for backwards compatibility: check if col3 is also used to specify flags
% % note: col3 would be a cluster flag which overrides the col2 flag
% for c = 1:size(crows,1)
%   dat = crows{c,3};
%   if isnumeric(dat) && ismember(dat,CONTROL_OPTIONS)
%     if dat ~= 0, crows{c,2} = 0; end
%   else
%     crows{c,3} = 0;
%   end
% end
% try
%   crows = crows(ismember(crows(:,1),data(:,1)),:);
% catch
%   crows = [];
% end
% 
% % remove from data: parms of functions w/ control flags = 0
% for c = 1:size(crows,1)
%   if ~any([crows{c,2:3}])
%     lose = strmatch(crows{c,1},data(:,1));
%     keep = setdiff(1:size(data,1),lose);
%     data = data(keep,:);
%   end
% end
% 
% % convert appropriate strings to cells, arrays...
% data = parsecell(data,'all');  
% 
% % make parms
% for n = 1:size(data,1)
%   parms.(data{n,1}).(data{n,2}) = data{n,3};
% end
% 
% %%% END OF FUNCTION
% 
% 



