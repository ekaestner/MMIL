function Y = ts_nan_std(X,flag,dim)
%function Y = ts_nan_std(X,flag,dim)
%
% Purpose: calculate standard deviation of values, not considering NaN values
%
% Usage: same as std
%
% Created:  12/09/13 by Don Hagler
% Last Mod: 05/04/15 by Don Hagler
%
% Based on nan_std by Arnaud Delorme from eeglab, created 09/04/03
%

% Author: Arnaud Delorme, CNL / Salk Institute, Sept 2003

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: nan_std.m,v $
% Revision 1.2  2004/11/23 02:20:52  hilit
% no more divide by zero problem
%
% Revision 1.1  2003/09/04 00:57:11  arno
% Initial revision
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('flag','var') || isempty(flag), flag = 0; end;
if ~exist('dim','var'), dim = []; end;

if isempty(dim)
  if size(X,1) ~= 1
    dim = 1;
  elseif size(X,2) ~= 1
    dim = 2;
  else 
    dim = 3; 
  end;
end;

% set NaN's to 0
ind_nans = find(isnan(X));
X(ind_nans) = 0;

% mask of non-NaN values
N = ones(size(X));
N(ind_nans) = 0;

% number of entries with non-NaN values
n = sum(N,dim);

% row/col with all NaNs
ind_all_nans = find(n==0);
N(ind_all_nans) = NaN;

if flag
  denom = n;
else
  denom = n-1;
end;

Y = sqrt((sum(X.^2,dim)-sum(X,dim).^2./n)./denom);
Y(ind_all_nans) = NaN;

% make sure result is not complex
Y = real(Y);

