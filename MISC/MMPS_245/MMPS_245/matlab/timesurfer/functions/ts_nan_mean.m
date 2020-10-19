function Y = ts_nan_mean(X,dim)
%function Y = ts_nan_mean(X,dim)
%
% Purpose: average values, not considering NaN values
%
% Usage: same as mean
%
% Created:  02/09/12 by Don Hagler
% Last Mod: 12/09/13 by Don Hagler
%
% Based on nan_std by Arnaud Delorme from eeglab, created 10/16/02
%

% Author: Arnaud Delorme, CNL / Salk Institute, 16 Oct 2002

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: nan_mean.m,v $
% Revision 1.3  2004/11/23 02:23:02  hilit
% no more divide by zero problem
%
% Revision 1.2  2002/10/17 18:43:16  arno
% debugging dim
%
% Revision 1.1  2002/10/17 02:34:52  arno
% Initial revision
%

if ~mmil_check_nargs(nargin,1), return; end;

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
tmpX = X;
tmpX(find(isnan(X(:)))) = 0;
denom = sum(~isnan(X),dim);
denom(find(~denom)) = nan;
Y = sum(tmpX, dim) ./ denom;

