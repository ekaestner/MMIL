function [nr,nc,ns] = mmil_serinfo_dims(serinfo)
%function [nr,nc,ns] = mmil_serinfo_dims(serinfo)
%
% Required Input:
%   serinfo: series info struct created by mmil_classify_dicoms
%     must be a single struct, not an array
%
% Output:
%   nr: number of rows
%   nc: number of columns
%   ns: number of slices
%
% Created:  09/11/12 by Don Hagler
% Last Mod: 12/21/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
nr = []; nc = []; ns = [];

% for Siemens mosaic:
%   'nr' has correct number of rows (# of phase-encoding lines)
% otherwise:
%   'Rows' is correct

ns = serinfo.nslices;
nr = mmil_getfield(serinfo,'nr',serinfo.Rows);
nc = mmil_getfield(serinfo,'nc',serinfo.Columns);

if isempty(nr), nr = serinfo.Rows; end;
if isempty(nc), nc = serinfo.Columns; end;

ns = double(ns);
nr = double(nr);
nc = double(nc);
