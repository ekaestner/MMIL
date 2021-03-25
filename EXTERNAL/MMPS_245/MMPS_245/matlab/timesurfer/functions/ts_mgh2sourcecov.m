function R = ts_mgh2sourcecov(fname,varargin)
% function R = ts_mgh2sourcecov(fname,[options])
%
% Usage:
%  R = ts_mgh2sourcecov(fname,'key1',value1,...);
%
% Required input:
%  fname - full or relative path name of mgh file (FreeSurfer surface data)
%
% Optional parameters:
%  'decdips': vector of 0's and 1's indicating which vertices should be included
%    in source covariance matrix; if empty, assume all
%    {default = []}
%  'thresh_flag': [0|1] binarize input values with a threshold
%       and set values for subthreshold and suprathreshold elements
%     otherwise, use input values as they are (scaled so max is maxvar)
%    {default = 1}
%  'thresh': threshold value
%    {default = 0}
%  'absflag': [0|1] apply threshold to absolute values
%    {default = 1}
%  'maxvar': value assigned to diagonal element of covariance matrix for
%    suprathreshold vertices
%    {default = 0.9}
%  'minvar': value assigned to diagonal for subthreshold vertices
%    {default = 0.09}
%  'ncomponents': number of components per dipole (e.g. x,y,z)
%    {default = 3}
%
%  Created:  03/01/11 by Don Hagler
%  Last Mod: 10/23/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse options
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'decdips',[],[],...
  'thresh_flag',true,[false true],...
  'thresh',0,[],...
  'absflag',true,[false true],...
  'maxvar',0.9,[],...
  'minvar',0.09,[],...
  'ncomponents',3,[],...
});
R=[];
if ~exist(fname,'file'), error('file %s not found',mfilename); end

% load file
w = mmil_rowvec(fs_load_mgh(fname));

% open dec file
if isempty(parms.decdips)
  parms.decdips = ones(length(w));
end
n_dips = length(parms.decdips);
v_dec = find(parms.decdips);
n_decdips = length(v_dec);

if parms.thresh_flag
  % threshold values
  w = ts_thresh_weights(w,[],parms.thresh,parms.absflag);
  w(w>0) = parms.maxvar;
  w(w==0) = parms.minvar;
else
  w = parms.maxvar * abs(w) / max(abs(w));
  w(w<parms.minvar) = parms.minvar;
end;

% reduce to dec dips
w_dec = w(v_dec);

% create sparse matrix with ncomponents per dip
W = repmat(w_dec,[parms.ncomponents,1]);
r = reshape(W,[1,prod(size(W))]);
R = diag(sparse(r));

