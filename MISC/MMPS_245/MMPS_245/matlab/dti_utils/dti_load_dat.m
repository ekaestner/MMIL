function vol = dti_load_dat(fname,volsz,permvec,revsliceflag,datatype)
%function vol = dti_load_dat(fname,volsz,[permvec],[revsliceflag],[datatype])
%
% Purpose: load raw data dat file (e.g. from DTIstudio)
%
% Required Input:
%   fname: full or relative path of dat file
%   volsz: vector of dimension sizes for volume in dat file
%
% Optional Input:
%   permvec: vector of integers defining how volume should be
%            permuted (e.g. [2,3,4,1])
%   datatype: data type of values (e.g. 'float', 'uint8', etc.);
%     (default = 'float')
%   revsliceflag: reverse slice sequencing (after applying permvec)
%
% Created:  10/29/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

vol = [];
if ~mmil_check_nargs(nargin,2), return; end;

if ~exist('permvec','var'), permvec=[]; end;
if ~exist('datatype','var'), datatype='float'; end;
if ~exist('revsliceflag','var'), revsliceflag=0; end;

if exist(fname,'file')
  fid=fopen(fname,'r');
  vol=fread(fid,datatype);
  fclose(fid);
else
  error('file %s not found',fname);
end;

if length(vol)==prod(volsz)
  vol=reshape(vol,volsz);
else
  error('%d elements read from file %s but %d indicted by volsz',...
    length(vol),fname,prod(volsz));
end;

if ~isempty(permvec), vol = permute(vol,permvec); end;
if revsliceflag
  vol = vol(:,:,end:-1:1,:);
end;

