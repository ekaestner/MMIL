function dti_save_dat(vol,fname,permvec,revsliceflag,datatype)
%function dti_save_dat(vol,fname,[permvec],[revsliceflag],[datatype])
%
% Purpose: save volume as raw data file
%
% Input:
%   vol: volume (3D or 4D)
%   fname: full or relative path of output file
%   permvec: vector of integers defining how volume should be
%            permuted (e.g. [4,2,1,3])
%            if omitted or empty, no permutation
%            {default = []}
%   revsliceflag: toggle [0|1] whether to reverse slice order of input volume
%            if omitted or empty, no reversal
%            {default = []}
%   datatype: data type of values (e.g. 'float', 'uint8', etc.);
%     (default = 'float')
%
% Created:  08/09/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('permvec','var'), permvec=[]; end;
if ~exist('revsliceflag','var'), revsliceflag=[]; end;
if isempty(revsliceflag), revsliceflag = 0; end;
if ~exist('datatype','var'), datatype='float'; end;

volsz = size(vol);
ndims = length(volsz);
if revsliceflag
  if ndims==4
    vol = vol(:,:,end:-1:1,:);
  else
    vol = vol(:,:,end:-1:1);
  end;
end;
if ~isempty(permvec), vol = permute(vol,permvec); end;

fid=fopen(fname,'w');
if fid<=0
  error('failed to open file %s for writing',fname);
end;
count=fwrite(fid,vol,datatype);
fclose(fid);

