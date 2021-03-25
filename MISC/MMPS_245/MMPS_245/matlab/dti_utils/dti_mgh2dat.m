function dti_mgh2dat(fname,permvec,revsliceflag,datatype,frames)
%function dti_mgh2dat(fname,[permvec],[revsliceflag],[datatype],[frames])
%
% Purpose: convert mgh file to raw data file
%
% Required Input:
%   fname: full or relative path of mgh file
%
% Optional Input:
%   permvec: vector of integers defining how volume should be
%            permuted (e.g. [4,2,1,3])
%            if omitted or empty, no permutation
%     {default = []}
%   revsliceflag: toggle [0|1] whether to reverse slice order of input volume
%            if omitted or empty, no reversal
%     {default = 0}
%   datatype: data type of values (e.g. 'float', 'uint8', etc.);
%     {default = 'float'}
%   frames: vector of frames to convert
%     {default = []}
%
% Early Mod: 03/04/07 by Don Hagler
% Last Mod:  10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('permvec','var'), permvec=[]; end;
if ~exist('revsliceflag','var') || isempty(revsliceflag), revsliceflag=0; end;
if ~exist('datatype','var') || isempty(datatype), datatype='float'; end;
if ~exist('frames','var'), frames = []; end;

[tmp_path,tmp_stem,tmp_ext] = fileparts(fname);
if ~ismember(tmp_ext,{'.mgh','.mgz'})
  error('file %s is not an mgh or mgz file',fname);
end;
if isempty(tmp_path)
  tmp_path = pwd;
end;

if ~exist(fname,'file')
  error('file %s not found',fname);
end;
[vol,M] = fs_load_mgh(fname,[],frames);

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

fname_out=sprintf('%s/%s.dat',tmp_path,tmp_stem);
fid=fopen(fname_out,'w');
if fid<0
  error('failed to open file %s for writing',fname_out);
end;
count=fwrite(fid,vol,datatype);
fclose(fid);

