function dti_dat2mgh(fname,volsz,M,permvec,revsliceflag,datatype)
%function dti_dat2mgh(fname,volsz,[M],[permvec],[revsliceflag],[datatype])
%
% Purpose: convert raw data dat file to mgh file
%
% Required Input:
%   fname: full or relative path of dat file
%   volsz: vector of dimension sizes for volume in dat file
%
% Optional Input:
%   M: 4x4 transformation matrix defining slice orientation, offsets, etc.
%     {default = eye(4)}
%   permvec: vector of integers defining how volume should be
%            permuted (e.g. [2,3,4,1])
%     {default = []}
%   revsliceflag: toggle [0|1] whether to reverse slice order of input volume
%            if omitted or empty, no reversal
%     {default = 0}
%   datatype: data type of values (e.g. 'float', 'uint8', etc.);
%     (default = 'float')
%
% Early Mod: 03/04/07 by Don Hagler
% Last Mod:  10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

if ~exist('M','var') || isempty(M), M=eye(4); end;
if ~exist('permvec','var'), permvec=[]; end;
if ~exist('revsliceflag','var'), revsliceflag=[]; end;
if isempty(revsliceflag), revsliceflag = 0; end;
if ~exist('datatype','var'), datatype='float'; end;

k=findstr(fname,'.dat');
if ~isempty(k)
  fstem = fname(1:k(end)-1);
else
  error('file %s is not dat file',fname);
end;

if exist(fname,'file')
  fid=fopen(fname,'r');
  vol=fread(fid,datatype);
  vol=reshape(vol,volsz);
  fclose(fid);
else
  error('file %s not found',fname);
end;

if ~isempty(permvec), vol = permute(vol,permvec); end;
if revsliceflag
  vol = vol(:,:,end:-1:1,:);
end;

fname_out = [fstem '.mgh'];
fs_save_mgh(vol,fname_out,M);

