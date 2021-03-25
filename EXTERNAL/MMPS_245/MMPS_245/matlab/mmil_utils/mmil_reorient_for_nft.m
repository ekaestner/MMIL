function vol = mmil_reorient_for_nft(fname,tmpdir,forceflag)
%function vol = mmil_reorient_for_nft(fname,tmpdir,forceflag)
%
% Created:  09/27/13 by Don Hagler
% Last Mod: 02/17/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('tmpdir','var') || isempty(tmpdir), tmpdir = pwd; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

mmil_mkdir(tmpdir);

[tmp,fstem] = fileparts(fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reorient to orientation expected by NFT
%   use fs_reorient so it is easily reversible (unlike mri_convert)
fname_out = sprintf('%s/%s_PSR.mgh',tmpdir,fstem);
[vol,M] = fs_load_mgh(fname);
[vol,M] = fs_reorient(vol,M,'PSR');
fs_save_mgh(vol,fname_out,M);
fname = fname_out;

% convert to analyze
fname_out = sprintf('%s/%s_PSR.img',tmpdir,fstem);
fs_mri_convert(fname,fname_out,...
  'options','-ot spm','forceflag',forceflag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reorient to axial
[vol,pixdim,dtype] = segm_readanalyze(fname_out);
[K,L,M] = size(vol);
tmp = vol; clear vol;
for i = 1:M
  vol(:,i,:) = (reshape(tmp(:,:,i), K, L));
end
[K,L,M] = size(vol);
for i = 1:M
  vol(:,:,i) = rot90(vol(:,:,i),3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

