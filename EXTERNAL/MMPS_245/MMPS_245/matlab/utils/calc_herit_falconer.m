function [vec_herit,vec_corr_Mz,vec_corr_Dz] = calc_herit_falconer(data_Mz,data_Dz);
%function [vec_herit,vec_corr_Mz,vec_corr_Dz] = calc_herit_falconer(data_Mz,data_Dz);
%
% Purpose: calculate Falconer heritability
%
% Required Input:
%   data_Mz: matrix of data for monozygotic twins
%     size must be = [nvals,npairs,2]
%   data_Dz: matrix of data for dizygotic twins
%     size must be = [nvals,npairs,2]
%
% Created:  10/17/07 by Don Hagler
% Last Mod: 10/24/07 by Don Hagler
%

outvec = [];
if nargin<2
  help(mfilename);
  return;
end;

if length(size(data_Mz))~=3 | length(size(data_Dz))~=3
  error('data_Mz and data_Dz have wrong dimensions');
end;
if size(data_Mz,3)~=2 | size(data_Dz,3)~=2
  error('last dimension of data_Mz and data_Dz must have size=2');
end;

nvals = size(data_Mz,1);
if size(data_Dz,1) ~= nvals
  error('mismatch in nvals between data_Mz and data_Dz');
end;

tmp1 = sum(data_Mz,3);
tmp1(tmp1>0)=1;
tmp2 = sum(tmp1,2);
ind_nonzero = find(tmp2>size(tmp1,2)/2);
num_nonzero = length(ind_nonzero);
fprintf('%s: %d non zero values...\n',mfilename,num_nonzero);

vec_corr_Mz = zeros(nvals,1);
vec_corr_Dz = zeros(nvals,1);
for j=1:num_nonzero
  i=ind_nonzero(j);
  tmp_MzA = squeeze(data_Mz(i,:,1))';
  tmp_MzB = squeeze(data_Mz(i,:,2))';
  tmp_DzA = squeeze(data_Dz(i,:,1))';
  tmp_DzB = squeeze(data_Dz(i,:,2))';
  corr_Mz = corr(tmp_MzA,tmp_MzB);
  corr_Dz = corr(tmp_DzA,tmp_DzB);
  vec_corr_Mz(i) = corr_Mz;
  vec_corr_Dz(i) = corr_Dz;
end;

vec_corr_Mz(isnan(vec_corr_Mz))=0;
vec_corr_Dz(isnan(vec_corr_Dz))=0;
vec_herit = (vec_corr_Mz - vec_corr_Dz)*2;


