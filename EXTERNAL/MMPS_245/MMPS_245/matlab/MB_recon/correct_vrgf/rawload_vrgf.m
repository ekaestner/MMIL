function [ramp_flt, uid] = rawload_vrgf(nx_ramp, nx_pres, fname)
%
% function [ramp_flt, uid] = rawload_vrgf(nx_ramp, [nx_pres, fname])
%
% Load the GE vrgf.dat file for the ramp sampling filter in EPI.
%
% Inputs:
%   nx_ramp  - Size in FE, with ramp sampling on.
%   nx_pres  - Prescribed size in FE.
%   fname    - File name. Default: 'vrgf.dat'
%
% Output:
%   ramp_flt - The filter for correcting the ramp sampling, size: [nx_pres, nx_ramp].
%              ksp_corrected(nx_pres-by-ny) = ramp_flt * ksp(nx_ramp-by-ny).
%   uid      - the uinique id from the corresponding p-file (if available)
%              (will be empty if the file does not have a uid header)
%
% (c) Kangrong Zhu,     Stanford University     June 2012

if ~exist('nx_pres', 'var')
    nx_pres = [];
end
if ~exist('fname', 'var')
    fname = 'vrgf.dat';
end

% Open file
fip = fopen(fname, 'r', 'l');
if fip == -1
  error('File %s not found\n', fname);
end

% Attempt to read the 32-byte UID header
uid = fread(fip, 32, 'unsigned char')';
% Check for the magic 7 bytes that the 32-byte UIDs always begin with
if ~all(uid(1:7) == [43 59 149 27 34 71 42])
    uid = [];
    frewind(fip);
end

% Read data
if isempty(nx_pres)
    [ramp_flt, nread] = fread(fip, Inf, 'float32');
    nx_pres = floor(nread/nx_ramp);
    nexpected = nx_ramp * nx_pres;
else
    nexpected = nx_ramp * nx_pres;
    [ramp_flt, nread] = fread(fip, nexpected, 'float32');
end
fclose(fip);
if nread == nexpected
    ramp_flt = permute( reshape(ramp_flt, nx_ramp, nx_pres), [2,1] );
else
    error('%d of %d expected data points are read.', nread, nexpected);
end

return
