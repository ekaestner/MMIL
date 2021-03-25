function psi_mtx = noise_file_to_cov(fname, nc)
%
% function psi_mtx = noise_file_to_cov(fname, nc)
%
% Calculates coil noise covariance matrix from noise data file.
%
% Inputs
%   fname   - Filename of noise data file (_PXXXXX.7_noise.dat).
%   nc      - Number of physical coils.
%
% Output
%   psi_mtx - Coil noise covariance matrix. Dim: [nc, nc].
%
% (c) Kangrong Zhu, Hua Wu      Stanford University     Dec 2014

NUM_DATA_PTS = 4096;

fid = fopen(fname, 'r', 'l');
if fid == -1
    error('Can not open file %s.\n', fname);
end
noise = fread(fid, inf, 'double');
fclose(fid);

if length(noise) ~= (2 * NUM_DATA_PTS * nc)
    error('Number of data points in noise.dat is wrong.');
end

noise = reshape(noise, 2, NUM_DATA_PTS, nc);
noise = noise(1, :, :) + 1i * noise(2, :, :);

psi_mtx = calc_coil_noise_cov(noise, 3);

return