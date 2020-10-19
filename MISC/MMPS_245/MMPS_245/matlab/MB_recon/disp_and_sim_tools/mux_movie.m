function mux_movie(mux_im, p, fname, map, varargin)
% Syntax:
%   mux_movie(mux_im, p[, fname, map])
%   mux_movie(mux_im, p, fname, map, param, value, param, value...)
%
% Generate a '.avi' movie file for the multiplexed time series.
%
% Inputs
%   mux_im   - Reconstructed SOS combined multiplexed images.
%              Dim: [FE, PE, Echo, Slice, Coil, TemporalPhase].
%   p        - The structure with all parameters of the scan and for the
%              reconstruction. Please refer to mux_epi_params.m for the
%              definition of each field in this structure. The following
%              fields in 'p' are used in this function: mux, FE_DIM, PE_DIM,
%              EC_DIM, SL_DIM, C_DIM, T_DIM.
%   fname    - File name for the .avi movie file. An extension of .avi will
%              be added if the file name does not contain any extensions.
%              Default: mux_movie.avi.
%   map      - The colormap for the movie. Default: gray(256).
%   varargin - Input parameters for the function movie2avi. 
%              For the available parameters, see function movie2avi.
%
% Example
%   mux_movie(dat, p, 'mux_movie.avi', gray(256), 'fps', 3);
%
% (c) Kangrong Zhu,     Stanford University	July 2012

% -- Parse inputs
if ~exist('fname','var') || isempty(fname)
    fname = 'mux_movie.avi';
end
if ~exist('map','var') || isempty(map)
    map = gray(256);
end

% -- Constants
INDEXED_IMAGE = 1;       % The input images to function 'immovie' is an 
                         % m-by-n-by-1-by-k array if the images are indexed 
                         % images(i.e. not true color images).
MAX_INT_FOR_UINT8 = 255; % Maximum integer value for data type 'uint8'

% -- Data type conversion 
% (immovie converts non uint8 data into uint8 data, but the non uint8 data
% must be within the range [0, 255]. The colormap usually only has 256 or
% fewer levels, so this should have little impact on the final movie quality.)
mov_im = mux_im./max(abs(mux_im(:)));
mov_im = uint8(mov_im .* MAX_INT_FOR_UINT8);

% -- Make the movie
datsz = get_dat_sz(mov_im, p);

mov_im = reshape(mov_im, [datsz.x, datsz.y, datsz.ec, datsz.sl/p.mux, p.mux, datsz.c, datsz.t]);
mov_im = permute(mov_im, [1,3,4,6,2,5,7]);

mov_im = reshape(mov_im, [datsz.x * datsz.ec * (datsz.sl/p.mux) * datsz.c, ...
    datsz.y * p.mux, INDEXED_IMAGE, datsz.t]);

mov = immovie(mov_im, map);

if nargin < 5
    movie2avi(mov, fname);
else
    movie2avi(mov, fname, varargin{:});
end
