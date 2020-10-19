function dat = homodyne(dat, ny_part, ntran, niter, in_dom, out_dom)
%
% function dat = homodyne(dat, ny_part, [ntran=2], [niter=4], [in_dom='KxKy'], [out_dom=in_dom])
%
% Homodyne detection reconstruction of partially acquired k-space.
%
% Inputs
%   dat     - Partially acquired k-space, with full matrix size. The
%             acquired k-space data is FOLLOWED by zeros in the PE direction.
%             Dim: [FE, PE, Echo, Slice, Coil, TemporalPhase].
%   ny_part - # of acquired ky lines in the partial k-space acquisition.
%   ntran   - Number of points in the transition band of the filters.
%   niter   - Number of iterations.
%             0: Use zero-filling.
%             1,2,3...: Use Homodyne. Iterative when niter > 1.
%   in_dom  - A string specifying the domain of the input 'dat'. 'xy', 'xKy', 'KxY' or 'KxKy'(uppercase or lowercase doesn't matter).
%   out_dom - A string specifying the domain of the output 'dat'.
%
% Output
%   dat     - Reconstructed k-space data, having the same size and in the same domain as the input 'dat'.
%
% From GE company. Based on Homodyne.cpp.
% Modified by Robert Dougherty and Kangrong Zhu, Stanford University, July 2014

if ~exist('niter', 'var') || isempty(niter)
    niter = 4;
end
if niter < 0
    error('Number of iterations must be no smaller than 0.');
end
if ~exist('ntran', 'var') || isempty(ntran)
    ntran = 2;
end
if ~exist('in_dom', 'var') || isempty(in_dom)
    in_dom = 'KxKy';
end
if ~exist('out_dom', 'var') || isempty(out_dom)
    out_dom = in_dom;
end
in_dom = lower(in_dom);
out_dom = lower(out_dom);

FE_DIM = 1;
PE_DIM = 2;

% Data size
datsz = size(dat);
ny_pres = datsz(PE_DIM);
if length(datsz) < 7
    datsz(end+1 : 7) = 1;        % To make sure datsz(7) exists when it is called below
end

% Transform input data into (x/Kx, Ky)
switch in_dom
    case 'xy'
        dat = fftc(dat, PE_DIM); % (x, Ky)
    case 'xky'
        % Do nothing             % (x, Ky)
    case 'kxy'
        dat = fftc(dat, PE_DIM); % (Kx, Ky)
    case 'kxky'
        % Do nothing             % (Kx, Ky)
    otherwise
        error('Input in_dom string must be ''xy'' or ''xky'' or ''kxy'' or ''kxky''.');
end

% Make sure unacquired points are set to zero
dat(:, ny_part+1:end, :, :, :, :, :) = 0;

% If no Homodyne(just zero-filling), output the data in the specified domain.
if niter == 0
    if strcmp(out_dom, 'xky') || strcmp(out_dom, 'kxky')
        % Do nothing along phase encoding direction
    end
    if strcmp(out_dom, 'xy') || strcmp(out_dom, 'kxy')
        dat = ifftc(dat, PE_DIM);
    end
    if (strcmp(in_dom, 'xy') && strcmp(out_dom, 'kxy')) || (strcmp(in_dom, 'xky') && strcmp(out_dom, 'kxy')) || ...
            (strcmp(in_dom, 'xy') && strcmp(out_dom, 'kxky')) || (strcmp(in_dom, 'xky') && strcmp(out_dom, 'kxky'))
        dat = fftc(dat, FE_DIM);
    end
    if (strcmp(in_dom, 'kxy') && strcmp(out_dom, 'xy')) || (strcmp(in_dom, 'kxky') && strcmp(out_dom, 'xy')) || ...
            (strcmp(in_dom, 'kxy') && strcmp(out_dom, 'xky')) || (strcmp(in_dom, 'kxky') && strcmp(out_dom, 'xky'))
        dat = ifftc(dat, FE_DIM);
    end
    
    return
end

% Keep a copy of the input data in (x, Ky)
dat_input = dat;                                   % If in_dom is 'xy' or 'xky', dat_input is in (x, Ky); If in_dom is 'kxy' or 'kxky', dat_input is in (Kx, Ky)
if strcmp(in_dom, 'kxy') || strcmp(in_dom, 'kxky')
    dat_input = ifftc(dat_input, FE_DIM);          % dat_input is in (x, Ky)
end

% Generate filters
[lowPassFilter, highPassFilter, mergeWindow] = generate_filters(ntran, ny_part, ny_pres);

lowPassFilter(ny_part+1 : ny_pres) = 0.0;
highPassFilter(ny_part+1 : ny_pres) = 0.0;

lowPassFilter = repmat( lowPassFilter,[size(dat,1) 1 size(dat,3) size(dat,4) size(dat,5)]);
highPassFilter = repmat(highPassFilter,[size(dat,1) 1 size(dat,3) size(dat,4) size(dat,5)]);

% Calculate high and low pass images in (x/Kx, Ky)
lowPass = zeros(datsz);
for i7 = 1 : datsz(7)
    for i6 = 1 : datsz(6)
        % Get low pass views
        lowPass(:, :, :, :, :, i6, i7)  = dat(:, :, :, :, :, i6, i7) .* lowPassFilter;
        % Get high pass views
        dat(:, :, :, :, :, i6, i7) = dat(:, :, :, :, :, i6, i7) .* highPassFilter;
    end
end

% Inverse FT both high and low pass images into (x, y)
if strcmp(in_dom, 'xy') || strcmp(in_dom, 'xky')        % dat in (x, Ky)
    lowPass = ifftc(lowPass, PE_DIM);
    dat = ifftc(dat, PE_DIM);
else if strcmp(in_dom, 'kxy') || strcmp(in_dom, 'kxky') % dat in (Kx, Ky)
        lowPass = ifft2c(lowPass);
        dat = ifft2c(dat);
    end
end

% Apply Phase Correction to high pass images - use real part only
scalar = 0.000001;
phaseCorrection = abs(lowPass) ./ (lowPass + scalar);
lpImagePhase = (lowPass + scalar) ./ abs(lowPass);
dat = dat .* phaseCorrection;
dat = real(dat); % real part

for iter = 1 : niter
    % Reinsert phase
    dat = dat .* lpImagePhase;
    
    % Transform high-pass image back to (x, Ky)
    dat = fftc(dat, PE_DIM);
    
    % Here the output HP image is as follows:
    % Orig*Merge + HP*(1-Merge) = (Orig-Hp)*Merge + HP;
    for i1 = 1 : datsz(1)
        for i3 = 1 : datsz(3)
            for i4 = 1 : datsz(4)
                for i5 = 1 : datsz(5)
                    for i6 = 1 : datsz(6)
                        for i7 = 1 : datsz(7)
                            dat(i1, :, i3, i4, i5, i6, i7) = (dat_input(i1, :, i3, i4, i5, i6, i7) - dat(i1, :, i3, i4, i5, i6, i7)) .* mergeWindow + dat(i1, :, i3, i4, i5, i6, i7);
                        end
                    end
                end
            end
        end
    end
    
    % Transform high-pass image back to (x,y)
    dat = ifftc(dat, PE_DIM);
    
    % Correct phase
    dat = dat .* phaseCorrection;
    dat = real(dat);
end

% Load output
switch out_dom
    case 'xy'
        % Do nothing
    case 'xky'
        dat = fftc(dat, PE_DIM);
    case 'kxy'
        dat = fftc(dat, FE_DIM);
    case 'kxky'
        dat = fft2c(dat);
end

return

function [lowPassFilter, highPassFilter, mergeWindow] = generate_filters(ntran, nacq, npres)
%
% function [lowPassFilter, highPassFilter, mergeWindow] = generate_filters(ntran, nacq, npres)
%
% Inputs
%   ntran          - Number of points in transition band.
%   nacq           - Number of partial-k acquired lines.
%   npres          - Number of total prescribed lines.
%
% Outputs
%   lowPassFilter  - Low-pass filter.
%   highPassFilter - High-pass filter.
%   mergeWindow    - Window for data merging.
%
% Generate filters for Homodyne detection reconstruction.
%

symmetricSamples = nacq - npres/2;
transitionWidth = ntran;
mergeWindowSize = npres;
filterSize = nacq;
dataCenter = filterSize - symmetricSamples;
dataOffset = mergeWindowSize/2 - dataCenter;

windowCenter2 = filterSize - (2*transitionWidth);
windowCenter1 = dataCenter - (windowCenter2 - dataCenter);

exponential1 = generate_exponential(filterSize, windowCenter1, transitionWidth);
exponential2 = generate_exponential(filterSize, windowCenter2, transitionWidth);

lowPassFilter = exponential2 - exponential1;
highPassFilter = exponential2 + exponential1;

mergeWindow = zeros(1, mergeWindowSize);
mergeWindow((dataOffset : dataOffset + length(exponential2) -1) + 1) = exponential2;

return

function vector = generate_exponential(vecSize,winWidth,transWidth)
%
% function vector = generate_exponential(vecSize,winWidth,transWidth)
%
% Generate exponentials for Homodyne detection reconstruction.
%

lowThreshold = -87.0;
highThreshold = 87.0;

vector = 1:vecSize;
vector = vector - winWidth;
vector = vector ./ transWidth;

vector(vector>highThreshold) = highThreshold;
vector(vector<lowThreshold) = lowThreshold;

vector = exp(vector) + 1.0;
vector =  1.0./vector;

return