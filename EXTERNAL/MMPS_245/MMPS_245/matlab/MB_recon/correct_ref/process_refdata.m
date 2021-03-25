function pha_coe = process_refdata(refData, p)
%
% function pha_coe = process_refdata(refData, p)
%
% Calculates EPI x-ky phase correction coefficients from reference scan k-space data.
%
% Inputs
%   refData - Reference scan k-space data. Dim: [Ky(=ny), Kx(=nx), Slice(=nsl), Coil(=nc)].
%   p       - Parameter structure with fields:    
%             debug         - True: Display per slice per coil static phase correction coefficients. Default: False.
%             pccoil        - 0(Default):      No averaging across coil; 
%                             >=1 && <=ncoils: Use one of the coils' coefficients for all coils; 
%                             -1:              Average across coils.
%             pcslice       - 0(Default):      No averaging across slice;
%                             >=1 && <=nslics: Use one of the slices' coefficients for all slices;
%                             -1:              Average across slices.
%             nshot         - Number of shots.
%             do_quad_final - True: Smooth the least squares fitted coefficients with a quadratic smoothing function. Default: True.
%
% Output
%   pha_coe - EPI x-ky phase correction coefficients. Dim: [2(0thOrderPhase, 1stOrderPhase), Ky(=ny), Slice(=nsl), Coil(=nc)]
% 
% Original code                                                                 Suchandrima Banerjee, GE Healthcare     Aug 2013
% Modified According to the rc_pc_xyfit.cpp function from GE Healthcare         Kangrong Zhu, Stanford University       Jan 2014

%% Parameters
if ~exist('p', 'var') || isempty(p)
    p = struct();
end

if ~isfield(p, 'debug') || isempty(p.debug)
    p.debug = false;
end

if ~isfield(p, 'pccoil') || isempty(p.pccoil)
    p.pccoil = 0;
end

if ~isfield(p, 'pcslice') || isempty(p.pcslice)
    p.pcslice = 0;
end

if ~isfield(p, 'nshot') || isempty(p.nshot)
    p.nshot = 1;
end

if ~isfield(p, 'do_quad_final') || isempty(p.do_quad_final)
    p.do_quad_final = true;
end

NUM_COE_PER_KY_LINE = 2;
DIV_BY              = 10^(-8);
SCALE_FACTOR        = 0.25;

[ny_tot, nx, nslices, ncoils] = size(refData);
ny = ny_tot / p.nshot;
if round(ny) ~= ny
    error('Number of phase encoding lines per shot is not an integer.');
end

%% Calculate x-ky phase correction coefficients
weightAcrossSliceCoil = zeros(nslices, ncoils);
linearCoefAll = zeros(ny, nslices, ncoils);
constCoefAll = zeros(ny, nslices, ncoils);
for coil = 1:ncoils
    for slice = 1:nslices
        % Ref data of the current slice and coil
        curRefData = refData(:, :, slice, coil);                     % Dim: [Ky(=ny), Kx(=nx)]
        
        % Bandpass assymetry correction
        % Not done in current implementation

        % Row Fourier Transform
        curRefData = ifftc(curRefData, 2);                           % Dim: [Ky(=ny), X(=nx)]

        % Weights used to average the resulting coefs between slices and coils
        weightAcrossSliceCoil(slice, coil) = sum(abs(curRefData(:)));
        
        % For multishot ref scan, averge along the 'shot' axis
        if p.nshot > 1
            tmp = zeros(ny, nx, p.nshot);
            for shot = 1 : p.nshot
                tmp(:, :, shot) = curRefData(shot : p.nshot : end, :);
            end
            curRefData = mean(tmp, 3);
        end
        
        % Calculate the phase difference between each two adjacent frames
        phaseDiffData = curRefData(1:(end-1), :) .* conj(curRefData(2:end, :));
        phaseDiffData(2:2:end, :) = conj(phaseDiffData(2:2:end, :)); % This will give sign alteration to create 'uniform' phase in kx-ky
        phaseDiff = atan2(imag(phaseDiffData), real(phaseDiffData)); % This effectively does phase(line+1) - phase(line) and also maps to [-pi, pi]
        
        % Magnitude of the phase difference data
        ap_ref_mag = abs(phaseDiffData);
        
        % Remove linear phase estimated by Ahn's method
        tmp = phaseDiffData(:, 1:end-1) .* conj(phaseDiffData(:, 2:end));
        tmp = sum(tmp(:));
        linearPhaseParaAhn = - atan2(imag(tmp), real(tmp));
        
        linearPhaseParaAhn = repmat(linearPhaseParaAhn, [ny-1, 1]);
        fit_x = -nx/2 : (nx/2-1);
        phaseDiffLinearTermRemoved = phaseDiff - linearPhaseParaAhn * fit_x;
        
        % Get the angle which is the constant phase
        tmp = ap_ref_mag .* exp(1i * phaseDiffLinearTermRemoved);
        tmp = sum(tmp(:));
        ap_const_term = atan2(imag(tmp(:)), real(tmp(:)));           % Map the remaining phase to [-pi, pi]
        
        % Subtract const term from the phase
        ap_ref_phase = phaseDiffLinearTermRemoved - ap_const_term;   % Phase of the data with the linear and constant phase components removed
        tmp = exp(1i * ap_ref_phase);
        ap_ref_phase = atan2(imag(tmp), real(tmp));                  % Unwrapped phase of the data with the linear and constant phase components removed
        
        % Alternative to the above atan2, use quadsmooth to smooth and interpolate across low signal regions
        % Not done in current implementation
        
        % Scale down the magnitude vector
        ap_ref_mag = ap_ref_mag .* DIV_BY;
        
        % Add the linear and const terms back to ref_phase
        ap_ref_phase = ap_ref_phase + linearPhaseParaAhn*fit_x + ap_const_term;
        
        % Least-squares fitting
        ap_const_term = zeros(ny-1, 1);
        ap_linear_term = zeros(ny-1, 1);
        for ind = 1 : ny-1
            fit_wt = ap_ref_mag(ind, :);                             % The weight vector is square of the ref magnitude vector within the 'lsqfit' function
            fit_y  = ap_ref_phase(ind, :);
            [ap_linear_term(ind), ap_const_term(ind)] = lsqfit(fit_x, fit_y, fit_wt);
        end
        
        % Interpolate constant and linear coefficients to ny points.
        ap_bo_coef = SCALE_FACTOR .* (ap_const_term(2:end) + ap_const_term(1:end-1));
        ap_lin_coef = SCALE_FACTOR .* (ap_linear_term(2:end) + ap_linear_term(1:end-1));
        if ny > 2
            ap_bo_coef = [ap_bo_coef(1); ap_bo_coef; ap_bo_coef(end)];
            ap_lin_coef = [ap_lin_coef(1); ap_lin_coef; ap_lin_coef(end)];
        else                                                         % If ny == 2, lsqfit() is only called once
            ap_bo_coef = [ap_const_term; ap_const_term] * 0.5;
            ap_lin_coef = [ap_linear_term; ap_linear_term] * 0.5;
        end
        
        % Now need to smooth coefficients, after interpolate to ny points.
        if p.do_quad_final && (ny > 2)
            tmpx = (0 : 1 : ny-1).';
            tmpwt = ones(ny, 1);
            kernel = 15;
            
            ap_bo_coef = quadsmooth1(tmpx, ap_bo_coef, tmpwt, kernel);            
            ap_bo_coef = quadsmooth1(tmpx, ap_bo_coef, tmpwt, kernel);
                        
            ap_lin_coef = quadsmooth1(tmpx, ap_lin_coef, tmpwt, kernel);            
            ap_lin_coef = quadsmooth1(tmpx, ap_lin_coef, tmpwt, kernel);
        end
        
        % Negate half of the coefficients
        ap_bo_coef(2:2:end) = - ap_bo_coef(2:2:end);
        ap_lin_coef(2:2:end) = - ap_lin_coef(2:2:end);
        
        % Output coefficients
        linearCoefAll(:, slice, coil) = ap_lin_coef;
        constCoefAll(:, slice, coil) = ap_bo_coef;        
    end
end

%% Display results
% Display per slice per coil static PC coefficients
if p.debug    
    tmp = linearCoefAll;
    tmp(2:2:end,:,:) = -tmp(2:2:end,:,:);                            % display the sign alternated version, so a smooth map is expected
    figure; imshowALL(tmp);
    
    tmp = constCoefAll;
    tmp(2:2:end,:,:) = -tmp(2:2:end,:,:);
    figure; imshowALL(tmp);
end

%% Average across coils and/or slices when the corresponding flags are turned on
wt = weightAcrossSliceCoil;                                          % save the original weights because this section can do averaging on them
if p.pccoil == 0                                                     % no averaging across coil
    linearCoef = linearCoefAll;
    constCoef = constCoefAll;
elseif p.pccoil>=1 && p.pccoil<=ncoils                               % use one of the coils' coef for all coils
    linearCoef = linearCoefAll(:, :, p.pccoil);
    constCoef = constCoefAll(:, :, p.pccoil);
    wt = wt(:, p.pccoil);
elseif p.pccoil == -1                                                % average across coil
    linearCoef = sum(permute(repmat(wt,[1,1,ny]),[3,1,2]).*linearCoefAll,3)...
        ./sum(permute(repmat(wt,[1,1,ny]),[3,1,2]),3);
    constCoef = sum(permute(repmat(wt,[1,1,ny]),[3,1,2]).*constCoefAll,3)...
        ./sum(permute(repmat(wt,[1,1,ny]),[3,1,2]),3);
    wt = sum(wt, 2);
else
    error('Invalid value for p.pccoil.');
end

if p.pcslice == 0                                                    % no averaging across slice
    % do nothing in this option
elseif p.pcslice>=1 && p.pcslice<=nslices                            % use one of the slices' coef for all slices
    linearCoef = linearCoef(:, p.pcslice, :);
    constCoef = constCoef(:, p.pcslice, :);
elseif p.pcslice == -1                                               % average across slice
    linearCoef = sum(permute(repmat(wt,[1,1,ny]),[3,1,2]).*linearCoef,2)...
        ./sum(permute(repmat(wt,[1,1,ny]),[3,1,2]),2);
    constCoef = sum(permute(repmat(wt,[1,1,ny]),[3,1,2]).*constCoef,2)...
        ./sum(permute(repmat(wt,[1,1,ny]),[3,1,2]),2);
else
    error('Invalid value for p.pcslice.');
end

%% Output results
if p.pccoil == 0                                                     % No averaging across coil
    nc = ncoils;
else
    nc = 1;                                                          % Use one of the coils' coef for all coils OR average across coil
end
if p.pcslice == 0                                                    % No averaging across slice
    nsl = nslices;
else                                                                 % Use one of the slices' coef for all slices OR average across slice
    nsl = 1;
end
pha_coe = zeros(ny, nsl, nc, NUM_COE_PER_KY_LINE);
pha_coe(:, :, :, 1) = -constCoef;
pha_coe(:, :, :, 2) = -linearCoef;
pha_coe = permute(pha_coe, [4, 1, 2, 3]);                            % Dim: [2(0thOrderPhase, 1stOrderPhase), Ky(=ny), Slice(=nsl), Coil(=nc)]

if p.nshot > 1                                                       % For multi-shot reference data, the same coefficients are used for all frames in a group
    pha_coe_fullsz = zeros(NUM_COE_PER_KY_LINE, ny_tot, nsl, nc);
    for shot = 1 : p.nshot
        pha_coe_fullsz(:, shot : p.nshot : end, :, :) = pha_coe;
    end
    pha_coe = pha_coe_fullsz;
end

function [slope, intercept] = lsqfit(x, y, wt)
%
% Linear least squares fit: wt(ind) * y(ind) = intercept + slope * x(ind).
%
% Inputs
%   x         - Independent data.
%   y         - Dependent data.
%   wt        - Weighting function
%
% Outputs
%   slope     - Fit result, slope.
%   intercept - Fit result, intercept.
%
% This is the Matlab version of the lsqfit function in rc_pc_xyfit.cpp from GE Healthcare.
% Kangrong Zhu    Stanford university   Jan 2014

x = x(:);
y = y(:);
wt = wt(:);

n = length(x);
if length(y)~=n || length(wt)~=n
    error('Lengths of input vector mismatch.');
end

if n < 2
    slope = [];
    intercept = [];
    return;
end

yi = y .* wt;
sx0 = sum(wt);
sx = sum(x .* wt);
sx2 = sum(x.^2 .* wt);
sy = sum(yi);
sxy = sum(x .* yi);

r = (sx0 * sx2) - (sx * sx);
if r == 0
    r = 10^(-15);
end
intercept = ((sy * sx2) - (sxy * sx)) /r;
slope = ((sx0 * sxy) - (sx * sy)) /r;

return

function [yout, dyout] = quadsmooth1(xvals, yvals, mag, kernel)
%
% Smooths a 2D set of points to lsq fit to a quadratic equation. Smooths
% each column in 2D array.
%
% Inputs
%   xvals  - Input x values. Dim: [xres, yres].
%   yvals  - Input y values. (x, y) values are to be parameterized to y =
%            ax^2 + bx + c.  Dim: [xres, yres].
%   mag    - Magnitude for weighting. Dim: [xres, yres].
%   kernel - Kernel size for smoothing each column of the 2D array. Dim: a scalar.
%
% Outputs
%   yout   - Output y values at the input xvals locations. yout(:,ind) = a*xvals(:,ind).^2 + b*xvals(:,ind) + c,
%            where a, b and c are the fitted parameters. Dim: [xres, yres].
%   dyout  - 1st derivative of the fitted function at the input xvals locations. dyout(:,ind) = 2*a*xvals(:,ind)+b. Dim: [xres, yres].
%
% This is the Matlab version of the quadsmooth1 function in rc_pc_xyfit.cpp from GE Healthcare.
% Kangrong Zhu    Stanford University     Jan 2014

[xres, yres] = size(xvals);
if (size(yvals, 1) ~= xres) || (size(yvals, 2) ~= yres)
    error('Dimension of the input ''xvals'' and ''yvals'' mismatch.');
end
if (size(mag, 1) ~= xres) || (size(mag, 2) ~= yres)
    error('Dimension of the input ''xvals'' and ''mag'' mismatch.');
end

yout = zeros(xres, yres);
dyout = zeros(xres, yres);

sumx   = 0;
sumx2  = 0;
sumx3  = 0;
sumx4  = 0;
sumy   = 0;
sumy2  = 0;
sumxy  = 0;
sumx2y = 0;
nsum   = 0;
for y = 0 : yres-1
    for x = 0 : xres+kernel-1
        x1 = x;
        x2 = round(x - kernel);
        x3 = round(x - floor(kernel/2));
        if x1 < xres
            xdatd = xvals(x1+1, y+1);
            ydatd = yvals(x1+1, y+1);
            magd  = mag(x1+1, y+1);
            
            if magd > 0.0
                wt     = magd^2;
                sumx   = sumx  + xdatd   * wt;
                sumx2  = sumx2 + xdatd^2 * wt;
                sumx3  = sumx3 + xdatd^3 * wt;
                sumx4  = sumx4 + xdatd^4 * wt;
                sumy   = sumy  + ydatd   * wt;
                sumy2  = sumy2 + ydatd^2 * wt;
                sumxy  = sumxy + xdatd*ydatd * wt;
                sumx2y = sumx2y + xdatd^2*ydatd * wt;
                nsum   = nsum + 1.0 * wt;
            end
        end
        
        if (x2 >= 0) && (x2 < xres)
            xdatd = xvals(x2+1, y+1);
            ydatd = yvals(x2+1, y+1);
            magd  = mag(x2+1, y+1);
            
            if magd > 0.0
                wt     = magd^2;
                sumx   = sumx  - xdatd   * wt;
                sumx2  = sumx2 - xdatd^2 * wt;
                sumx3  = sumx3 - xdatd^3 * wt;
                sumx4  = sumx4 - xdatd^4 * wt;
                sumy   = sumy  - ydatd   * wt;
                sumy2  = sumy2 - ydatd^2 * wt;
                sumxy  = sumxy - xdatd*ydatd * wt;
                sumx2y = sumx2y - xdatd^2*ydatd * wt;
                nsum   = nsum - 1.0 * wt;
            end
        end
        
        if (x3 >=0) && (x3 <xres)
            xdatd = xvals(x3+1, y+1);
            
            if nsum > 0.0
                D1 = sumx*sumx2 - nsum*sumx3;
                D2 = nsum*sumx2 - sumx*sumx;
                D3 = nsum*sumx2y - sumx2*sumy;
                D4 = nsum*sumxy - sumx*sumy;
                D5 = nsum*sumx3 - sumx*sumx2;
                
                a_denominator = nsum*D2*sumx4 - D2*sumx2*sumx2 + nsum*D1*sumx3 - D1*sumx*sumx2;
                c_denominator = nsum;
                
                % Check for divide by zero
                if a_denominator == 0.0
                    a_denominator = 10^(-13);
                end
                if D2 == 0.0
                    D2 = 10^(-13);
                end
                if c_denominator == 0.0
                    c_denominator = 10^(-13);
                end
                
                a = (D2*D3 + D1*D4) / a_denominator;
                b = (D4 - a*D5) / D2;
                c = (sumy - a*sumx2 - b*sumx) / c_denominator;
                
                yout(x3+1, y+1) = a*xdatd^2 + b*xdatd + c;
                dyout(x3+1, y+1) = 2*a*xdatd + b;
            end
        end % if x3
        
    end % x
end % y

return           