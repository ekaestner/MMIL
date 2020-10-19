% An interface for offline reconstruction of slice-multiplexed EPI data.
%
% Use:
%  1. In this script, manually set values in the 'Set these parameters manually' section.
%  2. (TODO: save the parameters into the p-file header correctly)
%     In mux_epi_params.m, manually set value for: ychop(Usually this should be true)
%  3. Pick the p-file for the actual mux scan when asked so.
%  4. (Skipped if using internal calibration) Pick the p-file for the external calibration scan when asked so.
%
% Results:
%  1. The reconstructed images are saved in matrix 'd'.
%     Dimension [FE, PE, Slice(=total # of slices), temporalPhases, Coil].
%
% (c) Kangrong Zhu, Stanford University     Aug 2012

clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%% Start: Set these parameters manually %%%%%%%%%%%%%%%%%%
slices      = 1;      % Muxed slices to reconstruct. Don't pass this parameter or pass an empty matrix to reconstruct all slices.
nt_to_recon = [];     % Number of time points to reconstruct, excluding the first few mux phase cycling time points. Don't pass or pass an empty matrix to reconstruct all time points.
n_vcoils    = 16;     % Number of virtual coils in coil compression. Pass an empty matrix for no coil compression.
debug       = 1;      % True: Calculate intermediate results, print out messages.
recon_method= '1Dgrappa'; % '1Dgrappa' or 'sense' or 'slice-grappa' or 'split-slice-grappa'.
apply_fermi = 0;      % True: Apply Fermi filter.
use_homodyne = 1;     % True: Use homodyne for partial ky acquisition; False: use zero-filling.
notch_thresh = 0;     % If non-zero, a denotching filter will be applied. Any point in k-space with a value lower than notch_thresh will be replaced by an adjacent time point.
save_res    = 0;      % True: Save results to a file.

if exist('/home/', 'dir') % For Linux
    mux_epi_recon_dir = '~/Dropbox/mux_epi_recon_current/mux_epi_recon/';
else                      % For Windows
    mux_epi_recon_dir = './';
end
%%%%%%%%%%%%%%%%%%%  End: Set these parameters manually %%%%%%%%%%%%%%%%%%

addpath(mux_epi_recon_dir, [mux_epi_recon_dir 'pfile_reading'], [mux_epi_recon_dir 'disp_and_sim_tools'], '-begin');
addpath([mux_epi_recon_dir 'utils'], '-end');

% Select the mux scan.
[fname_mux, dir_mux, idx_mux] = uigetfile('P*.7', 'Pick an Accelerated P-file');
if idx_mux == 0;        error('Must select the mux scan.');     end % User pressed Cancel.

% Check whether internal or external calibration is needed.
[mux_header, rhuser] = rawheadX([dir_mux fname_mux]);
internal_cal = rhuser(6);
num_mux_cycle = rhuser(8);
internal_cal = internal_cal && (num_mux_cycle > 0);

% Select the calibration scan when external calibration is needed.
if internal_cal
    dir_cal = [];
else 
    [tmp, dir_cal, idx_cal] = uigetfile('P*.7', 'Pick a Calibration P-file');
    if idx_cal == 0
        error('Must select the calibration scan because external calibration is needed.');
    end
end

% The actual recon routine.
if save_res
    time = get_current_time();
    outfile = ['res_', get_foldername(dir_mux), '_started_', time{1}, time{2}, time{3}, '_', time{4}, time{5}, time{6}, '.mat'];
else
    outfile = [];
end

if debug;       tic;    end

[d, snr_d] = mux_epi_main([dir_mux, fname_mux], outfile, dir_cal, slices, nt_to_recon, n_vcoils, debug, recon_method, apply_fermi, use_homodyne, notch_thresh);

if debug
    time_used = toc;

    fprintf('dir_mux \t %s \ndir_cal \t %s \n', dir_mux, dir_cal);
    fprintf('Time used = %d s.\n', time_used);
end

beep;