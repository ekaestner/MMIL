function [pha_coe, descend_acq_ref, dat] = rawload_ref_pfile(frames, slices, coils, refpfile, pha_coe_1stpass, ramp_flt, pccoil, do_quad_final, debug)
%
% function [pha_coe, descend_acq_ref, dat] = rawload_ref_pfile([frames, slices, coils, refpfile, pha_coe_1stpass, ramp_flt, pccoil, do_quad_final, debug])
%
% Calls rawload_ref_data.m to calculate EPI x-ky phase correction coefficients from reference scan pfile.
%
% Inputs
%   frames          - Frames we want the phase coefficients of. Default: All frames.
%   slices          - Slices(normally ordered indices, not pfile indices) we want the phase coefficients of. Default: All slices.
%   coils           - Coils we want the phase coefficients of. Default: All coils.
%   refpfile        - Reference scan pfile. Default: 'refpfile.7'
%   pha_coe_1stpass - First-pass phase correction coefficients. See rawload_ref_data.m for details. Default: [].
%   ramp_flt        - First-pass ramp sampling correction filter. See rawload_ref_data.m for details. Default: [].
%   pccoil          - Whether to average the phase correction coefficients across coil. See rawload_ref_data.m for details. Default: Set to rdb_hdr.pccoil in pfile header.
%   do_quad_final   - True(Default): Conduct quadsmooth along ky direction after the x-ky phase correction coefficients are calculated by least squares fitting.
%   debug           - True: Display per slice per coil static phase correction coefficients and print messages. Default: False.
%
% Outputs
%   pha_coe         - EPI x-ky phase correction coefficients. Dim: [2(0thOrderPhase, 1stOrderPhase), Ky, Slice, Coil]
%   descend_acq_ref - True: the slices in the reference scan were acquired in descending order(from higher frequency to lower frequency in the slice select direction); False: acquired in ascending order.
%   dat             - k-space data loaded from the reference scan pfile. Dim: [Kx, Ky, Echo(=1), Slice, Coil]
%
% (c) Kangrong Zhu      Stanford University     Sep 2013

%% Set defaults
if ~exist('frames', 'var');          frames = [];          end;
if ~exist('slices', 'var');          slices = [];          end;
if ~exist('coils', 'var');           coils  = [];          end;
if ~exist('pha_coe_1stpass', 'var'); pha_coe_1stpass = []; end;
if ~exist('ramp_flt', 'var');        ramp_flt = [];        end;
if ~exist('do_quad_final', 'var');   do_quad_final = true; end;

if ~exist('refpfile', 'var') || isempty(refpfile)
    refpfile = 'refpfile.7';
end
if ~exist('debug', 'var') || isempty(debug)
    debug = false;
end

%% Parameter structure, similar to mux_epi_params.m, but only including needed fields. 
p.FE_DIM = 1;
p.PE_DIM = 2;
p.EC_DIM = 3;
p.SL_DIM = 4;
p.C_DIM  = 5;
p.T_DIM  = 6;
p.KZ_DIM = 7;
p.KEEP_ORIG_SZ        = 1;
p.NUM_COE_PER_KY_LINE = 2;

% For ky traversal direction
p.TOP_DOWN   = 0;
p.CENTER_OUT = 1;
p.BOTTOM_UP  = 2;

hdr = read_MR_headers(refpfile, 'all', 'raw');                  % Full header

if ~exist('pccoil', 'var') || isempty(pccoil)
    pccoil = hdr.rdb_hdr.pccoil;                                % Coil averaging when calculating x-ky phase correction coefficients. 0: No averaging across coil; >=1 && <=ncoils: Use one of the coils' coefficients for all coils; -1: average across coils.
end

p.pcslice      = hdr.rdb_hdr.pcspacial;                         % Slice averaging when calculating x-ky phase correction coefficients. 0: No averaging across slice; >=1 && <=nslics: Use one of the slices' coefficients for all slices; -1: average across slices.
p.nshot        = hdr.rdb_hdr.ileaves;                           % Number of shots (equals number of interleaves)
p.kydir        = hdr.rdb_hdr.kydir;                             % 0=top down; 1=center out; 2=bottom up; 3=seq center out; 4=rev center out; 5 = quasi centric
p.num_slices   = hdr.image.slquant;                             % Number of slices
p.debug        = debug;
p.frames       = frames;
p.coils        = coils;

switch p.num_slices
    case hdr.rdb_hdr.nslices/hdr.rdb_hdr.npasses
        p.acq_order = 'interleaved';                                       % In this case, hdr.data_acq_tab.pass_number is always 0.
        p.sl_acq_order = hdr.data_acq_tab.slice_in_pass(1 : p.num_slices); % Acquisition order of prescribed slices.
    case hdr.rdb_hdr.nslices
        p.acq_order = 'sequential';                                        % In this case, hdr.data_acq_tab.slice_in_pass is always 1.
        p.sl_acq_order = hdr.data_acq_tab.pass_number(1 : p.num_slices) + 1;
    otherwise
        error('Inconsistent slice information in header.');
end

if hdr.series.start_loc > hdr.series.end_loc
    p.descend_acq = true;
else
    p.descend_acq = false;
end
descend_acq_ref = p.descend_acq;                                % for output

% For reference scan data saved by RDS
p.rds_data_order = ( (hdr.rdb_hdr.recon>=9000) || (hdr.rdb_hdr.user16==1));

% Check values
if (p.pcslice ~= 0) && debug
    fprintf('p.pcslice = %.2f\n', p.pcslice);
end
if (pccoil ~= 0) && debug
    fprintf('pccoil = %.2f\n', pccoil);
end

%% Load k-space data
if isempty(slices)
    slices = 1 : p.num_slices;
end
slices_in_pfile = p.sl_acq_order(slices);                       % Pfile slice indices

if ~p.rds_data_order                                            % Pfile saved by the scanner
    dat = rawloadX(refpfile, [], [], slices_in_pfile, p.coils); % Dim: [Kx, Ky, Echo(=1), Slice, Coil]. Load only the desired slices and coils because each slice and each coil will be calculated individually.
else                                                            % Pfile saved by RDS
    p.inplane_R = hdr.rdb_hdr.ileaves;
    p.kissoff_views = hdr.rdb_hdr.kissoff_views;
    p.num_passes = hdr.rdb_hdr.npasses;
    p.mux = hdr.rdb_hdr.user6;
    p.num_mux_cycle = hdr.rdb_hdr.user7;
    p.num_tphases_to_recon = -(p.mux*p.num_mux_cycle-1);
    p.tphases_to_load = 1;
    p.slice_shuffling = 0;
    p.pfile_header.version = hdr.rdb_hdr.rdbm_rev;
    p.pfile_header.npasses = p.num_passes ;
    p.pfile_header.nslices = hdr.rdb_hdr.nslices;
    p.pfile_header.nframes = hdr.rdb_hdr.nframes;
    p.pfile_header.nechoes = hdr.rdb_hdr.nechoes;
    p.pfile_header.hnover  = hdr.rdb_hdr.hnover;
    p.pfile_header.frsize  = hdr.rdb_hdr.frame_size;
    p.pfile_header.ncoils  = hdr.rdb_hdr.dab(2)-hdr.rdb_hdr.dab(1)+1;
    p.pfile_header.ptsize  = hdr.rdb_hdr.point_size;
    p.pfile_header.rawhdrsize  = hdr.rdb_hdr.off_data;
    
    dat = load_raw_tseries(refpfile, p, slices_in_pfile, p.tphases_to_load );
end

pha_coe = rawload_ref_data(dat, p, pha_coe_1stpass, ramp_flt, pccoil, do_quad_final, false);

return