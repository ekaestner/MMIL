function [pha_coe, uid] = rawload_ref(ny, nsl, nc, frames, slices, coils, fname)
%
% function [pha_coe, uid] = rawload_ref(ny, nsl, nc[, frames, slices, coils, fname])
%
% Load the GE ref.dat file for phase correction in EPI.
%
% Inputs
%   ny     - Total # of ky lines in the acquisition corresponding to the
%            ref.dat file to be read.
%   nsl    - Total # of slices in the acquisition.
%   nc     - Total # of coils in  the acquisition.
%   frames - Frames we want the phase coefficients of. Default: all frames.
%   slices - Slices we want the phase coefficients of. Default: all slices.
%   coils  - Coils we want the phase coefficients of. Default: all coils.
%   fname  - File name. Default: 'ref.dat'.
%
% Output
%   pha_coe - Coefficients for x-ky phase correction. Dim: [2(0th order, 1st order), ny, nsl, nc].
%   uid     - the uinique id from the corresponding p-file (if available)
%             (will be empty if the file does not have a uid header)
%
% (c) Kangrong Zhu,     Stanford University     July 2012

% -- Parse inputs.
if ~exist('frames', 'var');      frames = [];            end;
if ~exist('slices', 'var');      slices = [];            end;
if ~exist('coils', 'var');       coils  = [];            end;
if ~exist('fname', 'var');       fname = 'ref.dat';      end;

if (isempty(frames));            frames = 1:ny;          end;
if (isempty(slices));            slices = 1:nsl;         end;
if (isempty(coils));             coils  = 1:nc;          end;

% -- Specify some constants. (Guessed from the data loaded from the ref.dat file.)
MAX_NUM_FRAMES = 512; % Maximum # of frames which could be stored in the ref.dat file.
                      % Only the first ny points in a 512 point block are non-zero.
PHA_ORDER = 1;        % Order for the phase term.

% -- Open file.
fid = fopen(fname,'r','l');
if fid == -1
    error('Can not open file %s.\n',fname);
end

% Attempt to read the 32-byte UID header
uid = fread(fid, 32, 'unsigned char')';
% Check for the magic 7 bytes that the 32-byte UIDs always begin with
if ~all(uid(1:7) == [43 59 149 27 34 71 42])
    uid = [];
    frewind(fid);
end

% -- Read data
nexpected = MAX_NUM_FRAMES * (PHA_ORDER+1) * nc * nsl; % Expected # of points
[pha_coe, nread] = fread(fid, nexpected, 'float32');
fclose(fid);
if nread == nexpected
    pha_coe = reshape(pha_coe, MAX_NUM_FRAMES, PHA_ORDER+1, nc, nsl);
    pha_coe = permute(pha_coe, [2, 1, 4, 3]);
    pha_coe = pha_coe( :, frames, slices, coils );
else
    error('%d of %d expected data points are read.', nread, nexpected);
end

