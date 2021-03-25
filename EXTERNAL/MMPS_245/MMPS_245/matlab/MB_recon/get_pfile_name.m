function [pfile, ref, vrgf, refp, noise_fname, pfile_name] = get_pfile_name(d)
% function [pfile, ref, vrgf, refp, noise_fname, pfile_name] = get_pfile_name(d)
%
% Get the full paths to the pfile and the associated ref.dat, vrgf.dat and reference scan pfile.
%
% Input
%   d          - Filename of the pfile, or the directory containing a pfile.
%
% Outputs
%   pfile      - Filename of the p-file.
%   ref        - Filename of the ref.dat file associated with the p-file.
%   vrgf       - Filename of the vrgf.dat file associated with the p-file.
%   refp       - Filename of reference scan pfile associated with the p-file.
%   noise_fname- Filename of the noise.dat (noise data) file associated with the p-file.
%   pfile_name - Name of pfile, PXXXXX.7.
%
% (c) Kangrong Zhu      Stanford University     Oct 2013

if ~exist('d', 'var') || isempty(d)
    error('No pfile specified!');
end

% P-file
if isdir(d)
    pfile_dir = d;
    pfiles_in_dir = dir(fullfile(pfile_dir, 'P*.7'));
    for file_idx = 1 : length(pfiles_in_dir)
        if length(pfiles_in_dir(file_idx).name) == 8 % PXXXXX.7
            pfile_name = pfiles_in_dir(file_idx).name;
        end
    end
    pfile = fullfile(pfile_dir, pfile_name);
else
    pfile = d;
    [pfile_dir, name, ext] = fileparts(pfile);
    pfile_name = [name ext];
end

% ref.dat file
ref = get_pfile_associated_file(pfile_dir, pfile_name, 'ref.dat');

% vrgf.dat file
vrgf = get_pfile_associated_file(pfile_dir, pfile_name, 'vrgf.dat');

% Reference scan pfile
refp = get_pfile_associated_file(pfile_dir, pfile_name, 'refscan.7');

% noise.dat file (Noise data)
noise_fname = get_pfile_associated_file(pfile_dir, pfile_name, 'noise.dat');

return

function file = get_pfile_associated_file(pfile_dir, pfile_name, associated_file_type)
%
% Get filename of the ref.dat file or vrgf.dat file or reference scan pfile associated with a pfile.
%
% Inputs:
%   pfile_dir            - Directory of the p-file.
%   pfile_name           - P-file name, PXXXXX.7.
%   associated_file_type - String 'ref.dat' or 'vrgf.dat' or 'refscan.7' or 'noise.dat'.
%
% Output
%    file                - File name of the associated file.
%
% (c) Kangrong Zhu      Stanford University     Oct 2013

file = fullfile(pfile_dir, [pfile_name '_' associated_file_type]);
if ~exist(file, 'file')
    file = fullfile(pfile_dir, ['_' pfile_name '_' associated_file_type]);
end
if ~exist(file, 'file')
    file = fullfile(pfile_dir, associated_file_type);
end
if ~exist(file, 'file')
    files_in_dir = dir(fullfile(pfile_dir, ['*' associated_file_type]));
    if ~isempty(files_in_dir)
        file = files_in_dir(1).name;
        file = fullfile(pfile_dir, file);
    end
end
if ~exist(file, 'file')
    file = '';
end

return
