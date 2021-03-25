function foldername = get_foldername(fpath)
%
% function foldername = get_foldername(fpath)
%
% Get the folder name from the path to a file or the path to a folder.
% Example
%   fn = get_foldername('~\path\to\data\PXXXXX.7') => fn = 'data';
%   fn = get_foldername('~\path\to\data\textfile') => fn = 'data';
%   fn = get_foldername('~\path\to\data\')         => fn = 'data';
%   fn = get_foldername('~\path\to\data')          => fn = 'to';
%
% Input
%   fpath      - String for the path to a file or the path to a folder.
%
% Output
%   foldername - If 'fpath' is the path to a file, this is the name of
%                the folder containing the file; If 'fpath' is the path to
%                a folder, this is the name of that folder.
%
% (c) Kangrong Zhu,     Stanford University     June 2013

pathstr = fileparts(fpath);

for idx = length(pathstr) : -1 : 1
    if strcmp(pathstr(idx), '\') || strcmp(pathstr(idx), '/')
        foldername = pathstr(idx+1 : end);
        break;
    end
end

if ~exist('foldername', 'var')
    foldername = pathstr;
end