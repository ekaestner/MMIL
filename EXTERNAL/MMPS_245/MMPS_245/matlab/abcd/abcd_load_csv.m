function abcd_info = abcd_load_csv(fname_csv,outdir)
%function abcd_info = abcd_load_csv(fname_csv,[outdir])
%
% Purpose: load info from csv, write to mat file for future access
%
% Required Input:
%   fname_csv: name of csv info file
%
% Optional Input:
%   outdir: cache output path
%     if  relative path, relative to path of fname_csv
%    {default = 'cache'}
%
% Created:  01/06/16 by Don Hagler
% Last Mod: 02/06/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('outdir','var') || isempty(outdir)
  outdir = 'cache';
end;

if ~exist(fname_csv,'file')
  error('file %s not found',fname_csv);
end;

[fpath,fstem] = fileparts(fname_csv);
if mmil_isrelative(outdir)
  outdir = [fpath '/' outdir];
end;
mmil_mkdir(outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_cache = sprintf('%s/%s_%s.mat',...
  outdir,fstem,datestr(now,'yyyymmdd'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load info
if exist(fname_cache)
  fprintf('%s: loading info from cache file %s...\n',...
    mfilename,fname_cache);
  tic;
  load(fname_cache);
  toc;
else
  mmil_mkdir(outdir);
  fprintf('%s: loading info from csv file %s...\n',...
    mfilename,fname_csv);
  tic;
  abcd_info = mmil_csv2struct(fname_csv);
  toc;
  fprintf('%s: saving info to cache file %s...\n',...
    mfilename,fname_cache);
  tic;
  save(fname_cache,'abcd_info');
  toc;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

