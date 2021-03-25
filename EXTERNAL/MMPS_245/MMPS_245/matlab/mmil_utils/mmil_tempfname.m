function [tempfname] = mmil_tempfname(prefix,outdir)
% function [tempfname] = mmil_tempfname(prefix,outdir)
%
% Optional Input:
%   prefix: prefix of new temp filename
%   outdir: preferred output location
%
% Output:
%   tempfname: Unique temp filename
% Created:                09/05/17 by Feng Xue
% Prev Mod:               10/20/17 by Don Hagler
% Last Mod:               11/2/17 by Feng Xue
%

  [~,tmp_stem,~] = fileparts(tempname);

%  rand('twister',sum(100*clock));
%% NOTE: twister is deprecated in later versions
%%       but RandStream not found for 2007
  stream = RandStream.create('mt19937ar','Seed',sum(100*clock));
%  RandStream.setDefaultStream(stream);
  RandStream.setGlobalStream(stream);
  if exist('prefix','var') && ~isempty(prefix)
    tmp_stem = [prefix '_' tmp_stem '_' num2str(round(rand(1)*10000000))];
  else
    tmp_stem = [tmp_stem '_' num2str(round(rand(1)*10000000))];
  end

  if ~exist('outdir','var') || isempty(outdir)
    tmpdir = '/scratch';
  else
    tmpdir = outdir;
  end

  % test if we can write to output directory
  tmpfile = [tmpdir '/' tmp_stem]; 
  cmd = sprintf('touch %s',tmpfile);
  [status,result] = unix(cmd);
  if ~status % can write to this directory
    if exist(tmpfile,'file'), warning('off','all');delete(tmpfile);warning('on','all'); end;
  else
    tmpdir = '/tmp';
%    tmpdir = '.';
  end;
  if ~exist(tmpdir,'dir')
    error('tmpdir %s not found',tmpdir);
  end;

  tempfname = regexprep(sprintf('%s/%s',tmpdir,tmp_stem),'//','/');

return;
