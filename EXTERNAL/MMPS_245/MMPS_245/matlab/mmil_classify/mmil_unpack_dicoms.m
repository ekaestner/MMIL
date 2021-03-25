function [SeriesInfo,errcode,msg] = mmil_unpack_dicoms(indir,varargin)
%function [SeriesInfo,errcode,msg] = mmil_copy_dicoms(indir,[options])
%
% Required Input:
%   indir: full path of input directory containing dicoms
%     recursive search will be used to find dicoms in subdirectories
%
% Optional Parameters:
%   'outdir': output directory; subdirectories will be created for each series
%     {default = pwd}
%   'linkflag': [0|1] create symbolic links instead of copying
%     {default = 1}
%   'batch_limit': maximum number of lines per unix call
%     {default = 250}
%   'forceflag': [0|1] overwrite existing output files
%     {default = 0}
%
% Output:
%   errcode: returns 1 if error, 0 if successful
%   msg: error message
%
% Created:  04/15/11 by Don Hagler
% Last Mod: 05/31/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'outdir',pwd,[],...
  'linkflag',true,[false true],...
  'batch_limit',250,[],...
  'forceflag',false,[false true],...
});
SeriesInfo = [];
errcode = 0;
msg = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);
fname_info = sprintf('%s/OrigSeriesInfo.mat',parms.outdir);
if ~exist(fname_info,'file') || parms.forceflag
  fprintf('%s: sorting dicoms...\n',mfilename);
  tic;
  [SeriesInfo,errcode,msg] = mmil_sort_dicoms(indir,...
    'batch_limit',parms.batch_limit);
  toc;
  if errcode, return; end;
  save(fname_info,'SeriesInfo');
else
  load(fname_info);
end;

fname_info = sprintf('%s/SeriesInfo.mat',parms.outdir);
if ~exist(fname_info,'file') || parms.forceflag
  fprintf('%s: copying dicoms...\n',mfilename);
  tic;
  [SeriesInfo,errcode,msg] = mmil_copy_dicoms(SeriesInfo,'outdir',parms.outdir,...
    'batch_limit',parms.batch_limit,...
    'linkflag',parms.linkflag,'forceflag',parms.forceflag);
  toc;
  if errcode, return; end;
  save(fname_info,'SeriesInfo');
else
  load(fname_info);
end;
