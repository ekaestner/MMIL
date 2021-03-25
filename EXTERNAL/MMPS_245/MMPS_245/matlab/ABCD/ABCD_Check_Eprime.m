function fname_out = ABCD_Check_Eprime(ProjID,varargin)
%function fname_out = ABCD_Check_Eprime(ProjID,[options])
%
% Purpose: identify eprime files for each task fMRI series
%
% Optional Input: ('key', value pairs)
%   'rootdir': root input directory containing container type  dirs
%     {default = '/space/syn05/1/data/MMILDB'}
%   'indir': root input directory containing eprime files
%     if empty, will be in {rootdir}/aux_incoming
%     {default = []}
%   'fname_info': spreadsheet containing PC info for each series
%     if empty, will use pcinfo from MetaData
%     {default = []}
%   'infix': string attached to input pcinfo (if fname_info is empty)
%     {default = 'pcinfo'}
%   'tasknames': cell array of task names
%     {default = {'MID','SST','nBack'}}
%   'outdir': output directory
%     if not  specified, output in MetaData/{ProjID}
%     {default = []}
%   'outstem': output file stem
%     if not  specified, {ProjID}
%     {default = []}
%   'outfix': string attached to output file names
%     {default = 'eprime'}
%   'scanner_flag': [0|1] find eprime files only for tasks done in scanner
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% Created:  11/30/16 by Don Hagler
% Prev Mod: 07/19/17 by Feng Xue
% Last Mod: 08/02/17 by Don Hagler,  Octavio Ruiz (2017aug-oct02, 16)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'rootdir','/space/syn05/1/data/MMILDB',[],...
  'indir',[],[],...
  'fname_info',[],[],...
  'tasknames',{'MID','SST','nBack'},{'MID','SST','nBack'},...
  'outdir',[],[],...
  'outstem',[],[],...
  'infix','pcinfo',[],...
  'outfix','eprime',[],...
  'scanner_flag',false,[false true],...
  'forceflag',false,[false true],...
  'sers_eprime_time_diff_max',[],[], ...   % Octavio Ruiz (2017oct02)
...
  'tags',{'indir','fname_info','tasknames','outdir',...
          'outstem','outfix','scanner_flag','forceflag', ...
          'sers_eprime_time_diff_max'   % Octavio Ruiz (2017oct02)
          },[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

homedir = getenv('HOME');
parms.indir = sprintf('%s/%s/aux_incoming',parms.rootdir,ProjID);
if ~exist(parms.indir,'dir')
  error('input directory %s not found',parms.indir);
end;
if isempty(parms.fname_info)
  parms.fname_info = sprintf('%s/MetaData/%s/%s_%s.csv',...
    homedir,ProjID,ProjID,parms.infix);
end;
if isempty(parms.outdir)
  parms.outdir = sprintf('%s/MetaData/%s',homedir,ProjID);
end;
if isempty(parms.outstem)
  parms.outstem = ProjID;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check eprime files
args = mmil_parms2args(parms,parms.tags);
fname_out = abcd_check_eprime_all(args{:});
